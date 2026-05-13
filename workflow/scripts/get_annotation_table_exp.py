# Adapted from get_annotation_table_env.py

import yaml
import argparse
import re
import os
import pandas as pd
from Bio import SeqIO
from ete4 import NCBITaxa

class TaxonCache:
    """Handles NCBI lookups with local caching to minimize database hits."""
    def __init__(self, ttaxnames=None):
        self.ncbi = NCBITaxa()
        self.cache_dom = {}             # Domain cache (taxon_name -> domain)
        self.cache_tid = {}         # Target name cache (taxid -> target_name)
        self.name_to_id_cache = {}  # Name resolution cache (name -> taxid)
        self.id_to_name_cache = {}  # Name resolution cache (taxid -> name)
        self.ttaxnames = ttaxnames

    def get_domain(self, taxon_name):
        if taxon_name in self.cache_dom:
            return self.cache_dom[taxon_name]

        try:
            # Get taxon ID from name
            name_map = self.ncbi.get_name_translator([taxon_name])
            if not name_map:
                return "Unknown (Taxon not found)"
            
            taxid = name_map[taxon_name][0]
            lineage = self.ncbi.get_lineage(taxid)
            rank_dict = self.ncbi.get_rank(lineage)
            
            # Map IDs to names to find the domain
            id_to_name = self.ncbi.get_taxid_translator(lineage)
            
            # Iterate through lineage to find the 'domain' rank
            domain = "Unknown (Domain not in lineage)"
            for node_id in lineage:
                if rank_dict.get(node_id) == "domain":
                    domain = id_to_name[node_id]
                    break
            
            self.cache_dom[taxon_name] = domain
            return domain

        except Exception as e:
            return f"Error: {str(e)}"

    def get_targettaxname(self, query):
        """
        Retrieves the scientific name for the highest available rank 
        (phylum, superkingdom, or domain) for a given taxid.
        """
        # 1. Handle Input & Resolve Names
        if isinstance(query, str):
            # Taxid input
            if query.isdigit():
                taxid = int(query)
                if taxid in self.id_to_name_cache:
                    name = self.id_to_name_cache[taxid]
                else:
                    try:
                        name = self.ncbi.get_taxid_translator([taxid])[taxid][0]
                        self.id_to_name_cache[taxid] = name
                    except:
                        return 'No_taxon_annotation'
            # taxname input
            else:
                # Check the name resolution cache first
                name = query
                if name in self.name_to_id_cache:
                    taxid = self.name_to_id_cache[name]
                else:
                    # Resolve name to ID via NCBI
                    try:
                        taxid = self.ncbi.get_name_translator([name])[name][0]
                        self.name_to_id_cache[query] = taxid
                    except:
                        raise ValueError(f"Name {name} could not be translated to taxid")
                    
        if taxid in self.cache_tid:
            return self.cache_tid[taxid]
        try:
            if self.ttaxnames is not None:
                lineage = self.ncbi.get_lineage(taxid)
                names = self.ncbi.get_taxid_translator(lineage)
                # Check if target tax names are in lineage
                for _, tn in names.items():
                    if tn in self.ttaxnames:
                        self.cache_tid[taxid] = tn
                        return tn

                # Invert ranks to search by rank name
                ranks = self.ncbi.get_rank(lineage)
                rank_to_id = {rank: tid for tid, rank in ranks.items()}
                # Priority order for highest available rank
                for rank_name in ['phylum', 'superkingdom', 'domain']:
                    if rank_name in rank_to_id:
                        tid = rank_to_id[rank_name]
                        tn = self.ncbi.get_taxid_translator([tid])[int(tid)]
                        self.cache_tid[taxid] = tn
                        return tn
                tn = 'No_taxon_annotation'
                self.cache_tid[taxid] = tn
            else:
                tn = name
            return tn
        except Exception:
            tn = 'No_taxon_annotation'
            self.cache_tid[taxid] = tn
            return tn

def process_files(input_files, fmt, taxon_lookup, dict_exp_info, dict_gene_substrate):
    all_data = []
    
    # Regex pattern: TaxonName_GeneName.fasta
    regex = fmt.format(exp="(?P<exp>[^/]+)", gene="(?P<gene>[^/]+)")
    pattern = re.compile(regex)

    for filepath in input_files:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Warning: File {filepath} not found. Skipping.")

        match = pattern.search(filepath)
        
        if not match:
            raise ValueError(f"Warning: Filename {filepath} does not match expected format. Skipping.")

        exp = match.group('exp')
        gene = match.group('gene')
        taxname = dict_exp_info[exp]['taxname']
        substrate = dict_gene_substrate[gene]
        taxon_name = taxon_lookup.get_targettaxname(taxname)
        domain = taxon_lookup.get_domain(taxon_name)

        # Parse FASTA sequences
        try:
            for record in SeqIO.parse(filepath, "fasta"):
                all_data.append({
                    "Sequence_ID": record.id,
                    "Taxon": taxon_name,
                    "Gene": gene,
                    "Domain": domain,
                    "Domain": substrate,
                    "Source": "Experiment"
                })
        except Exception as e:
            print(f"Error parsing sequences in {filepath}: {e}")

    return all_data

def main():
    parser = argparse.ArgumentParser(description="Process FASTA files into an annotated table.")
    parser.add_argument("-i", "--input", nargs='+', required=True, help="Input FASTA file(s)")
    parser.add_argument("-fi", "--fmt_input", required=True, help="Format for input filename wildcards 'exp' and 'gene'")
    parser.add_argument("-ei", "--exp_info", required=True, help="Input table with (rcsb_id,substrate,taxname)")
    parser.add_argument("-gt", "--gene_table", required=True, help="Table with (gene_name,substrate)")
    parser.add_argument("-tn", "--ttaxnames", nargs='+', required=True, help="List of taxnames to group things by.")
    parser.add_argument("-o", "--output", required=True, help="Output filename (.csv or .tsv)")
    args = parser.parse_args()

    # Initialize NCBI database and cache
    print("Initializing NCBI Database...")
    lookup_service = TaxonCache(ttaxnames=args.ttaxnames)

    # Get info 
    with open(args.exp_info, 'r') as f:
        dict_exp_info = yaml.safe_load(f)

    df = pd.read_csv(args.gene_table)
    dict_gene_substrate = dict(zip(df['gene_name'].values, df['substrate'].values))

    # Process files
    print(f"Processing {len(args.input)} files...")
    data = process_files(
        args.input, args.fmt_input, 
        lookup_service, dict_exp_info, dict_gene_substrate
    )
    if not data:
        print("No data extracted. Check your file formats and paths.")
        return

    # Aggregate with Pandas
    df = pd.DataFrame(data)

    # Determine delimiter based on file extension
    if args.output.lower().endswith('.csv'):
        df.to_csv(args.output, index=False)
    else:
        df.to_csv(args.output, sep='\t', index=False)

    print(f"Successfully saved results to {args.output}")

if __name__ == "__main__":
    main()