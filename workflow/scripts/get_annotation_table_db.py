# Adapted from get_annotation_table_env.py

import argparse
import re
import os
import pandas as pd
from Bio import SeqIO
from ete4 import NCBITaxa
import gzip



class TaxonCache:
    """Handles NCBI lookups with local caching to minimize database hits."""
    def __init__(self):
        self.ncbi = NCBITaxa()
        self.cache = {}
        self.cache_tid = {}

    def get_domain(self, taxon_name):
        if taxon_name in self.cache:
            return self.cache[taxon_name]

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
            
            self.cache[taxon_name] = domain
            return domain

        except Exception as e:
            return f"Error: {str(e)}"
        
    def get_targettaxname(self, taxid, ttnames):
        """
        Retrieves the scientific name for the highest available rank 
        (phylum, superkingdom, or domain) for a given taxid.
        """
        if taxid in self.cache_tid:
            return self.cache_tid[taxid]
        try:
            lineage = self.ncbi.get_lineage(taxid)
            names = self.ncbi.get_taxid_translator(lineage)
            # Check if target tax names are in lineage
            for _, tn in names.items():
                if tn in ttnames:
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
            return tn
        except Exception:
            tn = 'No_taxon_annotation'
            self.cache_tid[taxid] = tn
            return tn


def main():
    parser = argparse.ArgumentParser(description="Process FASTA files into an annotated table.")
    parser.add_argument("-i", "--input", required=True, help="Path to CSV (target_name, source_file)")
    parser.add_argument("-u", "--uid2tax", required=True, help="Database taxon mapping")
    parser.add_argument("--ttaxnames", nargs='+', required=True, help="List of taxnames to group things by.")
    parser.add_argument("-o", "--output", required=True, help="Output filename (.csv or .tsv)")
    args = parser.parse_args()

    # Initialize NCBI database and cache
    print("Initializing NCBI Database...")
    lookup_service = TaxonCache()

    # Get mapping
    dict_marfmmdb_tax = {}
    with gzip.open(args.uid2tax, 'rt') as f:
        header = next(f)
        # for _ in range(10):
            # line = next(f)
        for line in f:
            _, header, tax, _ = line.split()
            dict_marfmmdb_tax[header] = tax

    # Process file
    df = pd.read_csv(args.input, sep='\t')
    pattern = re.compile(r"^(?P<gene>[^_]+)_hmmsearch")
    all_data = []
    try:
        for header, source in df[['target_name','source_file']]:
            tid = dict_marfmmdb_tax.get(header)
            if tid is not None:
                s = os.path.basename(source)
                gene = pattern.search(s).group('gene')
                taxon_name = lookup_service.get_targettaxname(tid, args.ttaxnames)
                domain = lookup_service.get_domain(taxon_name)
                all_data.append({
                    "Sequence_ID": header,
                    "Taxon": taxon_name,
                    "Gene": gene,
                    "Domain": domain,
                    "Source": "Database"
                })
    except Exception as e:
        print(f"Error parsing sequences in {args.input}: {e}")

    if not all_data:
        print("No data extracted. Check your file formats and paths.")
        return

    # Aggregate with Pandas
    df = pd.DataFrame(all_data)

    # Determine delimiter based on file extension
    if args.output.lower().endswith('.csv'):
        df.to_csv(args.output, index=False)
    else:
        df.to_csv(args.output, sep='\t', index=False)

    print(f"Successfully saved results to {args.output}")

if __name__ == "__main__":
    main()