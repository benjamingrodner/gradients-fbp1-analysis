# Generated with Gemini 3 Flash, operating in the Free tier 
# 29 Apr 2026
# using the following prompt:

####
# "Write a Python script that processes a list of FASTA files into a structured table (CSV or TSV). The script should include a main() function and use argparse for command-line inputs (-i for input files, -o for output).
# Requirements:

# Filename Parsing: Use re.search to parse filenames (assumed format: TaxonName_GeneName.fasta) to extract the Taxon and Gene name.
# Sequence Parsing: Use Biopython’s SeqIO to iterate through the FASTA files. For every sequence header, create a row containing the sequence ID.
# Taxonomic Annotation: Use the ETE4 NCBITaxa toolkit to look up the Domain of the extracted Taxon name (not superkingdom). Implement a cache for these lookups to ensure efficiency.
# Constant Metadata: Add a column named Source with the constant value "Environmental metatranscriptome" for all rows.
# Output: Use Pandas to aggregate the data and save it to the specified output filename. If the filename ends in .csv, save as a comma-separated file; otherwise, default to tab-separated (TSV).
# Please include error handling for missing files or taxa not found in the NCBI database."
####

import argparse
import re
import os
import pandas as pd
from Bio import SeqIO
from ete4 import NCBITaxa

class TaxonCache:
    """Handles NCBI lookups with local caching to minimize database hits."""
    def __init__(self):
        self.ncbi = NCBITaxa()
        self.cache = {}

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

def process_files(input_files, taxon_lookup):
    all_data = []
    
    # Regex pattern: TaxonName_GeneName.fasta
    pattern = re.compile(r"^taxgene_(?P<taxon>.+)_(?P<gene>[^_]+)_hmmsearch")

    for filepath in input_files:
        if not os.path.exists(filepath):
            print(f"Warning: File {filepath} not found. Skipping.")
            continue

        filename = os.path.basename(filepath)
        match = pattern.search(filename)
        
        if not match:
            print(f"Warning: Filename {filename} does not match expected format. Skipping.")
            continue

        taxon_name = match.group('taxon')
        gene_name = match.group('gene')
        domain = taxon_lookup.get_domain(taxon_name)

        # Parse FASTA sequences
        try:
            for record in SeqIO.parse(filepath, "fasta"):
                all_data.append({
                    "Sequence_ID": record.id,
                    "Taxon": taxon_name,
                    "Gene": gene_name,
                    "Domain": domain,
                    "Source": "Environmental metatranscriptome"
                })
        except Exception as e:
            print(f"Error parsing sequences in {filename}: {e}")

    return all_data

def main():
    parser = argparse.ArgumentParser(description="Process FASTA files into an annotated table.")
    parser.add_argument("-i", "--input", nargs='+', required=True, help="Input FASTA file(s)")
    parser.add_argument("-o", "--output", required=True, help="Output filename (.csv or .tsv)")
    args = parser.parse_args()

    # Initialize NCBI database and cache
    print("Initializing NCBI Database...")
    lookup_service = TaxonCache()

    # Process files
    print(f"Processing {len(args.input)} files...")
    data = process_files(args.input, lookup_service)

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