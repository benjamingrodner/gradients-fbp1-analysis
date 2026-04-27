#!/usr/bin/env python3

# Generated from Gemini 3 Flash, operating in the Free tier 20 April 2026

##############################
# "Write a Python script that integrates metadata from a CSV table with genomic mappings in a Parquet file using Ibis with a DuckDB backend.
# 1. Resource Management & CLI:

# Use argparse to accept:
# --table1: Path to a CSV containing target_name and source_file.
# --parquet: Path to a Parquet file containing target_id, contig_name_6tr, and taxid.
# --threads: Integer for CPU core allocation.
# --memory: String for RAM allocation (e.g., '16GB').
# --b_string: A string for conditional naming logic.
# Configure the Ibis DuckDB connection using the provided threads and memory arguments.
# 2. Data Integration:

# Perform an inner join using Ibis between the Parquet file and the CSV table, joining on contig_name_6tr == target_name.
# Execute the join and bring only the relevant subset into a Pandas DataFrame for complex row-wise processing.
# 3. Taxonomy Memoization (Two-Pass Logic):

# Use ete3.NCBITaxa for taxonomy lookups.
# First Pass: Create a lookup dictionary (dict_tt_ttaxname) by calling a custom function get_targettaxname for every unique taxid found in the joined data.
# The get_targettaxname function should:
# Retrieve the NCBI lineage for a taxid.
# Traverse the lineage to find the scientific name for the highest available rank among 'phylum', 'superkingdom', or 'domain'.
# Return 'No_taxon_annotation' as a fallback.
# 4. Mapping & Multi-hit Resolution:

# Iterate through each unique target_id in the joined data.
# Multi-mapping resolution: If a target_id has multiple associated rows, select the row where the taxid has the longest NCBI lineage.
# Naming Logic: If 'PA' is in the b_string, use the contig_name_6tr for the output; otherwise, use target_id.
# 5. Output Structure:

# Populate a multilevel nested dictionary where:
# Level 1 Key: The resolved taxon name (ttax).
# Level 2 Key: The source_file mapped from Table 1.
# Value: A list of the processed contig names.
# Save each list of contigs to a separate file defined by the ttax and source_file keys from the dict"
##############################

# Edited during testing 

import argparse
import pandas as pd
import ibis
from ete4 import NCBITaxa
import os

def get_targettaxname(taxid, ttnames, ncbi):
    """
    Retrieves the scientific name for the highest available rank 
    (phylum, superkingdom, or domain) for a given taxid.
    """
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        # Check if target tax names are in lineage
        for _, tn in names.items():
            if tn in ttnames:
                return tn

        # Invert ranks to search by rank name
        ranks = ncbi.get_rank(lineage)
        rank_to_id = {rank: tid for tid, rank in ranks.items()}
        # Priority order for highest available rank
        for rank_name in ['phylum', 'superkingdom', 'domain']:
            if rank_name in rank_to_id:
                tid = rank_to_id[rank_name]
                return ncbi.get_taxid_translator([tid])[int(tid)]
        
        return 'No_taxon_annotation'
    except Exception:
        return 'No_taxon_annotation'

def main():
    parser = argparse.ArgumentParser(description="Genomic Metadata Integration via Ibis & DuckDB")
    parser.add_argument("--table1", required=True, help="Path to CSV (target_name, source_file)")
    parser.add_argument("--parquet", required=True, help="Path to Parquet (target_id, contig_name_6tr, taxid)")
    parser.add_argument("--ttaxnames", nargs='+', required=True, help="List of taxnames to group things by.")
    parser.add_argument("--threads", type=int, default=4, help="CPU cores")
    parser.add_argument("--memory", default="16GB", help="RAM allocation (e.g. 16GB)")
    parser.add_argument("--b_string", required=True, help="String for conditional naming logic")
    parser.add_argument("--dir_out", required=True, help="Output directory")
    parser.add_argument("--wc_prefix", required=True, help="Wildcard prefix for parsing by future rules")
    args = parser.parse_args()

    # 1. Resource Management & Ibis/DuckDB Setup
    con = ibis.duckdb.connect(threads=args.threads, memory_limit=args.memory)
    
    table_csv = ibis.read_csv(args.table1)
    table_parquet = ibis.read_parquet(args.parquet)

    # 2. Data Integration
    # Inner join on contig_name_6tr == target_name
    jkey = 'contig_name_6tr' if 'PA' in args.b_string else 'target_id'  # Since PA seqs are 6tr and NS best frame have already been picked
    joined = table_parquet.join(
        table_csv, 
        table_parquet[jkey] == table_csv.target_name
    )
    
    # Materialize to Pandas for row-wise processing and taxonomy lookups
    df = joined.execute()

    # 3. Taxonomy Memoization
    ncbi = NCBITaxa()
    unique_taxids = df['taxid'].unique()
    dict_tt_ttaxname = {tid: get_targettaxname(tid, args.ttaxnames, ncbi) for tid in unique_taxids}
    # 4. Mapping & Multi-hit Resolution
    results_tree = {} # Level 1: ttax, Level 2: source_file, Value: list of contigs

    for target_id, group in df.groupby('target_id'):
        # Multi-hit resolution: Select taxid with the longest lineage
        if len(group) > 1:
            lineages = {}
            for tid in group['taxid'].unique():
                l = 0
                if int(tid) > 0:
                    l = len(ncbi.get_lineage(tid)) 
                lineages[tid] = l
            best_taxid = max(lineages, key=lineages.get)
            selected_row = group[group['taxid'] == best_taxid].iloc[0]
        else:
            selected_row = group.iloc[0]

        # Naming Logic
        final_contig_name = (
            selected_row['contig_name_6tr'] if 'PA' in args.b_string else target_id
        )

        # Retrieve resolved taxon name from memoized dict
        ttax = dict_tt_ttaxname.get(selected_row['taxid'], 'No_taxon_annotation')
        source_file = selected_row['source_file']

        # 5. Populate Nested Structure
        if ttax not in results_tree:
            results_tree[ttax] = {}
        if source_file not in results_tree[ttax]:
            results_tree[ttax][source_file] = []
        
        results_tree[ttax][source_file].append(final_contig_name)
    # Save outputs to separate files
    os.makedirs(args.dir_out, exist_ok=True)
    for ttax, sources in results_tree.items():
        for source_file, contigs in sources.items():
            # Sanitize filename
            sf = os.path.splitext(os.path.split(source_file)[1])[0]
            filename = f"{args.dir_out}/{args.wc_prefix}{ttax}_{sf}.txt"
            # if len(contigs) > 0:
            with open(filename, "w") as f:
                f.write("\n".join(contigs))

    print(f"Processing complete. Files saved in {args.dir_out}.")

if __name__ == "__main__":
    main()