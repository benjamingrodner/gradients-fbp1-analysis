#!/usr/bin/env python3

# Generated from Gemini 3 Flash, operating in the Free tier 20 April 2026


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

import argparse
import pandas as pd
import ibis
from ete3 import NCBITaxa

def get_targettaxname(taxid, ncbi):
    """
    Traverses the NCBI lineage to find the highest available rank 
    among phylum, superkingdom, or domain.
    """
    if pd.isna(taxid) or taxid == 0:
        return 'No_taxon_annotation'
    
    try:
        lineage = ncbi.get_lineage(taxid)
        names = ncbi.get_taxid_translator(lineage)
        ranks = ncbi.get_rank(lineage)
        
        # Rank priority mapping
        rank_priority = ['phylum', 'superkingdom', 'kingdom'] 
        
        # Inverted lineage to check from lower to higher or specific ranks
        # Create a map of rank_name -> scientific_name
        rank_map = {ranks[tid]: names[tid] for tid in lineage}
        
        for rank in rank_priority:
            if rank in rank_map:
                return rank_map[rank]
                
        return 'No_taxon_annotation'
    except (ValueError, KeyError):
        return 'No_taxon_annotation'

def main():
    parser = argparse.ArgumentParser(description="Genomic metadata integration script.")
    parser.add_argument("--table1", required=True, help="Path to CSV (target_name, source_file)")
    parser.add_argument("--parquet", required=True, help="Path to Parquet (target_id, contig_name_6tr, taxid)")
    parser.add_argument("--threads", type=int, default=4, help="CPU core allocation")
    parser.add_argument("--memory", default="16GB", help="RAM allocation (e.g. '16GB')")
    parser.add_argument("--b_string", default="", help="Conditional naming logic string")
    args = parser.parse_args()

    # 1. Resource Management & Ibis Configuration
    con = ibis.duckdb.connect()
    con.execute(f"SET threads TO {args.threads};")
    con.execute(f"SET memory_limit TO '{args.memory}';")

    # 2. Data Integration
    t1 = con.read_csv(args.table1)
    pq = con.read_parquet(args.parquet)

    # Inner join on contig_name_6tr == target_name
    joined = pq.join(t1, pq.contig_name_6tr == t1.target_name)
    
    # Materialize to Pandas for row-wise logic and ETE3 integration
    df = joined.execute()

    # 3. Taxonomy Memoization (Two-Pass Logic)
    ncbi = NCBITaxa()
    unique_taxids = df['taxid'].unique()
    
    # Pass 1: Build the lookup dictionary
    dict_tt_ttaxname = {
        tid: get_targettaxname(tid, ncbi) for tid in unique_taxids
    }

    # 4. Mapping & Multi-hit Resolution
    results = {} # Nested Dictionary
    
    # Grouping by target_id to resolve multi-mappings
    for tid, group in df.groupby('target_id'):
        
        # Resolve multi-hit: select row with longest lineage
        if len(group) > 1:
            # We calculate lineage length on the fly for resolution
            group = group.copy()
            group['lineage_len'] = group['taxid'].apply(lambda x: len(ncbi.get_lineage(x)) if x != 0 else 0)
            target_row = group.loc[group['lineage_len'].idxmax()]
        else:
            target_row = group.iloc[0]

        # Taxonomic name from memoized dict
        ttax = dict_tt_ttaxname.get(target_row['taxid'], 'No_taxon_annotation')
        source_file = target_row['source_file']
        
        # Naming logic
        out_name = target_row['contig_name_6tr'] if 'PA' in args.b_string else target_row['target_id']

        # 5. Populate Multilevel Nested Dictionary
        if ttax not in results:
            results[ttax] = {}
        if source_file not in results[ttax]:
            results[ttax][source_file] = []
            
        results[ttax][source_file].append(out_name)

    # Summary Output
    print("\n--- Summary of Dictionary Contents ---")
    for taxon, sources in results.items():
        total_contigs = sum(len(c_list) for c_list in sources.values())
        print(f"Taxon: {taxon} | Sources: {len(sources)} | Total Contigs: {total_contigs}")

if __name__ == "__main__":
    main()