#!/usr/bin/env python3

# Adapted from group_env_hitnames.py


import argparse
import pandas as pd
import yaml
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="Genomic Metadata Integration via Ibis & DuckDB")
    parser.add_argument("--fn_best_hit", required=True, help="Path to CSV (target_name, source_file)")
    parser.add_argument("--fn_config", required=True, help="Config filename with dict_gene_hmmprofile mapping gene name to hmmer filename")
    parser.add_argument("--fmt_out", required=True, help="Format string to write to the output, must contain wildcard 'gene'")
    args = parser.parse_args()

    df = pd.read_csv(args.fn_best_hit)

    with open(args.fn_config, 'r') as f:
        config = yaml.safe_load(args.fn_config)
    
    dict_gene_hmmprofile = config['dict_gene_hmmprofile']
    dict_hmmprofile_gene = {h:g for g, h in dict_gene_hmmprofile.items()}

    dict_gene_hitnames = defaultdict(list)
    for hitname, hmmprofile in df[['target_name','source_file']].values:
        gene = dict_hmmprofile_gene[hmmprofile]
        dict_gene_hitnames[gene].append(hitname)

    # Save outputs to separate files
    for gene, hitnames in dict_gene_hitnames.items():
        filename = args.fmt_out.format(gene=gene)
        with open(filename, "w") as f:
            f.write("\n".join(hitnames))
        print(f"Grouped hitnames for {gene} saved to {filename}.")

if __name__ == "__main__":
    main()