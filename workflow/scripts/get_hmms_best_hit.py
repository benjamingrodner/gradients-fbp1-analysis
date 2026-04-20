#!/usr/bin/env python3

# Generated from Gemini 3 Flash, operating in the Free tier 20 April 2026

# "Please write a Python script that processes a list of HMMER tblout files to find the best hit for each target.

# Requirements:

# Input: Use argparse to accept a list of input files (via -i or --input) and an output file path (via -o or --output).

# Logic: * Iterate through each HMMER table.

# Parse the target_name (column 1) and the score (column 6).

# For each unique target_name, keep only the entry with the highest numerical score across all provided files.

# Track which specific filename provided that highest score.

# Parsing: Ensure the script handles HMMER’s metadata header lines (which start with #) and deals with variable whitespace delimiters.

# Output: Save the final result as a TSV containing the target_name, the max_score, and the source_file.

# Efficiency: Use a dictionary-based approach to track best hits to keep memory usage low, rather than loading all files into one large dataframe at once."

# Unedited...could be buggy

import argparse
import sys

def parse_hmmer_tblout(input_files, output_path):
    # Dictionary structure: {target_name: (max_score, source_file)}
    best_hits = {}

    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    # Skip metadata header lines
                    if line.startswith('#'):
                        continue
                    
                    # HMMER tblout files use variable whitespace
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    
                    target_name = parts[0]
                    try:
                        # Column 6 is the score
                        score = float(parts[5])
                    except ValueError:
                        continue

                    # Update dictionary if target is new or if current score is higher
                    if target_name not in best_hits or score > best_hits[target_name][0]:
                        best_hits[target_name] = (score, file_path)
        
        except FileNotFoundError:
            print(f"Warning: File {file_path} not found. Skipping.", file=sys.stderr)

    # Write results to TSV
    with open(output_path, 'w') as out:
        # Header
        out.write("target_name\tmax_score\tsource_file\n")
        for target, (score, source) in best_hits.items():
            out.write(f"{target}\t{score}\t{source}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Filter multiple HMMER tblout files for the best hit per target."
    )
    parser.add_argument(
        "-i", "--input", 
        nargs='+', 
        required=True, 
        help="Space-separated list of HMMER tblout files."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path for the output TSV file."
    )

    args = parser.parse_args()
    parse_hmmer_tblout(args.input, args.output)
    print(f"Done. Best hits saved to {args.output}")

if __name__ == "__main__":
    main()