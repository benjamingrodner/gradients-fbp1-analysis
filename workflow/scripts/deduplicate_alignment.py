# Generated with Gemini 3 Flash, operating in the Free tier 27 Apr 2026
# using the following prompt:

# "Please generate a Snakemake rule and a companion Python script to deduplicate an MSA (Multiple Sequence Alignment) while preserving gaps and mapping duplicates.
# 1. Snakemake Rule Requirements:

# Rule Name: deduplicate_alignment
# Input: A masked FASTA alignment (work/queries_masked.fasta).
# Output: A deduplicated FASTA alignment (work/queries_dedup.fasta) and a JSON mapping file (work/dedup_map.json).
# Environment: Use a Conda environment file at ../envs/python.yaml.
# Resources: Allocate 2GB of RAM.
# Command: Execute the script using a shell command, passing the input, output, and map file as arguments. Redirect stderr to logs/deduplicate.log.
# 2. Python Script Requirements (scripts/deduplicate_msa.py):

# Framework: Use argparse for command-line inputs and wrap the logic in a main() function.
# Library: Use Bio.AlignIO to read and write the sequences.
# Logic: * The script must treat the entire sequence string (including gaps -) as the unique identifier to ensure identical alignment patterns are grouped.
# It should keep the first occurrence of a sequence as the 'representative.'
# It must generate a dictionary where the keys are the IDs of the representative sequences and the values are lists of IDs of all sequences (including the representative) that share that exact alignment string.
# Output: Write the representative sequences to a FASTA file and the dictionary to a JSON file.
# 3. Format:

# Provide the Snakemake rule code and the full Python script separately.
# Ensure the code is robust enough to handle large alignments by including basic error handling in the try/except block."

import argparse
import json
import sys
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def deduplicate_msa(input_file, output_fasta, output_json):
    """
    Deduplicates an MSA based on identical sequence strings (including gaps).
    """
    # map: { sequence_string: [record_id_1, record_id_2, ...] }
    seq_to_ids = defaultdict(list)
    # map: { representative_id: record_object }
    representative_records = {}
    
    try:
        alignment = AlignIO.read(input_file, "fasta")
        
        for record in alignment:
            seq_str = str(record.seq)
            seq_to_ids[seq_str].append(record.id)
            
            # Store the first occurrence as the representative
            if seq_str not in representative_records:
                representative_records[seq_str] = record

        # Build the final mapping dictionary: { RepresentativeID: [AllIDs] }
        final_mapping = {}
        deduplicated_records = []

        for seq_str, ids in seq_to_ids.items():
            rep_record = representative_records[seq_str]
            rep_id = rep_record.id
            final_mapping[rep_id] = ids
            deduplicated_records.append(rep_record)

        # Write deduplicated FASTA
        with open(output_fasta, "w") as f_out:
            AlignIO.write(AlignIO.MultipleSeqAlignment(deduplicated_records), f_out, "fasta")

        # Write JSON mapping
        with open(output_json, "w") as j_out:
            json.dump(final_mapping, j_out, indent=4)

    except Exception as e:
        sys.stderr.write(f"Error during MSA deduplication: {str(e)}\n")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Deduplicate MSA and map IDs.")
    parser.add_argument("-i", "--input", required=True, help="Input masked FASTA alignment")
    parser.add_argument("-o", "--output", required=True, help="Output deduplicated FASTA")
    parser.add_argument("-m", "--map", required=True, help="Output JSON mapping file")
    
    args = parser.parse_args()
    
    deduplicate_msa(args.input, args.output, args.map)

if __name__ == "__main__":
    main()