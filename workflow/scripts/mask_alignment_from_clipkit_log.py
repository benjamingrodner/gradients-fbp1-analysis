# Generated with Gemini 3 Flash, operating in the Free tier
# using the following prompt:

# "Please generate a complete Snakemake pipeline for phylogenetic placement based on the following specifications.
# 1. Directory Structure:

# Snakefile in the root.
# Conda environments in ../envs/ (python.yaml, raxml.yaml, hmmer.yaml, seqkit.yaml).
# External Python scripts in scripts/.
# Logs in logs/ and intermediate data in work/.
# Output files written to results/
# 2. Rules to Include:

# merge_seqs_to_place: Uses a Python input function to glob a directory, filter by a target_genes list in the config, and cat them into one FASTA. Redirect stderr to log.
# build_alignment_hmm: Runs hmmbuild using hmmer.yaml. Assign threads,  RAM, and runtime with reference to the config.
# align_seqs_to_place: Runs hmmalign with --trim and --mapali (referencing the original MSA). Use hmmer.yaml, assign threads, and RAM from the config.
# trim_alignment_by_index: Uses seqkit export fa piped into seqkit subseq to hard-trim the alignment based on 1-based start/end indices from the config. Use seqkit.yaml.
# mask_alignment_by_log: Calls a standalone Python script scripts/mask_alignment.py. Use python.yaml.
# append_alignment: Simple cat of the reference alignment and the masked query alignment.
# place_seqs: Runs raxml-ng --epa. Use raxml.yaml, and assign threads,  RAM, and  runtime from the config. Move the resulting .jplace to the final results folder.
# 3. Script Requirements:

# Provide the full code for scripts/mask_alignment.py.
# The script must use argparse, contain a main() function, and use Bio.AlignIO to mask columns with '-' based on indices found in a ClipKit log file.
# Include logic to handle potential 1-based to 0-based index conversion.
# 4. Global Requirements:

# All shell commands must redirect errors using 2> {log}.
# The all rule should target results/final_placement.jplace.
# Include a sample config.yaml structure at the top of the response."

import argparse
import sys
from Bio import AlignIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np

def mask_alignment(input_file, clip_log, output_file):
    df = pd.read_csv(clip_log, sep='\s', header=None)
    df.columns = ['colnum','keep_trim','reason','parsimony']
    bool_keep = (df['keep_trim'].values == 'keep')

    alignment = AlignIO.read(input_file, "fasta")
    
    for record in alignment:
        seq_arr = np.array(record.seq)[bool_keep]
        record.seq = Seq("".join(seq_arr))
        
    AlignIO.write(alignment, output_file, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Mask alignment columns via ClipKit log.")
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-l", "--log_file", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    
    try:
        mask_alignment(args.input, args.log_file, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

# import argparse
# from Bio import AlignIO
# from Bio.Seq import MutableSeq
# from Bio.Align import MultipleSeqAlignment

# def parse_clipkit_log(log_path):
#     """
#     Parses ClipKit log to find kept/trimmed sites.
#     Assumes log format contains indices of sites.
#     """
#     keep_indices = []
#     with open(log_path, 'r') as f:
#         for line in f:
#             # Example logic: extract site index from ClipKit 'keep' log
#             if "keep" in line.lower():
#                 parts = line.split()
#                 # Adjust based on specific ClipKit log version
#                 try:
#                     idx = int(parts[1]) 
#                     keep_indices.append(idx)
#                 except (ValueError, IndexError):
#                     continue
#     return keep_indices

# def main():
#     parser = argparse.ArgumentParser(description="Mask FASTA alignment based on ClipKit log.")
#     parser.add_argument("--input", required=True, help="Input FASTA alignment")
#     parser.add_argument("--log", required=True, help="ClipKit log file")
#     parser.add_argument("--output", required=True, help="Output masked FASTA")
#     args = parser.parse_args()

#     # Load alignment
#     alignment = AlignIO.read(args.input, "fasta")
    
#     # In a real scenario, you'd identify which columns were 'trimmed' 
#     # and replace them with '-' in the query sequences.
#     # Here we assume the log provides 1-based indices of sites to MASK or KEEP.
    
#     masked_records = []
#     for record in alignment:
#         seq_list = list(record.seq)
        
#         # Example: If the log lists sites that SHOULD have been removed,
#         # we mask them with '-' to maintain alignment length for EPA.
#         # This logic assumes we are matching indices 1-to-1.
        
#         # Placeholder for masking logic:
#         # for i in range(len(seq_list)):
#         #     if i + 1 not in keep_indices:
#         #         seq_list[i] = "-"
        
#         record.seq = "".join(seq_list)
#         masked_records.append(record)

#     new_alignment = MultipleSeqAlignment(masked_records)
#     AlignIO.write(new_alignment, args.output, "fasta")

# if __name__ == "__main__":
#     main()