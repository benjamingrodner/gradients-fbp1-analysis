# Generated from Gemini 3 Flash, operating in the Free tier 16 April 2026
#  with the following prompt:
# Please make a plan to convert the following lines of python to a script with input -hc for header_column, -h for a file with a list of headers, -i for input parquet file, -o for output fasta file:
#             t_seqs = ibis.read_parquet(fn_seqs)
#             fn_out = f'{dir_groups}/{gene}/Environmental_metatranscriptome/{ttax}_from_Environmental_metatranscriptome_{b}.faa'
#             os.makedirs(os.path.split(fn_out)[0], exist_ok=True)
#             t_seqs.filter(
#                 t_seqs[header_column].isin(headers)
#             ).select(
#                 fasta_format=(">" + t_seqs['contig_name'] + "\n" + t_seqs['aa_seq'])
#             ).to_csv(
#                 fn_out, delimiter="", quote="", header=False
#             )

import argparse
import os
import ibis

def main():
    parser = argparse.ArgumentParser(description="Filter a Parquet file and export to FASTA.")
    
    # Define CLI arguments
    parser.add_argument("-hc", "--header_column", required=True, help="Column name to filter headers by")
    parser.add_argument("-h", "--headers_file", required=True, help="Path to text file containing list of headers")
    parser.add_argument("-i", "--input", required=True, help="Input Parquet file path")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file path")
    
    args = parser.parse_args()

    # 1. Load the headers from the provided file
    if not os.path.exists(args.headers_file):
        print(f"Error: Headers file '{args.headers_file}' not found.")
        return

    with open(args.headers_file, 'r') as f:
        # Strips whitespace/newlines and ignores empty lines
        header_list = [line.strip() for line in f if line.strip()]

    # 2. Initialize ibis and read the parquet
    # Note: Ibis will use the default backend available (e.g., DuckDB or Polars)
    t_seqs = ibis.read_parquet(args.input)

    # 3. Create the output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # 4. Perform the transformation
    # We create the fasta_format column and export it
    print(f"Filtering {args.input}...")
    
    query = (
        t_seqs.filter(t_seqs[args.header_column].isin(header_list))
        .select(
            fasta_format=(">" + t_seqs['contig_name'] + "\n" + t_seqs['aa_seq'])
        )
    )

    # 5. Write to file
    # quote="" and delimiter="" ensures we get raw text output
    query.to_csv(args.output, delimiter="", quote="", header=False)
    
    print(f"Success! FASTA written to: {args.output}")

if __name__ == "__main__":
    main()