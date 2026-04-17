# Generated from Gemini 3 Flash, operating in the Free tier 17 April 2026
#  with the following prompt:

# please write a python script with biopython to download rcsb fasta files based on aguments: list of ids as input to a target output dir



import os
import argparse
import requests
from Bio import SeqIO

def download_fasta(pdb_ids, output_dir):
    """
    Downloads FASTA files from RCSB PDB for a list of IDs.
    """
    # Create target directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    base_url = "https://www.rcsb.org/fasta/entry/{}/download"
    
    success_count = 0
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.strip().upper()
        if not pdb_id:
            continue
            
        print(f"Downloading {pdb_id}...")
        try:
            response = requests.get(base_url.format(pdb_id))
            
            if response.status_code == 200:
                file_path = os.path.join(output_dir, f"{pdb_id}.fasta")
                with open(file_path, "w") as f:
                    f.write(response.text)
                print(f"  [OK] Saved to {file_path}")
                success_count += 1
            else:
                print(f"  [ERROR] ID {pdb_id} not found (Status {response.status_code})")
                
        except Exception as e:
            print(f"  [FAILED] Could not download {pdb_id}: {e}")

    print(f"\nFinished. Successfully downloaded {success_count}/{len(pdb_ids)} files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download FASTA files from RCSB PDB.")
    parser.add_argument("ids", nargs="+", help="Space-separated list of PDB IDs (e.g., 4HHB 1FAT)")
    parser.add_argument("-o", "--output", default="fasta_downloads", help="Target output directory")

    args = parser.parse_args()
    
    download_fasta(args.ids, args.output)