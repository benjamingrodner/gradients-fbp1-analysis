#!/bin/bash

# --- SLURM Settings ---
#SBATCH --job-name=build_table    # Job name
#SBATCH --partition=main          # Partition/Queue name
#SBATCH --nodes=1                     # Run on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per file
#SBATCH --mem=8G                      # Memory limit
#SBATCH --time=00:10:00               # Time limit (hrs:min:sec)
#SBATCH --output=slurmlog/build_table_%j.log           # Standard output and error log (%j = JobID)

# Given marferret mmdb hmmsearch hits for fbp1, what is their source entry id?
FNS_NAMES=(
    /Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/fre_hmmsearch_T30.names
)
    # /scratch/bgrodner/fbp1_paper/small_tree/results/hmmsearch/fbp1_hmmsearch_T30.names
    # /scratch/bgrodner/metaT_data/isip_hmmsearch/results/hmmsearch/isip1_hmmsearch_T30.names
    # /scratch/bgrodner/metaT_data/isip_hmmsearch/results/hmmsearch/isip2_hmmsearch_T30.names
    # /scratch/bgrodner/metaT_data/isip_hmmsearch/results/hmmsearch/isip3_hmmsearch_T30.names

FN_MERGE_NAMES='merge_contig_names-fre.txt'
FN_MERGE_NAMES_MFT='merge_contig_names-fre-mft.txt'
> $FN_MERGE_NAMES
for f in ${FNS_NAMES[@]}; do
    cat $f >> $FN_MERGE_NAMES
done

# cat /Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/* > $FN_MERGE_NAMES

grep -e '^mft' $FN_MERGE_NAMES > $FN_MERGE_NAMES_MFT

FN_AAINFO=MarFERReT.v1.1.1.proteins_info.tab.gz
FN_FBP1INFO=marfmmdb_contig_source-fre.grep.tsv
zgrep \
    -f "$FN_MERGE_NAMES_MFT" \
    "$FN_AAINFO" \
    > "$FN_FBP1INFO"

# # Get fbp1 cluster ids
# DIR_CLUSTERS=/scratch/bgrodner/iron_ko_contigs/sidero_receptors/fbp1/alignment/clustering/fbp1-group_source_tax/mmseqsi80c80covmode1
# grep \
#     -f "$FN_FBP1NAMES" \
#     "${DIR_CLUSTERS}"/source_Database*.tsv
#     > 


# cat $FN_FBP1INFO \
#     | while IFS=$'\t' read -r aa_id entry_id source_defline; do
#         echo 
#     done


