#!/bin/bash
# INDEX=(0 1 2 3 4 5 6)

# GENES=(
#     FBP1
#     ISIP1
#     ISIP2
#     ISIP3
#     FRE
#     ALDO
#     fldA
# )

# EXT=.Trinity.6tr.bf100.fasta.gz
ISO_DIR=/Users/benjamingrodner/work/armbrust/data/prevalence/clustering
DIR_OUT=/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2
mkdir -p $DIR_OUT
for fasta in ${ISO_DIR}/*.fasta; do 
    bn=$(basename $fasta)
    clust="${DIR_OUT}/${bn%.fasta}.mmseqs2"
    log="${DIR_OUT}/${bn%.fasta}.log"

    echo "$bn"
    echo "..."

    mmseqs easy-cluster \
        "$fasta" \
        "$clust" \
        "${DIR_OUT}/tmp" \
        --min-seq-id "0.95" \
        -c "0.50" \
        --cov-mode 1 \
        2> "$log"

    echo "...done"

done