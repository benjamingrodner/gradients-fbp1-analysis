#!/bin/bash
INDEX=(0 1 2 3 4 5 6)

GENES=(
    FBP1
    ISIP1
    ISIP2
    ISIP3
    FRE
    ALDO
    fldA
)

EXT=FBP1.fasta
ISO_DIR=/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch
DIR_OUT=/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch

for i in ${INDEX[@]}; do 
    gene=${GENES[$i]}
    ext="_${gene}.fasta"
    dir_out="${DIR_OUT}/${gene}"
    mkdir -p $dir_out
    for fasta in ${dir_out}/*${ext}; do

        clust="${fasta%${ext}}.mmseqs2"
        log="${fasta%${ext}}.mmseqs2.log"

        echo -e "gene $gene"
        echo -e "iso $iso"
        echo "fasta $fasta"
        echo "..."
        
        mmseqs easy-cluster \
            "$fasta" \
            "$clust" \
            "${dir_out}/tmp" \
            --min-seq-id "0.95" \
            -c "0.50" \
            --cov-mode 1 \
            2> "$log"
        echo "...done"

    done
done