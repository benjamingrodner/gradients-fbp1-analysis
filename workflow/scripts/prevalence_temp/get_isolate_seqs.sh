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

EXT=.Trinity.6tr.bf100.fasta.gz
ISO_DIR=/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/fastas
DIR_OUT=/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch

for fasta in ${ISO_DIR}/*${EXT}; do 
    iso=$(basename ${fasta%$EXT})
    for i in ${INDEX[@]}; do 
        gene=${GENES[$i]}
        hmm=${HMMS[$i]}

        dir_out="${DIR_OUT}/${gene}"
        mkdir -p $dir_out

        names="${dir_out}/${iso}_${gene}.names"
        seqs="${dir_out}/${iso}_${gene}.fasta"
        log="${dir_out}/${iso}_${gene}.log_seqkit"

        echo -e "gene $gene"
        echo -e "iso $iso"
        echo "fasta $fasta"
        
        seqkit grep \
            -f "$names" \
            "$fasta" \
            > "$seqs" \
            2> "$log"
    done
done