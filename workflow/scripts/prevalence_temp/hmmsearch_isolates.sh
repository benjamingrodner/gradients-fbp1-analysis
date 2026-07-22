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

HMMS=(
    /Users/benjamingrodner/work/armbrust/code/gradients-fbp1-analysis/resources/fbp1.hmm
    /Users/benjamingrodner/work/armbrust/data/prevalence/custom_hmms/ISIP1.trimal025.aln.selex.hmm
    /Users/benjamingrodner/work/armbrust/data/prevalence/custom_hmms/ISIP2.trimal025.aln.selex.hmm
    /Users/benjamingrodner/work/armbrust/data/prevalence/custom_hmms/ISIP3.trimal025.aln.selex.hmm
    /Users/benjamingrodner/work/armbrust/data/prevalence/KO_profiles/K00521.hmm
    /Users/benjamingrodner/work/armbrust/data/prevalence/KO_profiles/K01623.hmm
    /Users/benjamingrodner/work/armbrust/data/prevalence/KO_profiles/K03839.hmm
)
THRESH=30
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

        tblout="${dir_out}/${iso}_${gene}.tblout"
        out="${dir_out}/${iso}_${gene}.out"
        names="${dir_out}/${iso}_${gene}.names"
        log="${dir_out}/${iso}_${gene}.log"

        echo "i $i"
        echo -e "gene $gene"
        echo -e "iso $iso"
        echo "fasta $fasta"
        echo -e "hmm $hmm"
        echo -e "tblout $tblout"
        echo -e "..."

        hmmsearch \
            --tblout "$tblout" \
            -o "$out" \
            -T $THRESH \
            "$hmm" \
            "$fasta" \
            2> "$log"

        # Get headers from hmmtable
        # ignore Grep's exit status if it didn't find any matches
        (set +o pipefail; grep -v '^#' "$tblout" \
            | awk '{print $1}' \
            > "$names") \
            2>> "$log" 

        echo -e "\done"
        echo -e "\n"
    done

done


