rule get_exp_assemblies:
    output:
        fmt_exp_assembly,
    log:
        "logs/get_exp_assemblies/{exp}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        method = lambda w: DICT_EXP[w.exp]['method_assembly'],
        fn_link = lambda w: get_exp_assm_fn_or_link(w.exp),
    shell:
        """
        TMP={output:q}.tmp
        if [[ {params.method} == 'download' ]]; then
            echo "Downloading {params.fn_link}..." 2> {log:q}
            wget -O "$TMP" {params.fn_link:q} 2>> {log:q}
        elif [[ {params.method} == 'from_author' ]]; then
            echo "Copying {params.fn_link}..." 2> {log:q}
            cp {params.fn_link} "$TMP" 2>> {log:q}
        else
            echo "Error: method {params.method} is not defined"  2> {log:q}
            exit 1
        fi
        
        if file "$TMP" | grep -q 'gzip compressed data'; then
            echo "File is gzipped. Decompressing now..." 2>> {log:q}
            gunzip -c "$TMP" > {output:q} 2>> {log:q}
            rm "$TMP" 2>> {log:q}
        else
            echo "File is already unzipped. Renaming..." 2>> {log:q}
            mv "$TMP" {output:q} 2>> {log:q}
        fi
        """

rule rename_exp_assemblies:
    input:
        fmt_exp_assembly,
    output:
        fmt_exp_assembly_rename,
    log:
        "logs/rename_exp_assemblies/{exp}.log"
    conda:
        "../envs/seqkit.yaml"
    params:
        prefix = lambda w: DICT_EXP[w.exp]['prefix'],
        regex = lambda w: get_regex_seqname(w.exp)[0],
        replace = lambda w: get_regex_seqname(w.exp)[1],
    shell:
        """
        # Custom renaming with prefix for merging experiments
        # Also regex sub if you need to match the faa names with the author counts table
        # Also sub out illegal characters for RaxML
        seqkit replace \
            -p "{params.regex}" \
            -r "{params.prefix}_{params.replace}" \
            {input:q} \
            | seqkit replace -p "[:,\)\(\[\]\']" -r "_" \
            > {output:q} \
            2> {log:q}
        """



rule translate_exp_assemblies:
    input:
        fmt_exp_assembly_rename,
    output:
        fmt_exp_assembly_6tr,
    log:
        "logs/translate_exp_assemblies/{exp}.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        SEQTYPE=$(seqkit stats {input:q} \
            | awk 'NR==2 {{print $3}}')
        if [[ "$SEQTYPE" == DNA || "$SEQTYPE" == RNA ]]; then
            transeq -auto -sformat pearson -frame 6 \
                -sequence {input:q} \
                -outseq {output:q} \
                2> {log:q}
        elif [[ "$SEQTYPE" == Protein ]]; then
            cp {input:q} {output:q} \
                2> {log:q}
        else
            echo "Error: seqtype ${{SEQTYPE}} is unknown"  \
                2> {log:q}
            exit 1
        fi
        """


rule hmmsearch_exp:
    input:
        fn_hmm = lambda w: config['dict_gene_hmmprofile'][w.gene],
        fn_seqs = fmt_exp_assembly_6tr,
    output:
        fn_table_hmm = fmt_table_hmm_exp,
        fn_table_hmm_domain = fmt_table_hmm_domain_exp,
        fn_stdout_hmm = temp(fmt_stdout_hmm_exp),
        fn_hmm_hitnames = fmt_hmm_hitnames_exp,
    params:
        thresh_score = config['hmmsearch']['thresh_score']
    log:
        "logs/hmmsearch_exp/{exp}/{gene}.log"
    benchmark:
        "benchmarks/hmmsearch_exp/{exp}/{gene}.benchmark.txt"
    conda:
        "../envs/hmmer.yaml"
    shell:
        """
        hmmsearch \
            --tblout {output.fn_table_hmm:q} \
            --domtblout {output.fn_table_hmm_domain:q} \
            -o {output.fn_stdout_hmm:q} \
            -T {params.thresh_score} \
            {input.fn_hmm:q} \
            {input.fn_seqs:q} \
            2> {log:q}

        # Get headers from hmmtable
        # ignore Grep's exit status if it didn't find any matches
        (set +o pipefail; grep -v '^#' {output.fn_table_hmm} \
            | awk '{{print $1}}' \
            > {output.fn_hmm_hitnames}) \
            2>> {log:q} 
        """

rule merge_hmmsearch_hitnames_exp:
    input:
        [fmt_hmm_hitnames_exp.format(exp='{exp}', gene=gene) for gene in GENES]
    output:
        fn_hmm_hitnames_exp_all
    shell:
        "cat {input:q} > {output:q}" 


rule get_hmms_best_hit_exp:
    input:
        [fmt_table_hmm_exp.format(exp='{exp}', gene=gene) for gene in GENES]
    output:
        fmt_hmms_best_hit_exp,
    log:
        "logs/get_hmms_best_hit_exp/{exp}.log"
    params: 
        script=config['dir_scripts'] + "/get_hmms_best_hit.py",
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q} \
            2> {log:q}
        """
    
    