rule download_exp_assembly:
    output:
        directory(dir_exp_download_assembly),
    conda:
        "../envs/sra_tools.yaml"
    params:
        link = lambda w: DICT_EXP[w.exp_download_assembly]['link_assembly']
    shell:
        """
        wget -P {output:q} {params.link:q}
        """

rule merge_exp_download_assembly_dirnames:
    input:
        get_dirs_download_assembly,
    output:
        fn_download_assemblies_done,
    shell:
        """
        echo {input:q} > {output:q}
        """

rule hmmsearch_exp:
    input:
        fn_seqs = lambda w: get_fn_exp_assembly(w.exp),
        fn_hmm = lambda w: config['dict_gene_hmmprofile'][w.gene],
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
            --cpu {threads} \
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
        expand(fmt_hmm_hitnames_exp, gene=GENES)
    output:
        fn_hmm_hitnames_exp_all
    shell:
        "cat {input:q} > {output:q}" 


rule get_hmms_best_hit_exp:
    input:
        expand(fmt_table_hmm_exp, gene=GENES)
    output:
        fmt_hmms_best_hit_exp,
    log:
        "logs/get_hmms_best_hit_exp.log"
    params: 
        script=config['dir_scripts'] + "/get_hmms_best_hit.py",
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q} \
            2> {log:q}
        """
    
rule group_exp_hitnames:
    input:
        fmt_hmms_best_hit_exp,
    output:
        expand(fmt_exp_hitnames_grouped_gene),
    log:
        "logs/group_exp_hitnames.log"
    params:
        script = config['dir_scripts'] + "/group_exp_hitnames.py",
        fmt_out = lambda w: fmt_exp_hitnames_grouped_gene.format(
            exp=w.exp, gene='{gene}'
        ),
        fn_config = config['fn_config'],
    shell:
        """
        python {params.script:q} \
            {input:q} {params.fn_config:q} {params.fmt_out:q} \
            2> {log:q}
        """
    