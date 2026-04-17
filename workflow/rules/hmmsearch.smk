
rule get_fasta_from_parquet:
    input:
        fmt_seqs_parquet,
    output:
        temp(fmt_seqs_fasta),
    params: 
        script="scripts/parquet2fasta.py",
    log:
        "logs/get_fasta_from_parquet/{batch}.log"
    benchmark:
        "benchmarks/get_fasta_from_parquet/{batch}.benchmark.txt"
    conda: 
        "envs/python.yaml"
    threads:
        config['get_fasta_from_parquet']['threads'],
    resources:
        mem_mb=config['get_fasta_from_parquet']['mem_mb'],
        runtime=config['get_fasta_from_parquet']['runtime'],
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q} \
            -hc contig_name_6tr \
            -hf all \
            -m {resources.mem_mb} \
            -t (threads) \
            2> {log:q}
        """


rule merge_fns_to_search:
    input:
        get_fns_seqs_to_search,
    output:
        temp(fn_seqs_to_search),
    shell:
        "cat {input:q} > {output:q}" 


rule hmmsearch:
    input:
        fn_seqs = fn_seqs_to_search,
        fn_hmm = lambda w: config['dict_gene_hmmprofile'][w.gene],
    output:
        fn_table_hmm = fmt_table_hmm,
        fn_table_hmm_domain = fmt_table_hmm_domain,
        fn_stdout_hmm = temp(fmt_stdout_hmm),
        fn_hmm_hitnames = fmt_hmm_hitnames,
    params:
        thresh_score = config['hmmsearch']['thresh_score']
    log:
        "logs/hmmsearch.log"
    benchmark:
        "benchmarks/hmmsearch.benchmark.txt"
    threads:
        config['hmmsearch']['threads'],
    resources:
        mem_mb=config['hmmsearch']['mem_mb'],
        runtime=config['hmmsearch']['runtime'],
    conda:
        "envs/hmmer.yaml"
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
        (set +o pipefail; grep -v '^#' {input.fn_table_hmm} \
            | awk '{{print $1}}' \
            > {output.fn_hmm_hitnames}) \
            2>> {log:q} 
        """

rule merge_hmmsearch_headers:
    input:
        expand(fmt_hmm_hitnames, GENES)
    output:
        fn_hmm_hitnames_all
    shell:
        "cat {input:q} > {output:q}" 
