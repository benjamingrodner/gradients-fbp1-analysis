rule hmmsearch_iso:
    input:
        fn_hmm = lambda w: config['dict_gene_hmmprofile'][w.gene],
        fn_seqs = fmt_iso_assembly,
    output:
        fn_table_hmm = fmt_table_hmm_iso,
        fn_table_hmm_domain = fmt_table_hmm_domain_iso,
        fn_stdout_hmm = temp(fmt_stdout_hmm_iso),
        fn_hmm_hitnames = fmt_hmm_hitnames_iso,
    params:
        thresh_score = config['hmmsearch']['thresh_score']
    log:
        "logs/hmmsearch_iso/{iso}/{gene}.log"
    benchmark:
        "benchmarks/hmmsearch_iso/{iso}/{gene}.benchmark.txt"
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
