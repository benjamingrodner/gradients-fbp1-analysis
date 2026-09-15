rule merge_db_fns_to_search:
    input:
        get_db_fns_seqs_to_search,
    output:
        temp(fn_db_seqs_to_search),
    log:
        "logs/merge_fns_to_search.log"
    benchmark:
        "benchmarks/merge_fns_to_search.benchmark.txt"
    shell:
        """
        > {output:q} 2> {log:q}
        for fn in {input:q}; do
            if [[ $fn == *.gz ]]; then
                zcat "$fn" >> {output:q} 2>> {log:q}
            else
                cat "$fn" >> {output:q} 2>> {log:q}
            fi
        done
        """ 

rule hmmsearch_db:
    input:
        fn_seqs = fn_db_seqs_to_search,
        fn_hmm = lambda w: config['dict_gene_hmmprofile'][w.gene],
    output:
        fn_table_hmm = fmt_table_hmm.format(gene='{gene}', batch='db'),
        fn_table_hmm_domain = fmt_table_hmm_domain.format(gene='{gene}', batch='db'),
        fn_stdout_hmm = temp(fmt_stdout_hmm.format(gene='{gene}', batch='db')),
        fn_hmm_hitnames = fmt_hmm_hitnames.format(gene='{gene}', batch='db'),
    params:
        thresh_score = config['hmmsearch']['thresh_score']
    log:
        "logs/hmmsearch/{gene}/db.log"
    benchmark:
        "benchmarks/hmmsearch/{gene}/db.benchmark.txt"
    threads:
        config['hmmsearch']['threads'],
    resources:
        mem_mb=config['hmmsearch']['mem_mb'],
        runtime=config['hmmsearch']['runtime'],
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


rule merge_hmmsearch_db_headers:
    input:
        [fmt_hmm_hitnames.format(gene=g, batch='db') for g in GENES_TREE],
    output:
        fn_hmm_db_hitnames_all
    shell:
        "cat {input:q} > {output:q}" 


rule get_hmms_db_best_hit:
    input:
        [fmt_table_hmm.format(gene=g, batch='db') for g in GENES_TREE],
    output:
        fn_hmms_db_best_hit,
    log:
        "logs/get_hmms_db_best_hit.log"
    params: 
        script=config['dir_scripts'] + "/get_hmms_best_hit.py",
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q} \
            2> {log:q}
        """
    
