
rule annotate_custom_hmms:
    input:
        fn_ann_parquet = fmt_bigtable,
        fns_tbl = [fmt_table_hmm.format(gene=g, batch='{batch}') for g in config['genes_custom']]
    output:
        fn_custom_ann = fmt_custom_ann,
    conda: 
        "../envs/python.yaml"
    threads: config['get_fasta_from_parquet']['threads']
    resources:
        mem_mb=config['get_fasta_from_parquet']['mem_mb'],
        runtime=config['get_fasta_from_parquet']['runtime'],
    params:
        script=config['dir_scripts'] + "/annotate_custom_hmms.py",
        regex = fmt_table_hmm.format(gene=r'(?P<gene>.+)', batch=r'(?P<batch>.+)'),
        join_key = config['annotate_custom_genes']['join_key'],
        col_6tr = config['annotate_custom_genes']['col_6tr']
    shell:
        """
        python {params.script} \
            --ann-parquet {input.fn_ann_parquet} \
            --output-parquet {output.fn_custom_ann} \
            --tbl-files {input.fns_tbl} \
            --regex-pattern '{params.regex}' \
            --join-key '{params.join_key}' \
            --col-6tr '{params.col_6tr}' \
            --threads {threads} \
            --mem_mb {resources.mem_mb}
        """


rule custom_ann_tables_done:
    input:
        expand(fmt_custom_ann, batch=BATCHES)
    output:
        fn_custom_ann_tables_done,
    shell:
        "echo {input} > {output}"
