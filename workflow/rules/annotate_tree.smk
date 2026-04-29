rule get_annotation_table_env:
    input:
        aggregate_env_clusters,
    output:
        fn_table_annotate_env,
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_annotation_table_env.py'
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q}
        """

rule get_annotation_table_db:
    input:
        fn_db_rep_seqs,
    output:
        fn_table_annotate_db,
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_annotation_table_db.py',
        fn_db_tax = config['fn_db_tax'],
        ttaxnames = config['ttaxnames'],
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -u {params.fn_db_tax} \
            --ttaxnames {params.ttaxnames} \
            -o {output:q}
        """
    