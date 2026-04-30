rule get_annotation_table_env:
    input:
        aggregate_env_clusters,
    output:
        fn_table_annotate_env,
    log:
        "logs/get_annotation_table/env.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_annotation_table_env.py',
        t = config['fn_gene_info'],
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -t {params.t:q} \
            -o {output:q}
        """

rule get_annotation_table_db:
    input:
        hmm_table = fn_hmms_best_hit,
        seqs = fn_db_rep_seqs,
    output:
        fn_table_annotate_db,
    log:
        "logs/get_annotation_table/db.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_annotation_table_db.py',
        fn_db_tax = config['fn_db_tax'],
        tn = config['ttaxnames'],
        gt = config['fn_gene_info'],
    shell:
        """
        python {params.script:q} \
            -i {input.seqs:q} \
            -ht {input.hmm_table:q} \
            -u {params.fn_db_tax:q} \
            -gt {params.gt:q} \
            -tn {params.tn:q} \
            -o {output:q}
        """
    
rule get_annotation_table_crystal:
    input:
        expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
    output:
        fn_table_annotate_crystal,
    log:
        "logs/get_annotation_table/crystal.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_annotation_table_crystal.py',
        rt = config['fn_rcsb_table'],
        tn = config['ttaxnames']
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -rt {params.rt:q} \
            -tn {params.tn:q} \
            -o {output:q}
        """

rule merge_annotation_tables:
    input:
        fn_table_annotate_env,
        fn_table_annotate_db,
        fn_table_annotate_crystal,
        config['fn_ref_man_info'],
    output:
        fn_table_annotations,
    shell:
        """
        head -n 1 {input[0]:q} \
            > {output:q}
        tail -n +2 -q {input:q} \
            >> {output:q}
        """


rule get_itol_annotations:
    input:
        annot = fn_table_annotations,
        tree_done = fn_place_env_tree_done,
    output:
        fn_taxon_colorstrip,
        fn_domain_colorstrip,
        fn_substrate_treecolors,
        fn_source_treecolors,
        fn_fbp1_colorstrip,
        # fn_crystal_symbol,
        # fn_bootstrap_symbol,
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_itol_annotations.py',
        tree = fn_place_env_tree,
        d = dir_annot,
        y = config['fn_config']
    shell:
        """
        python {params.script:q} \
            -t {params.tree:q} \
            -a {input.annot:q} \
            -y {params.y:q} \
            -d {params.d:q}
        """
