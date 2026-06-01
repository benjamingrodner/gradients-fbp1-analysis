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
        seqs_rep = fn_db_rep_seqs_sub,
        seqs_out = fn_outgroup_db_seqs_sub,
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
            -i {input.seqs_rep:q} {input.seqs_out:q} \
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
        
rule get_annotation_table_exp:
    input:
        get_target_exp_rep_seqs,
    output:
        fn_table_annotate_exp,
    log:
        "logs/get_annotation_table/exp.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_annotation_table_exp.py',
        ei = config['fn_experiment_info'],
        gt = config['fn_gene_info'],
        tn = config['ttaxnames'],
        fi = lambda w: fmt_exp_rep_seqs,
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -fi {params.fi} \
            -ei {params.ei:q} \
            -gt {params.gt:q} \
            -tn {params.tn:q} \
            -o {output:q}
        """

rule merge_annotation_tables:
    input:
        fn_table_annotate_env,
        fn_table_annotate_exp,
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


# rule get_itol_annotations:
#     input:
#         annot = fn_table_annotations,
#         tree_done = fn_place_env_tree_done,
#     output:
#         fn_taxon_colorstrip,
#         fn_domain_colorstrip,
#         fn_substrate_treecolors,
#         fn_source_treecolors,
#         fn_fbp1_colorstrip,
#         # fn_crystal_symbol,
#         # fn_bootstrap_symbol,
#     conda:
#         "../envs/python.yaml"
#     params:
#         script = config['dir_scripts'] + '/get_itol_annotations.py',
#         tree = fn_place_env_tree,
#         d = dir_annot,
#         y = config['fn_config']
#     shell:
#         """
#         python {params.script:q} \
#             -t {params.tree:q} \
#             -a {input.annot:q} \
#             -y {params.y:q} \
#             -d {params.d:q}
#         """

# rule get_itol_annotations_exp:
#     input:
#         annot = fn_table_annotations,
#         tree_done = fn_place_exp_tree_done,
#     output:
#         fn_taxon_colorstrip_exp,
#         fn_domain_colorstrip_exp,
#         fn_substrate_treecolors_exp,
#         fn_source_treecolors_exp,
#         fn_fbp1_colorstrip_exp,
#         # fn_crystal_symbol,
#         # fn_bootstrap_symbol,
#     conda:
#         "../envs/python.yaml"
#     params:
#         script = config['dir_scripts'] + '/get_itol_annotations.py',
#         tree = fn_place_exp_tree,
#         d = dir_annot_exp,
#         y = config['fn_config']
#     shell:
#         """
#         python {params.script:q} \
#             -t {params.tree:q} \
#             -a {input.annot:q} \
#             -y {params.y:q} \
#             -d {params.d:q}
#         """

rule get_itol_annotations_target_pl:
    input:
        annot = fn_table_annotations,
        tree_done = fn_extra_ml_and_bootstraps_done,
    output:
        fn_target_pl_annot_done,
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/get_itol_annotations.py',
        tree = fn_place_target_tree,
        fn_tree_jplace = fn_place_target_jplace,
        support = fn_extra_tree_support,
        d = dir_annot_itol,
        y = config['fn_config']
    shell:
        """
        python {params.script:q} \
            -t {params.tree:q} \
            -tj {params.fn_tree_jplace} \
            -a {input.annot:q} \
            -y {params.y:q} \
            -d {params.d:q} \
            -s {params.support:q}
        echo "done" > {output:q}
        """


# rule add_bootstraps_to_jplace:
#     input:
#         tree_done = fn_full_tree_done,
#         place_done = fn_place_target_tree_done,
#     output:
#         jplace_out = fn_place_target_jplace_support,
#     conda:
#         "../envs/gappa.yaml"
#     params:
#         jplace = fn_place_target_jplace,
#         support = fn_tree_support,
#         bn_out = bn_target_pl_tree,
#     run:
#         import json
#         import re

#         # 1. Read the reference tree that contains the bootstrap values
#         with open(params.support, "r") as f:
#             # Strip whitespace/newlines common in Newick files
#             bootstrap_tree_string = f.read().strip()
        
#         # 2. Load the original .jplace file
#         with open(params.jplace, "r") as f:
#             jplace_data = json.load(f)
            
#         # Crucial Step: Verify or adapt edge labels if necessary.
#         # If your bootstrap tree completely lacks the edge labels (e.g., '[#1]'),
#         # you must make sure the topologies match exactly.
        
#         # 3. Swap out the backbone tree string
#         jplace_data["tree"] = bootstrap_tree_string
        
#         # 4. Write back the updated JSON structure
#         with open(output.jplace_out, "w") as f:
#             json.dump(jplace_data, f, indent=2)
#         print(f"Successfully injected bootstrap tree into {output_jplace}")