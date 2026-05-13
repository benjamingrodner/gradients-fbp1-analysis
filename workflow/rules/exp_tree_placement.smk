rule filter_exp_target_genes_to_place:
    input:
        backbone = fn_alignment_noenv_clip_drop,
        fns_target = get_target_exp_rep_seqs,
        crystal = expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
        manual = glob.glob(config['dir_ref_man'] + '/*'),
        aln = fn_env_aligned_mask_dedup,
    output:
        fn_exp_aligned_mask_dedup_tfilt,
    log:
        "logs/filter_exp_target_genes_to_place.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit seq -n -i \
            {input.backbone:q} \
            {input.fns_target:q} \
            {input.crystal:q} \
            {input.manual:q} \
            | seqkit grep -f - {input.aln:q} -o {output:q} \
            2> {log:q}
        """
rule filter_exp_alignment_very_short_long:
    input:
        fn_exp_aligned_mask_dedup_tfilt,
    output:
        fn_exp_aligned_mask_dedup_tfilt_filt,
    log:
        "logs/filter_exp_alignment_very_short_long.log",
    params:
        frac_thresh = config['filter_alignment']['frac_thresh_short_long'],
        script = config['dir_scripts'] + "/filter_alignment.py"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q} \
            -f {params.frac_thresh:q} \
            2> {log:q}
        """



rule place_exp_on_tree:
    input:
        msa = fn_exp_aligned_mask_dedup_tfilt_filt,
        tree_done = fn_full_tree_done
    output:
        fn_place_exp_tree_done
    log:
        "logs/place_exp_on_tree.log"
    benchmark:
        "benchmarks/place_exp_on_tree.benchmark.txt"
    threads: config["build_tree"]["threads"]
    resources:
        mem_mb = config["build_tree"]["mem_mb"],
        runtime = config["build_tree"]["runtime"]
    params:
        model = config["build_tree"]['model'],
        bn_out = bn_exp_tree,
        w_out = dir_exp_tree,
        bn_tree = bn_tree,
        dir_tree = dir_tree,
    shell:
        """
        CWD=$( pwd )
        DIR_OUT="$CWD"/{params.w_out:q}
        FN_TREE={params.dir_tree:q}/RAxML_bestTree.{params.bn_tree:q}
        raxmlHPC-PTHREADS-AVX \
            -f v \
            -w "$DIR_OUT" \
            -s {input.msa} \
            -t "$FN_TREE" \
            -m {params.model} \
            -T {threads} \
            -n {params.bn_out:q} \
            2> {log}

        echo "Done" > {output:q}
        """