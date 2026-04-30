
rule mask_env_alignment_by_db_clip_log:
    input:
        aln = fn_trim_clip_filt_dedup,
        log_file = fn_alignment_noenv_cliplog,
    output:
        masked = fn_env_aligned_mask
    log:
        "logs/mask_env_alignment_by_db_clip_log.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + "/mask_alignment_from_clipkit_log.py"
    shell:
        "python {params.script} "
        "--input {input.aln:q} "
        "--log {input.log_file:q} "
        "--output {output.masked} 2> {log}"


rule place_env_on_tree:
    input:
        msa = fn_env_aligned_mask,
        tree_done = fn_full_tree_done
    output:
        fn_place_env_tree_done
    log:
        "logs/place_env_on_tree.log"
    benchmark:
        "benchmarks/place_env_on_tree.benchmark.txt"
    threads: config["build_tree"]["threads"]
    resources:
        mem_mb = config["build_tree"]["mem_mb"],
        runtime = config["build_tree"]["runtime"]
    params:
        model = config["build_tree"]['model'],
        bn_out = bn_env_tree,
        w_out = dir_env_tree,
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

        echo "Done" > {output.tree_done:q}
        """