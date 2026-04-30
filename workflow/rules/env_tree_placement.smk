
rule mask_env_alignment_by_db_clip_log:
    input:
        aln = fn_trim_clip,
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

rule deduplicate_env_mask_alignment:
    input:
        fasta = fn_env_aligned_mask
    output:
        fasta = fn_env_aligned_mask_dedup,
        mapping = fn_env_aligned_mask_dedup_map
    log:
        err = "logs/deduplicate_env_mask_alignment.log"
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb = config['deduplicate_alignment']['mem_mb']
    params:
        script = config['dir_scripts'] + "/deduplicate_alignment.py"
    shell:
        """
        python {params.script:q} \
            --input {input.fasta:q} \
            --output {output.fasta:q} \
            --map {output.mapping:q} \
            2> {log.err:q}
        """

rule filter_env_target_genes_to_place:
    input:
        fns_target = get_target_env_rep_seqs,
        crystal = expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
        db = fn_db_rep_seqs_sub,
        manual = glob.glob(config['dir_ref_man'] + '/*'),
        aln = fn_env_aligned_mask_dedup,
    output:
        fn_env_aligned_mask_dedup_tfilt,
    log:
        "logs/filter_env_target_genes_to_place.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit seq -n -i {input.fns_target:q} {input.db:q} {input.crystal:q} {input.manual:q} \
            | seqkit grep -f - {input.aln:q} -o {output:q} \
            2> {log:q}
        """
    
rule filter_env_alignment_very_short_long:
    input:
        fn_env_aligned_mask_dedup_tfilt,
    output:
        fn_env_aligned_mask_dedup_tfilt_filt,
    log:
        "logs/filter_env_alignment_very_short_long.log",
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



rule place_env_on_tree:
    input:
        msa = fn_env_aligned_mask_dedup_tfilt_filt,
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

        echo "Done" > {output:q}
        """