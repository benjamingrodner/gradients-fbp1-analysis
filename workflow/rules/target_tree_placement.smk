rule mask_alignment_by_db_clip_log:
    input:
        aln = fn_trim_clip,
        log_file = fn_alignment_noenv_cliplog,
    output:
        masked = fn_env_aligned_mask
    log:
        "logs/mask_alignment_by_db_clip_log.log"
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

rule filter_target_genes_to_place:
    input:
        backbone = fn_alignment_noenv_clip_filt_drop,
        fns_env = get_target_env_rep_seqs,
        fns_exp = get_target_exp_rep_seqs,
        crystal = expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
        manual = glob.glob(config['dir_ref_man'] + '/*'),
        aln = fn_env_aligned_mask_dedup,
    output:
        fn_target_aligned_mask_dedup_tfilt,
    log:
        "logs/filter_exp_target_genes_to_place.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit seq -n -i \
            {input.backbone:q} \
            {input.fns_env:q} \
            {input.fns_exp:q} \
            {input.crystal:q} \
            {input.manual:q} \
            | seqkit grep -f - {input.aln:q} -o {output:q} \
            2> {log:q}
        """

rule filter_target_alignment_very_short:
    input:
        fn_target_aligned_mask_dedup_tfilt,
    output:
        fn_target_aligned_mask_dedup_tfilt_filt,
    log:
        "logs/filter_target_alignment_very_short.log",
    params:
        frac_thresh = config['filter_alignment']['frac_thresh_short'],
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


# TODO: Substitute in epa-ng, requires splitting new from old seqs
# TODO: add an option to add seqs to alignment with papara

rule place_target_on_tree:
    input:
        msa = fn_target_aligned_mask_dedup_tfilt_filt,
        tree_done = fn_extra_ml_and_bootstraps_done,
        # d = dir_raxml,
    output:
        d = directory(dir_target_pl_tree_raxml),
        done = fn_place_target_tree_done,
    log:
        "logs/place_target_on_tree.log"
    benchmark:
        "benchmarks/place_exp_on_tree.benchmark.txt"
    threads: config["tree_placement"]["threads"]
    resources:
        mem_mb = config["tree_placement"]["mem_mb"],
        runtime = config["tree_placement"]["runtime"]
    params:
        model = config["build_tree"]['model_old'],
        bn_out = bn_target_pl_tree,
        w_out = dir_exp_tree,
        fn_tree = fn_extra_tree_support,
    shell:
        """
        CWD=$( pwd )
        DIR_OUT="$CWD"/{output.d:q}
        mkdir -p "$DIR_OUT"
        raxmlHPC-PTHREADS-AVX \
            -f v \
            -w "$DIR_OUT" \
            -s {input.msa} \
            -t {params.fn_tree} \
            -m {params.model} \
            -T {threads} \
            -n {params.bn_out:q} \
            2> {log}

        echo "Done" > {output.done:q}
        """