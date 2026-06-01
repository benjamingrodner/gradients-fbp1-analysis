
rule get_nonenv_seqs_from_alignment:
    input: 
        aln = fn_trim_clip_dedup,
        db = fn_db_rep_seqs,
        outgroup = fn_outgroup_db_seqs_sub,
        crystal = expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
        manual = glob.glob(config['dir_ref_man'] + '/*'),
    output:
        fn_alignment_noenv,
    log:
        "logs/get_nonenv_seqs_from_alignment.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit seq -n -i {input.db:q} {input.crystal:q} {input.manual:q} {input.outgroup} \
            | seqkit grep -f - {input.aln:q} -o {output:q} \
            2> {log:q}
        """


rule remove_gappy_columns_noenv_alignment:
    input:
        fn_alignment_noenv,
    output:
        clipped = fn_alignment_noenv_clip,
        cliplog = fn_alignment_noenv_cliplog,
    log:
        "logs/remove_gappy_columns_noenv_alignment.log"
    conda:
        "../envs/clipkit.yaml"
    shell:
        """
        clipkit {input:q} \
            -m gappy \
            -g 1.0 \
            -o {output.clipped:q} \
            -l \
            2> {log:q}
        """

rule filter_alignment:
    input:
        fn_alignment_noenv_clip,
    output:
        fn_alignment_noenv_clip_filt,
    log:
        "logs/filter_alignment.log",
    params:
        frac_thresh = config['filter_alignment']['frac_thresh'],
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
   

rule get_boostrap_resamples:
    input:
        fn_alignment_noenv_clip_filt,
    output:
        expand(
            fmt_bootstrap_resample, 
            rep=range(config['fasttree_bootstraps']['n_bootstraps'])
        ),
    log:
        "logs/get_boostrap_resamples.log"
    conda:
        "../envs/goalign.yaml"
    params:
        bn = bn_bootstrap_resample
    shell:
        """
        goalign build seqboot \
            -i {input:q} \
            -o {params.bn:q} \
            -n 100 \
            -S \
            --seed 42 \
            2> {log:q}
        """

rule fasttree_bootstraps:
    input:
        fmt_bootstrap_resample,
    output:
        fmt_bootstrap_fasttree,
    log:
        "logs/fasttree_bootstraps/{rep}.log"
    benchmark:
        "benchmarks/fasttree_bootstraps/{rep}.benchmark.txt"
    resources:
        mem_mb=config['fasttree_bootstraps']['mem_mb'],
        runtime=config['fasttree_bootstraps']['runtime'],
    threads:
        config['fasttree_bootstraps']['threads'],
    shell:
        """
        export OMP_NUM_THREADS={threads}
        fasttreeMP -lg {input:q} > {output:q} \
            2> {log:q}
        """

rule merge_fasttree_bootstraps:
    input:
        expand(
            fmt_bootstrap_fasttree, 
            rep=range(config['fasttree_bootstraps']['n_bootstraps'])
        ),
    output:
        fn_bootstrap_fasttree_merged,
    shell:
        """
        cat {input:q} > {output:q}
        """

rule fasttree:
    input:
        fn_alignment_noenv_clip_filt,
    output:
        fn_fasttree,
    log:
        "logs/fasttree.log"
    benchmark:
        "benchmarks/fasttree.benchmark.txt"
    resources:
        mem_mb=config['fasttree']['mem_mb'],
        runtime=config['fasttree']['runtime'],
    threads:
        config['fasttree']['threads'],
    shell:
        """
        export OMP_NUM_THREADS={threads}
        fasttreeMP -lg -gamma {input:q} > {output:q} \
            2> {log:q}
        """
    
rule get_headers_crystal_manual:
    input:
        crystal = expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
        manual = glob.glob(config['dir_ref_man'] + '/*'),
    output:
        fn_headers_crystal_manual,
    log:
        "logs/get_headers_crystal_manual.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit seq -n -i {input.crystal:q} {input.manual:q} \
            > {output:q}
            2> {log:q}
        """

rule roguenarok:
    input:
        tree = fn_fasttree,
        boots = fn_bootstrap_fasttree_merged,
        exclude = fn_headers_crystal_manual,
    output:
        fn = fn_roguenarok,
        d = directory(dir_roguenarok),
    log:
        "logs/roguenarok.log"
    benchmark:
        "benchmarks/roguenarok.benchmark.txt"
    resources:
        mem_mb=config['roguenarok']['mem_mb'],
        runtime=config['roguenarok']['runtime'],
    conda:
        "../envs/roguenarok.yaml"
    params:
        bn = bn_tree,
    shell:
        """
        mkdir -p {output.d:q}
        RogueNaRok \
            -i {input.boots:q} \
            -t {input.tree:q} \
            -n {params.bn:q} \
            -w {output.d:q} \
            -c 50 \
            -x {input.exclude:q} \
            2> {log:q}
        """

rule pick_rogues_to_drop:
    input:
        fn_roguenarok,
    output:
        fn_rogues_to_drop,
    log:
        "logs/pick_rogues_to_drop.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + '/pick_rogues_to_drop.py',
        frac = config['roguenarok']['frac_of_improvement_to_use']
    shell:
        """
        python {params.script:q} {input:q} {output:q} \
            -f {params.frac} \
            2> {log:q}
        """

rule drop_seqs_from_alignment:
    input:
        rogues = fn_rogues_to_drop,
        aln = fn_alignment_noenv_clip_filt,
    output:
        fn_alignment_noenv_clip_filt_drop,
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -v \
            -f {input.rogues:q} \
            {input.aln:q} \
            -o {output:q}
        """


rule full_tree_analysis:
    input:
        fn_alignment_noenv_clip_filt_drop,
    output:
        f = fn_full_tree_done,
        # d = directory(dir_raxml),
    log:
        "logs/full_tree_analysis.log"
    benchmark:
        "benchmarks/full_tree_analysis.benchmark.txt"
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    threads:
        config['build_tree']['threads'],
    conda:
        "../envs/raxml_ng.yaml"
    params:
        d = dir_raxml,
        model = config['build_tree']['model'],
        bn = bn_tree,
        n_boot = lambda w: config['build_tree']['n_bootstraps'],
        seed = config['build_tree']['seed'],
    shell:
        """
        PREFIX={params.d:q}/{params.bn:q}
        raxml-ng \
            --all \
            --msa {input:q} \
            --model {params.model} \
            --threads {threads} \
            --workers auto \
            --prefix "$PREFIX" \
            --seed {params.seed} \
            --bs-trees {params.n_boot}
            2> {log:q}
        echo "Done" > {output.f:q}
        """


rule tree_extra_ml_and_bootstrapping:
    input:
        fn_alignment_noenv_clip_filt_drop,
    output:
        f = fn_extra_ml_and_bootstraps_done,
        # d = directory(dir_raxml),
    log:
        "logs/tree_extra_ml_and_bootstrapping.log"
    benchmark:
        "benchmarks/tree_extra_ml_and_bootstrapping.benchmark.txt"
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    threads:
        config['build_tree']['threads'],
    conda:
        "../envs/raxml_ng.yaml"
    params:
        d = dir_raxml,
        model = config['build_tree']['model'],
        bn = bn_extra_ml_and_bootstraps,
        n_boot = lambda w: config['build_tree']['n_bootstraps_extra'],
        trs = lambda w: config['build_tree']['n_trees_extra'],
        seed = config['build_tree']['seed_extra'],
    shell:
        """
        PREFIX={params.d:q}/{params.bn:q}
        raxml-ng \
            --all \
            --msa {input:q} \
            --model {params.model} \
            --threads {threads} \
            --workers auto \
            --prefix "$PREFIX" \
            --seed {params.seed} \
            --tree {params.trs} \
            --bs-trees {params.n_boot} \
            --bs-metric fbp,tbe
            2> {log:q}
        echo "Done" > {output.f:q}
        """
        # """
        # CWD=$( pwd )
        # DIR_TREE="$CWD"/{output.d:q}
        # raxmlHPC-PTHREADS-AVX \
        #     -s {input:q} \
        #     -w "$DIR_TREE" \
        #     -m {params.model} \
        #     -T {threads} \
        #     -n {params.bn:q} \
        #     -f a \
        #     -x 42 \
        #     -p 42 \
        #     -# {params.n_boot} \
        #     2> {log:q}
        # echo "Done" > {output.f:q}
        # """