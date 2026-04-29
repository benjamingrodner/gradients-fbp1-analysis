
rule get_nonenv_seqs_from_alignment:
    input: 
        aln = fn_trim_clip_filt_dedup,
        db = fn_db_rep_seqs,
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
        seqkit seq -n -i {input.db:q} {input.crystal:q} {input.manual:q} \
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

        

rule get_boostrap_resamples:
    input:
        fn_alignment_noenv_clip,
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
        fn_alignment_noenv_clip,
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
        fn_roguenarok,
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
        w = dir_roguenarok,
    shell:
        """
        mkdir -p {params.w:q}
        RogueNaRok \
            -i {input.boots:q} \
            -t {input.tree:q} \
            -n {params.bn:q} \
            -w {params.w:q} \
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
        aln = fn_alignment_noenv_clip
    output:
        fn_alignment_noenv_clip_drop,
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
        fn_alignment_noenv_clip_drop,
    output:
        fn_full_tree_done,
    log:
        "logs/full_tree_analysis.log"
    benchmark:
        "benchmarks/full_tree_analysis.benchmark.txt"
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    threads:
        config['build_tree']['threads'],
    params:
        model = config['build_tree']['model'],
        w = dir_tree,
        bn = bn_tree,
        n_boot = config['build_tree']['n_bootstraps']
    shell:
        """
        CWD=$( pwd )
        DIR_TREE="$CWD"/{params.w:q}
        raxmlHPC-PTHREADS-AVX \
            -s {input:q} \
            -w "$DIR_TREE" \
            -m {params.model} \
            -T {threads} \
            -n {params.bn:q} \
            -f a \
            -x 42 \
            -p 42 \
            -# {params.n_boot} \
            2> {log:q}
        cat "Done" > {output:q}
        """