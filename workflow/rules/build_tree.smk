
rule build_ml_tree:
    input:
        fn_trim_clip_filt,
    output:
        fn_ml_tree,
    log:
        "logs/build_ml_tree.log"
    benchmark:
        "benchmarks/build_ml_tree.benchmark.txt"
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    threads:
        config['build_tree']['threads'],
    shell:
        """
        DIR_TREE=$( dirname {output:q} )
        CWD=$( pwd )
        DIR_TREE="$CWD"/"$DIR_TREE"
        BN=$( basename {output:q} )
        BN="${{BN#RAxML_bestTree.}}"
        raxmlHPC-PTHREADS-AVX \
            -s {input:q} \
            -w "$DIR_TREE" \
            -m PROTGAMMAAUTO \
            -n "$BN" \
            -p 42 \
            -T {threads} \
            2> {log:q}
        """


rule get_tree_bootstraps:
    input:
        fn_trim_clip_filt,
    output:
        fn_bootstraps_done,
    log:
        "logs/get_tree_bootstraps.log"
    benchmark:
        "benchmarks/get_tree_bootstraps.benchmark.txt"
    threads:
        config['build_tree']['threads'],
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    params:
        model = config['build_tree']['model']
    shell:
        """
        DIR_TREE=$( dirname {output:q} )
        CWD=$( pwd )
        DIR_TREE="$CWD"/"$DIR_TREE"
        BN=$( basename {output:q} )
        BN="${{BN#RAxML_bootstrap.}}"
        BN="${{BN%.*}}"
        raxmlHPC-PTHREADS-AVX \
            -x 42 \
            -p 42 \
            -# autoMRE \
            -m {params.model} \
            -s {input:q} \
            -w "$DIR_TREE" \
            -n "$BN" \
            -T {threads} \
            2> {log:q}
        echo "Done" > {output:q}
        """

rule bootstraps_to_ml_tree:
    input:
        fn_ml_tree = fn_ml_tree,
        fn_bootstraps_done = fn_bootstraps_done,
    output:
        fn_ml_bootstrap,
    log:
        "logs/bootstraps_to_ml_tree.log"
    benchmark:
        "benchmarks/bootstraps_to_ml_tree.benchmark.txt"
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    params:
        model = config['build_tree']['model']
    shell:
        """
        DIR_TREE=$( dirname {output:q} )
        CWD=$( pwd )
        DIR_TREE="$CWD"/"$DIR_TREE"
        BN=$( basename {output:q} )
        BN="${{BN#RAxML_bipartitions.}}"
        FN_B={input.fn_bootstraps_done:q}
        FN_B="${{FN_B%.*}}"
        raxmlHPC-PTHREADS-AVX \
            -f b \
            -m {params.model} \
            -z "$FN_B" \
            -t {input.fn_ml_tree:q} \
            -w "$DIR_TREE" \
            -n "$BN" \
            2> {log:q}
        """

rule majority_rule_consensus_tree:
    input:
        fn_bootstraps_done = fn_bootstraps_done,
    output:
        fn_mrc,
    log:
        "logs/majority_rule_consensus_tree.log"
    benchmark:
        "benchmarks/majority_rule_consensus_tree.benchmark.txt"
    resources:
        mem_mb=config['build_tree']['mem_mb'],
        runtime=config['build_tree']['runtime'],
    params:
        model = config['build_tree']['model']
    shell:
        """
        DIR_TREE=$( dirname {output:q} )
        CWD=$( pwd )
        DIR_TREE="$CWD"/"$DIR_TREE"
        BN=$( basename {output:q} )
        BN="${{BN#RAxML_result.}}"
        FN_B={input.fn_bootstraps_done:q}
        FN_B="${{FN_B%.*}}"
        raxmlHPC-PTHREADS-AVX \
            -J MR \
            -m {params.model} \
            -z "$FN_B" \
            -w "$DIR_TREE" \
            -n "$BN" \
            2> {log:q}
        """


rule full_tree_analysis:
    input:
        fn_trim_clip_filt,
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
        model = config['build_tree']['model']
    shell:
        """
        DIR_TREE=$( dirname {output:q} )
        CWD=$( pwd )
        DIR_TREE="$CWD"/"$DIR_TREE"
        BN=$( basename {output:q} )
        BN="${{BN%.*}}"
        raxmlHPC-PTHREADS-AVX \
            -s {input:q} \
            -w "$DIR_TREE" \
            -m {params.model} \
            -T {threads} \
            -n "$BN" \
            -f a \
            -x 42 \
            -p 42 \
            -# autoMRE \
            2> {log:q}
        cat "Done" > {output:q}
        """