rule cluster_experiment_counts:
    input:
        wait = fn_quant_done,
        count = fmt_exp_counts_merge,
        clust = get_exp_clusters(fmt_exp_clusters),
    output: 
        fmt_exp_counts_clust,
    log:
        "logs/cluster_experiment_counts/{exp}.log"
    benchmark:
        "benchmarks/cluster_experiment_counts/{exp}.benchmark.txt"
    threads:
        config['cluster_counts']['threads'],
    resources:
        mem_mb=config['cluster_counts']['mem_mb'],
        runtime=config['cluster_counts']['runtime'],
    conda:
        "../envs/python.yaml"
    params:
        script=config['dir_scripts'] + "/cluster_experiment_counts.py",
        colname = lambda w: get_colname_contigs(w.exp),
    shell:
        """
        python3 {params.script:q} \
            -c {input.count:q} \
            -l {input.clust:q} \
            -n {params.colname} \
            -o {output:q} \
            -m {resources.mem_mb}M \
            -t {threads} \
            2> {log:q}
        """

rule get_experiment_metadata:
    input:
        wait = fn_read_and_count_downloads_done,
    output: 
        fmt_exp_meta,
    log:
        "logs/get_experiment_metadata/{exp}.log"
    conda:
        "../envs/python.yaml"
    params:
        script=config['dir_scripts'] + "/get_experiment_metadata.py",
        fn_exp_info = config['fn_experiment_info'],
        fn_sample_info = lambda w: get_fn_sample_info(w.exp),
    shell:
        """
        ARG_S="-s {params.fn_sample_info} "
        if [[ {params.fn_sample_info} == "None" ]]; do
            ARG_S=""
        fi
        python3 {params.script:q} \
            -g {params.fn_exp_info:q} \
            -e {wildcards.exp} \
            -o {output:q} \
            $ARG_S\
            2> {log:q}
        """

rule run_exp_deseq:
    input:
        counts = fmt_exp_counts_clust,
        meta = fmt_exp_meta,
    output: 
        dir_exp_deseq,
    log:
        "logs/run_exp_deseq/{exp}.log"
    benchmark:
        "benchmarks/run_exp_deseq/{exp}.benchmark.txt"
    threads:
        config['deseq']['threads'],
    resources:
        mem_mb=config['deseq']['mem_mb'],
        runtime=config['deseq']['runtime'],
    conda:
        "../envs/deseq.yaml"
    params:
        script=config['dir_scripts'] + "/run_exp_deseq.py",
        fn_exp_info = config['fn_experiment_info'],
        colname = lambda w: get_colname_contigs(w.exp),
        min_mean_counts = config['deseq']['min_sum_counts']
    shell:
        """
        python3 {params.script:q} \
            {input.counts:q} \
            {input.meta:q} \
            {params.fn_exp_info:q} \
            {wildcards.exp} \
            {output:q} \
            {threads:q} \
            2> {log:q}
        """

rule plot_exp_foldchange:
    input:
        counts = fmt_exp_counts_clust,
        meta = fmt_exp_meta,
        dir_stats = dir_exp_deseq,
    output:
        fmt_exp_plot,
    conda:
        "../envs/deseq.yaml"
    params:
        script=config['dir_scripts'] + "/plot_exp_foldchange.py",
        fn_config = config['fn_config'],
    shell:
        """
        python3 {params.script:q} \
            2> {log:q}
        """


rule deseq_done:
    input:
        expand(fmt_exp_plot, exp=DICT_EXP.keys())
    output:
        fn_deseq_done,
    shell:
        "echo {input:q} > {output:q}"