import os

rule get_db_seqs_to_cluster:
    input:
        fn_hmm_hitnames_all = fn_hmm_hitnames_all,
        fn_db_seqs = config['path_database'],
    output:
        fn_db_seqs_to_cluster,
    log:
        "logs/get_db_seqs_to_cluster.log"
    benchmark:
        "benchmarks/get_db_seqs_to_cluster.benchmark.txt"
    threads: 4
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -j {threads} \
            -f {input.fn_hmm_hitnames_all:q} \
            {input.fn_db_seqs:q} \
            > {output:q} \
            2> {log:q}
        """


rule cluster_db_seqs:
    input:
        fn_db_seqs_to_cluster,
    output:
        fn_db_rep_seqs = fn_db_rep_seqs,
        fn_db_clusters = fn_db_clusters,
    params:
        bn_cluster = lambda wildcards, input, output: output.fn_db_rep_seqs.replace('_rep_seq.fasta',''),
        dir_tmp = lambda wildcards, input, output: os.path.dirname(output.fn_db_rep_seqs),
        min_seq_id = config['cluster_db_seqs']['min_seq_id'],
        coverage = config['cluster_db_seqs']['coverage'],
        cov_mode = config['cluster_db_seqs']['cov_mode'],
    log:
        "logs/cluster_db_seqs.log"
    benchmark:
        "benchmarks/cluster_db_seqs.benchmark.txt"
    threads:
        config['cluster_db_seqs']['threads'],
    resources:
        mem_mb=config['cluster_db_seqs']['mem_mb'],
        runtime=config['cluster_db_seqs']['runtime'],
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs easy-cluster \
            {input:q} \
            {params.bn_cluster:q} \
            {params.dir_tmp:q} \
            --min-seq-id {params.min_seq_id} \
            -c {params.coverage} \
            --cov-mode {params.cov_mode} \
            2> {log:q}
        """

rule group_env_hitnames:
    input:
        fn_hmms_best_hit,
        fmt_bigtable = fmt_bigtable,
    output:
        fmt_grouping_done,
    log:
        "logs/group_env_hitnames.log"
    benchmark:
        "benchmarks/group_env_hitnames.benchmark.txt"
    threads:
        config['group_env_hitnames']['threads'],
    resources:
        mem_mb=config['group_env_hitnames']['mem_mb'],
        runtime=config['group_env_hitnames']['runtime'],
    conda:
        "../envs/python.yaml"
    params:
        script=config['dir_scripts'] + "/group_env_hitnames.py",
    shell:
        """
        DIR_OUT=$( basename {output:q} )
        python {params.script:q} \
            -bt {input.fmt_bigtable:q} \
            -s {input.fmt_seqs_parquet:q} \
            -d "$DIR_OUT" \
            -o {output:q} \
            -t {threads} \
            -m {resources.mem_mb}
        """
