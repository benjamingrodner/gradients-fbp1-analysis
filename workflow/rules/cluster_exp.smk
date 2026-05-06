rule get_exp_seqs_to_cluster:
    input:
        fn_hitnames = fmt_exp_hitnames_grouped_gene,
        fn_seqs = lambda w: get_fn_exp_assembly(w.exp),
    output:
        fmt_exp_seqs_to_cluster,
    log:
        "logs/get_exp_seqs_to_cluster.log"
    benchmark:
        "benchmarks/get_exp_seqs_to_cluster.benchmark.txt"
    threads: 4
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit grep -j {threads} \
            -f {input.fn_hitnames:q} \
            {input.fn_seqs:q} \
            > {output:q} \
            2> {log:q}
        """

rule cluster_exp_seqs:
    input:
        fmt_exp_seqs_to_cluster,
    output:
        fn_rep = fmt_exp_rep_seqs,
        fn_clust = fmt_exp_clusters,
    params:
        bn_cluster = lambda wildcards, input, output: output.fn_rep.replace('_rep_seq.fasta',''),
        dir_tmp = lambda wildcards, input, output: os.path.dirname(output.fn_rep),
        min_seq_id = config['cluster_exp_seqs']['min_seq_id'],
        coverage = config['cluster_exp_seqs']['coverage'],
        cov_mode = config['cluster_exp_seqs']['cov_mode'],
    log:
        "logs/cluster_exp_seqs.log"
    benchmark:
        "benchmarks/cluster_exp_seqs.benchmark.txt"
    threads:
        config['cluster_exp_seqs']['threads'],
    resources:
        mem_mb=config['cluster_exp_seqs']['mem_mb'],
        runtime=config['cluster_exp_seqs']['runtime'],
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
