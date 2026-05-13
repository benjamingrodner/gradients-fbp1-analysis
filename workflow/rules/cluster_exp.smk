checkpoint group_exp_hitnames:
    input:
        fmt_hmms_best_hit_exp,
    output:
        directory(dir_exp_hitnames_grouped_gene)
    log:
        "logs/group_exp_hitnames/{exp}.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + "/group_exp_hitnames.py",
        fmt_source_files = lambda w: fmt_table_hmm_exp.format(
            exp=w.exp, gene="{gene}"
        ),
        fmt_out = lambda w: fmt_exp_hitnames_grouped_gene.format(
            exp=w.exp, gene="{gene}"
        ),
    shell:
        """
        python {params.script:q} \
            --fn_best_hit {input:q} \
            --fmt_source_files {params.fmt_source_files:q} \
            --fmt_out {params.fmt_out} \
            2> {log:q}
        """


rule get_exp_seqs_to_cluster:
    input:
        fn_hitnames = fmt_exp_hitnames_grouped_gene,
        fn_seqs = fmt_exp_assembly_6tr_rename,
    output:
        fmt_exp_seqs_to_cluster,
    log:
        "logs/get_exp_seqs_to_cluster/{exp}_{gene}.log"
    benchmark:
        "benchmarks/get_exp_seqs_to_cluster/{exp}_{gene}.benchmark.txt"
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

rule cluster_exp_hitseqs:
    input:
        fmt_exp_seqs_to_cluster,
    output:
        fn_rep = fmt_exp_rep_seqs,
        fn_clust = fmt_exp_clusters,
    params:
        bn_cluster = lambda wildcards, input, output: output.fn_rep.replace('_rep_seq.fasta',''),
        dir_tmp = lambda wildcards, input, output: os.path.dirname(output.fn_rep),
        min_seq_id = config['cluster_exp_hitseqs']['min_seq_id'],
        coverage = config['cluster_exp_hitseqs']['coverage'],
        cov_mode = config['cluster_exp_hitseqs']['cov_mode'],
    log:
        "logs/cluster_exp_hitseqs/{exp}_{gene}.log"
    benchmark:
        "benchmarks/cluster_exp_hitseqs/{exp}_{gene}.benchmark.txt"
    threads:
        config['cluster_exp_hitseqs']['threads'],
    resources:
        mem_mb=config['cluster_exp_hitseqs']['mem_mb'],
        runtime=config['cluster_exp_hitseqs']['runtime'],
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
