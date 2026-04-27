import os
import shutil

checkpoint group_env_hitnames:
    input:
        fn_hmms_best_hit = fn_hmms_best_hit,
        fmt_bigtable = fmt_bigtable,
    output:
        temp(directory(fmt_env_hitnames_dir_tmp)),
    log:
        "logs/group_env_hitnames/{batch}.log"
    benchmark:
        "benchmarks/group_env_hitnames/{batch}.benchmark.txt"
    threads:
        config['group_env_stuff']['threads'],
    resources:
        mem_mb=config['group_env_stuff']['mem_mb'],
        runtime=config['group_env_stuff']['runtime'],
    conda:
        "../envs/python.yaml"
    params:
        script=config['dir_scripts'] + "/group_env_hitnames.py",
        ttaxnames = config['ttaxnames'],
        wc_prefix = tg_prefix,
    shell:
        """
        python {params.script:q} \
            --table1 {input.fn_hmms_best_hit:q} \
            --parquet {input.fmt_bigtable:q} \
            --ttaxnames {params.ttaxnames} \
            --threads {threads} \
            --memory {resources.mem_mb}MB \
            --b_string {wildcards.batch} \
            --dir_out {output:q} \
            --wc_prefix {params.wc_prefix} \
            2> {log:q}
        """

rule check_fns_hitnames:
    output:
        fmt_env_hitnames,
    params:
        fn_env_hitnames_tmp = fmt_env_hitnames_tmp,
    run:
        fn_in = params.fn_env_hitnames_tmp.format(
            batch=wildcards.batch, taxgene=wildcards.taxgene
        )
        if os.path.exists(fn_in):
            with open(fn_in, 'r') as fr, open(output[0], 'w') as fo:
                shutil.copyfileobj(fr, fo)
        else:
            Path(output[0]).touch()


rule get_env_hitseqs:
    input:
        fn_env_hitnames = fmt_env_hitnames,
        fn_seqs_parquet = fmt_seqs_parquet,
        # fn_env_hitseqs_dir = fmt_env_hitseqs_dir
    output:
        fmt_env_hitseqs,
    log:
        "logs/get_env_hitseqs/{batch}/{taxgene}.log"
    benchmark:
        "benchmarks/get_env_hitseqs/{batch}/{taxgene}.benchmark.txt"
    threads:
        config['group_env_stuff']['threads'],
    resources:
        mem_mb=config['group_env_stuff']['mem_mb'],
        runtime=config['group_env_stuff']['runtime'],
    conda:
        "../envs/python.yaml"
    params:
        script=config['dir_scripts'] + "/parquet2fasta.py",
    shell:
        """
        python {params.script:q} \
            -hc contig_name \
            -hf {input.fn_env_hitnames:q} \
            -i {input.fn_seqs_parquet:q} \
            -o {output:q} \
            -m {resources.mem_mb}MB \
            -t {threads} \
            2>> {log:q}
        """


rule group_env_hitseqs_batch:
    input:
        fns_env_hitseqs = lambda wildcards: [
            fmt_env_hitseqs.format(taxgene=wildcards.taxgene, batch=batch)
            for batch in BATCHES
        ],
        # fn_env_hitseqs_dir = expand(fmt_env_hitseqs_dir, batch=BATCHES)
    output:
        fmt_grouped_env_hitseqs,
    shell:
        """
        cat {input.fns_env_hitseqs:q} > {output:q}
        """

rule cluster_env_hitseqs:
    input:
        fmt_grouped_env_hitseqs,    
    output:
        fmt_clustered_env_hitseqs,
    params:
        min_seq_id = config['cluster_env_hitseqs']['min_seq_id'],
        coverage = config['cluster_env_hitseqs']['coverage'],
        cov_mode = config['cluster_env_hitseqs']['cov_mode'],
    log:
        "logs/cluster_env_hitseqs/{taxgene}.log"
    benchmark:
        "benchmarks/cluster_env_hitseqs/{taxgene}.benchmark.txt"
    threads:
        config['cluster_env_hitseqs']['threads'],
    resources:
        mem_mb=config['cluster_env_hitseqs']['mem_mb'],
        runtime=config['cluster_env_hitseqs']['runtime'],
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        if [[ -s {input:q} ]]; then
            bn_cluster={output:q}
            bn_cluster="${{bn_cluster%_rep_seq.fasta}}"
            dir_cluster=$( dirname {output:q} )
            mmseqs easy-cluster \
                {input:q} \
                "$bn_cluster" \
                "$dir_cluster"/tmp \
                --min-seq-id {params.min_seq_id} \
                -c {params.coverage} \
                --cov-mode {params.cov_mode} \
                2>> {log:q}
        else
            touch {output:q}
        fi
        """


rule merge_env_clusters:
    input:
        aggregate_env_clusters,
    output:
        fn_env_seqs_clust_cat,
    shell:
        """
        cat {input:q} > {output:q}
        """
    
    