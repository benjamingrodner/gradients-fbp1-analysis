import os
import glob

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
    threads:
        config['get_db_seqs_to_cluster']['threads'],
    resources:
        mem_mb=config['get_db_seqs_to_cluster']['mem_mb'],
        runtime=config['get_db_seqs_to_cluster']['runtime'],
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        seqkit grep \
            -f {input.fn_hmm_hitnames_all:q} \
            {input.fn_db_seqs:q} \
            > {output.fn_db_seqs_to_cluster:q} \
            2> {log:q}
        """


rule cluster_db_seqs:
    input:
        fn_db_seqs_to_cluster,
    output:
        fn_db_rep_seqs,
        fn_db_clusters,
    params:
        bn_cluster = lambda output: output.fn_db_rep_seqs.replace('_rep_seq.fasta',''),
        dir_tmp = lambda output: os.path.dirname(output.fn_db_rep_seqs),
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
        "envs/mmseqs.yaml"
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


rule get_crystal_struct_seqs:
    output:
        expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
    params:
        ids = config['rcsb_ids'],
        script = "scripts/download_rcsb_fastas.py",
        dir_out = config['dir_rcsb']
    log:
        "logs/get_crystal_struct_seqs.log"
    conda:
        "envs/python.yaml"
    shell:
        """
        python {params.script:q} \
            {params.ids} \
            -o {params.dir_out:q} \
        2> {log:q}
        """


rule merge_clusters_with_crystal_struct_seqs:
    input:
        fn_db_rep_seqs = fn_db_rep_seqs,
        fns_crystal_seqs = glob.glob(config['dir_rcsb'] + '/*')
    output: 
        fn_db_crystal_seqs,
    shell:
        """
        cat {input.fn_db_rep_seqs:q} > {output:q}
        cat {input.fns_crystal_seqs:q} >> {output:q}
        """ 


rule align_db_crystal_seqs:
    input:
        fn_db_crystal_seqs,
    output:
        fn_alignment,
    log:
        "logs/align_db_crystal_seqs.log"
    benchmark:
        "benchmarks/align_db_crystal_seqs.benchmark.txt"
    threads:
        config['align_db_crystal_seqs']['threads'],
    resources:
        mem_mb=config['align_db_crystal_seqs']['mem_mb'],
        runtime=config['align_db_crystal_seqs']['runtime'],
    conda:
        "envs/python.yaml"
    shell:
        """
        mafft --auto \
            --thread {threads} \
            {input:q} \
            > {output:q} \
            2> {log:q}
        """


rule trim_alignment_based_on_crystal_seqs:
    input:
        fn_alignment,
    output:
        fn_trim_crystal,
        fn_trim_crystal_startend,
    log:
        "logs/trim_alignment_based_on_crystal_seqs.log"
    threads:
        config['get_db_seqs_to_cluster']['threads'],
    resources:
        mem_mb=config['get_db_seqs_to_cluster']['mem_mb'],
        runtime=config['get_db_seqs_to_cluster']['runtime'],
    conda:
        "envs/seqkit.yaml"
    params:
        formatted_ids = lambda input: " -p ".join(config['rcsb_ids_to_trim_on']),
    shell:
        """        
        # Get start and end positions
        seqkit grep -r -p {params.formatted_ids} {input:q}  \
            | seqkit replace -p "\s.+" -r "" \
            | seqkit fx2tab \
            | awk '{ 
                match($2, /[^-]/); 
                s_pos=RSTART; 
                match($2, /.*[^-]/); 
                e_pos=RLENGTH; 
                if(min=="" || s_pos < min) min=s_pos; 
                if(e_pos > max) max=e_pos
            } END{print min "," max}' \
            > {output.fn_trim_crystal_startend:q} \
            2> {log:q}

        # Trim alignment file
        cat {output.fn_trim_crystal_startend:q} \
            | (IFS=',' read -r start end; 
                seqkit subseq -r "$start":"$end" "$FN_ALN" \
                > {output.fn_trim_crystal:q} \
                2>> {log:q}
            )
        """

rule trim_alignment_gappy_columns:


rule pick_aa_substitution_model:


rule build_tree: