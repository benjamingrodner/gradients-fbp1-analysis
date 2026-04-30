import glob

rule get_crystal_struct_seqs:
    output:
        expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
    params:
        ids = config['rcsb_ids'],
        script = config['dir_scripts'] + "/download_rcsb_fastas.py",
        dir_out = config['dir_rcsb']
    log:
        "logs/get_crystal_struct_seqs.log"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python {params.script:q} \
            {params.ids} \
            -o {params.dir_out:q} \
        2> {log:q}
        """


rule subset_db_clusts:
    input:
        fn_db_rep_seqs,
    output:
        fn_db_rep_seqs_sub,
    params:
        pct = float(config['subset_db_clusts']['pct']) / 100,
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit sample -p {params.pct} {input:q} > {output:q}
        """




rule merge_clusters_with_crystal_struct_seqs:
    input:
        fn_db_rep_seqs_sub = fn_db_rep_seqs_sub,
        fns_env_rep_seqs = get_env_rep_seqs,
        fns_crystal_seqs = expand(fmt_crystal_seqs, rcsb_id=config['rcsb_ids']),
        fns_manual_ref = glob.glob(config['dir_ref_man'] + '/*'),
    output: 
        fn_db_crystal_seqs,
    shell:
        """
        cat {input.fn_db_rep_seqs_sub:q} > {output:q}
        cat {input.fns_env_rep_seqs:q} >> {output:q}
        cat {input.fns_crystal_seqs:q} >> {output:q}
        cat {input.fns_manual_ref:q} >> {output:q}
        """ 


rule alignment:
    input:
        fn_db_crystal_seqs,
    output:
        fn_alignment,
    log:
        "logs/alignment.log"
    benchmark:
        "benchmarks/alignment.benchmark.txt"
    threads:
        config['alignment']['threads'],
    resources:
        mem_mb=config['alignment']['mem_mb'],
        runtime=config['alignment']['runtime'],
    conda:
        "../envs/python.yaml"
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
        fn_trim_crystal = fn_trim_crystal,
        fn_trim_crystal_startend = fn_trim_crystal_startend,
    log:
        "logs/trim_alignment_based_on_crystal_seqs.log"
    threads: 4
    conda:
        "../envs/seqkit.yaml"
    params:
        formatted_ids = " -p ".join(config['rcsb_ids']),
    shell:
        r"""        
        # Get start and end positions
        seqkit grep -j {threads} -r -p {params.formatted_ids} {input:q}  \
            | seqkit replace -p '\s.+' -r '' \
            | seqkit fx2tab \
            | awk '{{ \
                match($2, /[^-]/); \
                s_pos=RSTART; \
                match($2, /.*[^-]/); \
                e_pos=RLENGTH; \
                if(min=="" || s_pos < min) min=s_pos; \
                if(e_pos > max) max=e_pos \
            }} END{{print min "," max}}' \
            > {output.fn_trim_crystal_startend:q} \
            2> {log:q}

        # Trim alignment file
        cat {output.fn_trim_crystal_startend:q} \
            | (IFS=',' read -r start end; 
                seqkit subseq -r "$start":"$end" {input:q} \
                > {output.fn_trim_crystal:q} \
                2>> {log:q}
            )
        """

rule trim_alignment_gappy_columns:
    input:
        fn_trim_crystal,
    output:
        fn_trim_clip,
    log:
        "logs/trim_alignment_gappy_columns.log"
    threads:
        config['trim_alignment_gappy_columns']['threads'],
    resources:
        mem_mb=config['trim_alignment_gappy_columns']['mem_mb'],
        runtime=config['trim_alignment_gappy_columns']['runtime'],
    conda:
        "../envs/clipkit.yaml"
    params:
        gaps = config['trim_alignment_gappy_columns']['gaps'],
        mode = config['trim_alignment_gappy_columns']['mode'],
    shell:
        """
        clipkit \
            {input:q} \
            -o {output:q} \
            -m {params.mode} \
            --gaps {params.gaps} \
            -t {threads} \
            -l \
            2> {log:q}
        """

rule filter_alignment:
    input:
        fn_trim_clip,
    output:
        fn_trim_clip_filt,
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

rule deduplicate_alignment:
    input:
        fasta = fn_trim_clip_filt
    output:
        fasta = fn_trim_clip_filt_dedup,
        mapping = fn_trim_clip_filt_dedup_map
    log:
        err = "logs/deduplicate_alignment.log"
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

