# Copied from Chris Berthiaume's pipeline
rule salmon_idx:
    input:
        assm = fmt_exp_assembly_quant,
    output:
        idx=directory(dir_salmon_idx),
    params:
        out_dir=lambda wildcards, input, output: Path(output.idx).parent,
    log:
        "logs/salmon_idx/{exp_quant}.log",
    benchmark:
        "benchmarks/salmon_idx/{exp_quant}.benchmark.txt"
    threads:
        config["salmon"]["index"]["threads"]
    resources:
        mem_mb=config['salmon']["index"]['mem_mb'],
        runtime=config['salmon']["index"]['runtime'],
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        [[ ! -d {params.out_dir:q} ]] && mkdir -p {params.out_dir:q} 2> {log:q}
        salmon index \
            --threads {threads} \
            --kmerLen 31 \
            --transcripts {input.assm:q} \
            --index {output.idx:q} \
            >> {log:q} 2>&1
        """

# Copied from Chris Berthiaume's pipeline
rule salmon_counts:
    input:
        left = fmt_fastq_trim_r1,
        right = fmt_fastq_trim_r2,
        idx = dir_salmon_idx,
    output:
        quant=fmt_quant,
    params:
        out_dir=lambda wildcards, input, output: Path(output.quant).parent,
        quant_nogz=lambda wildcards, input, output: Path(output.quant).with_suffix(""),
    log:
        "logs/salmon/{exp_quant}/{sample}.log"
    benchmark:
        "benchmarks/salmon/{exp_quant}/{sample}.benchmark.txt"
    threads:
        config['salmon']['threads'],
    resources:
        mem_mb=config['salmon']['mem_mb'],
        runtime=config['salmon']['runtime'],
    conda:
        "../envs/salmon.yaml"
    shell:
        """
        salmon quant -i {input.idx:q} --libType A \
            -p {threads} --validateMappings \
            -1 {input.left:q} -2 {input.right:q} \
            -o {params.out_dir:q} >> {log:q} 2>&1
        
        pigz -p {threads} {params.quant_nogz:q} 2>> {log:q}        
        """

rule merge_salmon_counts:
    input:
        wait = fn_read_prep_done,
        fns = aggregate_exp_salmon_counts,
    output: 
        fmt_salmon_counts_merge,
    log:
        "logs/merge_salmon_counts/{exp_quant}.log"
    benchmark:
        "benchmarks/merge_salmon_counts/{exp_quant}.benchmark.txt"
    threads:
        config['salmon']['merge']['threads'],
    resources:
        mem_mb=config['salmon']['merge']['mem_mb'],
        runtime=config['salmon']['merge']['runtime'],
    conda:
        "../envs/python.yaml"
    params:
        script=config['dir_scripts'] + "/aggregate_counts.py",
        prefix = lambda w: DICT_EXP[w.exp_quant]['prefix'],
    shell:
        """
        python3 {params.script:q} \
            --jobs {threads} \
            --prefix {params.prefix} \
            {input.fns:q} {output:q} \
            2> {log:q}
        """

# rule salmon_quant_done:
#     input:
#         expand_exp_quant(fmt_salmon_counts_merge)
#     output:
#         fn_salmon_quant_done,
#     shell:
#         """
#         echo {input:q} > {output:q}
#         """

rule merge_counts_fromauthor_and_downloaded:
    input:
        download_done = fn_read_and_count_downloads_done,
    output:
        fmt_fromauthor_and_downloaded_counts_merge,
    log:
        "logs/merge_salmon_counts/{exp_auth}.log"
    benchmark:
        "benchmarks/merge_salmon_counts/{exp_auth}.benchmark.txt"
    threads:
        config['salmon']['merge']['threads'],
    conda:
        "../envs/python.yaml"
    params:
        script = lambda w: get_script_merge_counts_auth(w.exp_auth),
        colname_contigs = lambda w: get_colname_merge_counts_auth(w.exp_auth),
        glob_counts = lambda w: get_glob_counts_auth(w.exp_auth),
        prefix = lambda w: DICT_EXP[w.exp_auth]['prefix'],
    shell:
        """
        ARG_C="-c {params.colname_contigs}"
        if [[ -z "{params.colname_contigs}" ]]; then
            ARG_C=""
        fi
        FNS=$( ls {params.glob_counts} )
        python3 {params.script} \
            -j {threads} \
            -p {params.prefix} \
            $ARG_C \
            $FNS \
            {output:q} \
            2> {log:q}
        """

rule quant_done:
    input:
        expand_exp_quant(fmt_salmon_counts_merge),
        expand_exp_auth(fmt_fromauthor_and_downloaded_counts_merge),
    output:
        fn_quant_done,
    shell:
        "echo {input:q} > {output:q}"
