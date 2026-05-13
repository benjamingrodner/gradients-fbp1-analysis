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
        "../envs/salmon.yml"
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

rule exp_quant_done:
    input:
        gz_done = fn_agg_biosample_dirs,
        fns = aggregate_exp_salmon_counts,
    output:
        fmt_exp_salmon_quant_done,
    shell:
        """
        echo {input.fns:q} > {output:q}
        """


rule quant_done:
    input:
        expand_exp_quant(fmt_exp_salmon_quant_done)
    output:
        fn_salmon_quant_done,
    shell:
        """
        echo {input:q} > {output:q}
        """
