checkpoint merge_biosample_fastqs:
    input:
        d = lambda w: get_dir_fastq(w.exp_quant),
        srr_info = lambda w: get_srr_info(w.exp_quant),
    output:
        temp(directory(dir_biosamples_fastq_merged)),
    log:
        "logs/merge_biosample_fastqs/{exp_quant}.log",
    threads: config['get_bioproject_reads']['threads'],
    params:
        merge = lambda w: check_merge_biosample(w.exp_quant)
    shell:
        """
        mkdir -p {output:q}
        if [[ {params.merge} == "no" ]]; then
            find {input.d:q} -type f -print0 \
                | xargs -0 -P {threads} -I {{}} cp {{}} {output:q} \
                2> {log:q}
        else
            > {log:q}
            {{
            ## Clear output files
            while IFS=',' read -r srr _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ sample _ _ _ _ _ _ _ _; do
                if [[ "$srr" != "Run" ]]; then
                    > {output:q}/"$sample"_1.fastq 2>> {log:q}
                    > {output:q}/"$sample"_2.fastq 2>> {log:q}
                fi
            done
            }} < {input.srr_info:q}
            
            ## Merge the SRR fastqs
            {{
            while IFS=',' read -r srr _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ sample _ _ _ _ _ _ _ _; do
                if [[ "$srr" != "Run" ]]; then
                    zcat {input.d:q}/"$srr"_1.fastq.gz >> {output:q}/"$sample"_1.fastq 2>> {log:q}
                    zcat {input.d:q}/"$srr"_2.fastq.gz >> {output:q}/"$sample"_2.fastq 2>> {log:q}
                fi
            done
            }} < {input.srr_info:q}
        fi
        """


rule gzip_biosample_fastqs:
    input:
        dir_biosamples_fastq_merged,
    output:
        fmt_gzip_biosample_fastq_done,
    log:
        "logs/gzip_biosample_fastqs/{exp_quant}.log",
    threads: config['get_bioproject_reads']['threads'],
    conda:
        "../envs/sra_tools.yaml"
    params:
        merge = lambda w: check_merge_biosample(w.exp_quant)
    shell:
        """
        if [[ {params.merge} == "yes" ]]; then
            pigz -p {threads} {input:q}/*.fastq 2> {log:q}
        fi
        echo {input:q} > {output:q}
        """


rule trim_fastq:
    input:
        d = dir_biosamples_fastq_merged,
        gz_done = fmt_gzip_biosample_fastq_done,
        r1 = fmt_biosample_fastq_merged_r1,
        r2 = fmt_biosample_fastq_merged_r2,
    output:
        temp(fmt_fastq_trim_r1),
        temp(fmt_fastq_trim_r1_unpaired),
        temp(fmt_fastq_trim_r2),
        temp(fmt_fastq_trim_r2_unpaired),
    params:
        opts=lambda w: get_trim_opts(w.exp_quant),
    log:
        "logs/trim/{exp_quant}/{sample}.trim.log",
    benchmark:
        "benchmarks/trim/{exp_quant}/{sample}.benchmark.txt"
    conda:
        "../envs/trimmomatic.yaml"
    threads: config["trim"]["threads"]
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.r1:q} {input.r2:q} \
            {output:q} \
            {params.opts} \
            > {log:q} 2>&1
        """

# rule gzip_trim:
#     input:
#         temp(fmt_fastq_trim_r1_unz),
#         temp(fmt_fastq_trim_r1_unpaired_unz),
#         temp(fmt_fastq_trim_r2_unz),
#         temp(fmt_fastq_trim_r2_unpaired_unz),
#     output:
#         temp(fmt_fastq_trim_r1),
#         temp(fmt_fastq_trim_r1_unpaired),
#         temp(fmt_fastq_trim_r2),
#         temp(fmt_fastq_trim_r2_unpaired),
#     threads: config["trim"]["threads"]
#     shell:
#         """
#         pigz 
#         """


rule fastqc_raw:
    input:
        d = dir_biosamples_fastq_merged,
        gz_done = fmt_gzip_biosample_fastq_done,
        fn = fmt_biosample_fastq_merged,
    output:
        fmt_fastqc_raw,
    log:
        "logs/fastqc_raw/{exp_quant}/{sample}_{read}.log",
    conda:
        "../envs/fastqc.yaml"
    params:
        out_dir=lambda wildcards, input, output: Path(output[0]).parent,
    shell:
        """
        fastqc -o {params.out_dir:q} {input.fn:q} > {log:q} 2>&1
        """

rule fastqc_trimmed:
    input:
        fmt_fastq_trim,
    output:
        fmt_fastqc_trimmed,
    log:
        "logs/fastqc_trimmed/{exp_quant}/{sample}_{read}.log",
    conda:
        "../envs/fastqc.yaml"
    params:
        out_dir=lambda wildcards, input, output: Path(output[0]).parent,
    shell:
        """
        fastqc -o {params.out_dir:q} {input:q} > {log:q} 2>&1
        """

rule multiqc_raw:
    input:
        aggregate_fastqc_raw,
    output:
        report=fmt_multiqc_raw,
    params:
        in_dir=lambda wildcards, input, output: Path(input[0]).parent,
        out_dir=lambda wildcards, input, output: Path(output[0]).parent,
    log:
        "logs/multiqc_raw/{exp_quant}.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc -o {params.out_dir:q} {params.in_dir:q} > {log:q} 2>&1
        """

rule multiqc_trimmed:
    input:
        aggregate_fastqc_trimmed,
    output:
        report=fmt_multiqc_trimmed,
    params:
        in_dir=lambda wildcards, input, output: Path(input[0]).parent,
        out_dir=lambda wildcards, input, output: Path(output[0]).parent,
    log:
        "logs/multiqc_trimmed/{exp_quant}.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc -o {params.out_dir:q} {params.in_dir:q} > {log:q} 2>&1
        """

rule read_prep_done:
    input:
        expand_exp_quant(fmt_multiqc_raw),
        expand_exp_quant(fmt_multiqc_trimmed),
    output:
        fn_read_prep_done,
    shell:
        """
        echo {input:q} > {output:q}
        """