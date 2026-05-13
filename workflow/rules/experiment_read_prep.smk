checkpoint merge_biosample_fastqs:
    input:
        d = lambda w: get_dir_fastq(w.exp_quant),
        srr_info = lambda w: get_srr_info(w.exp_quant),
    output:
        temp(directory(dir_biosamples_fastq_merged)),
    log:
        "logs/merge_biosample_fastqs/{exp_quant}.log",
    shell:
        """
        ## Clear output files
        mkdir -p {output:q}
        {{
        while IFS=',' read -r srr _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ sample _ _ _ _ _ _ _ _; do
            if [[ "$srr" != "Run" ]]; then
                > {output:q}/"$sample"_1.fastq 2> {log:q}
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
        """
checkpoint gzip_biosample_fastqs:
    input:
        expand_exp_quant(dir_biosamples_fastq_merged),
    output:
        fn_agg_biosample_dirs,
    log:
        "logs/gzip_biosample_fastqs.log",
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        > {log:q}
        for d in {input:q}; do
            pigz "$d"/*.fastq 2>> {log:q}
        done
        echo {input:q} > {output:q}
        """



rule trim_fastq:
    input:
        gz_done = fn_agg_biosample_dirs,
        r1 = fmt_biosample_fastq_merged_r1,
        r2 = fmt_biosample_fastq_merged_r2,
    output:
        temp(fmt_fastq_trim_r1),
        temp(fmt_fastq_trim_r1_unpaired),
        temp(fmt_fastq_trim_r2),
        temp(fmt_fastq_trim_r2_unpaired),
    params:
        adapter_file=lambda w: get_adapter_file(w.exp_quant),
    log:
        "logs/trim/{exp_quant}/{sample}.trim.log",
    benchmark:
        "benchmarks/trim/{exp_quant}/{sample}.benchmark.txt"
    conda:
        "../envs/trimmomatic.yml"
    threads: config["trim"]["threads"]
    shell:
        """
        CLIP=""
        if [[ -n {params.adapter_file:q} ]]; then
            CLIP="ILLUMINACLIP:{params.adapter_file:q}:2:30:10:1:true"
        fi

        trimmomatic PE -threads {threads} {input:q} {output:q} \
            "$CLIP" MAXINFO:135:0.5 \
            LEADING:3 TRAILING:3 MINLEN:60 AVGQUAL:20 > {log:q} 2>&1
        """

rule fastqc_raw:
    input:
        gz_done = fn_agg_biosample_dirs,
        fn = fmt_biosample_fastq_merged,
    output:
        fmt_fastqc_raw,
    log:
        "logs/fastqc_raw/{exp_quant}/{sample}_{read}.log",
    conda:
        "../envs/fastqc.yml"
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
        "../envs/fastqc.yml"
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
        "../envs/multiqc.yml"
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
        "../envs/multiqc.yml"
    shell:
        """
        multiqc -o {params.out_dir:q} {params.in_dir:q} > {log:q} 2>&1
        """

