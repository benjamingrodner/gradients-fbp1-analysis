rule download_counts:
    output:
        d = directory(dir_download_counts),
        done = fmt_download_counts_done,
    log:
        "logs/download_counts/{exp_download_counts}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        link = lambda w: DICT_EXP[w.exp_download_counts]['link_counts']
    shell:
        """
        wget {params.link:q} -P {output.d:q} \
            2> {log:q}
        echo "Done" > {output.done:q}
        """

rule get_bioproject_info:
    output:
        bp_info = fmt_bioproject_info,
        srr_info = fmt_bioproject_srr_info,
        srr_list = fmt_bioproject_srr_list,
    log:
        "logs/get_bioproject_info/{exp_bioproject}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        bioproject = lambda w: DICT_EXP[w.exp_bioproject]['bioproject'],
        elink_target = lambda w: get_elink_target(w.exp_bioproject)
    shell:
        """
        esearch -db bioproject -query {params.bioproject} \
            | elink -target {params.elink_target} \
            | efetch \
            > {output.bp_info:q} \
            2> {log:q}

        esearch -db sra -query {params.bioproject} \
            | efetch -format runinfo \
            > {output.srr_info:q} \
            2>> {log:q}

        cat {output.srr_info:q} \
            | cut -d ',' -f 1 \
            | grep SRR \
            > {output.srr_list:q} \
            2>> {log:q}
        """

rule get_bioproject_reads:
    input:
        fmt_bioproject_srr_list,
    output:
        dir_out = temp(directory(dir_bioproject_fastq)),
        fn_done = fmt_bioproject_fastq_done,
    log:
        "logs/get_bioproject_reads/{exp_bioproject}.log"
    benchmark:
        "benchmarks/get_bioproject_reads/{exp_bioproject}.benchmark.txt"
    threads:
        config['get_bioproject_reads']['threads'],
    resources:
        mem_mb=config['get_bioproject_reads']['mem_mb'],
        runtime=config['get_bioproject_reads']['runtime'],
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        prefetch \
            --option-file {input:q} \
            --output-directory {output.dir_out:q} \
            2> {log:q}
        for file in {output.dir_out:q}/*/*.sra; do \
            fasterq-dump \
                --split-files \
                --outdir {output.dir_out:q} \
                --temp {output.dir_out:q}/tmp \
                --threads {threads} \
                --mem {resources.mem_mb}MB \
                "$file" \
                2>> {log:q}
        done
        pigz -p {threads} {output.dir_out:q}/*.fastq
        echo "Done" > {output.fn_done:q}
        """

rule get_biosample_info:
    output:
        bp_info = fmt_biosample_info,
        srr_info = fmt_biosample_srr_info,
        srr_list = fmt_biosample_srr_list,
    log:
        "logs/get_biosample_info/{exp_biosample}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        biosamples = lambda w: DICT_EXP[w.exp_biosample]['biosamples']
    shell:
        """
        > {output.bp_info:q}
        > {output.srr_info:q}
        > {output.srr_list:q}
        for bs in {params.biosamples}; do
            esearch -db biosample -query "$bs" \
                | efetch \
                >> {output.bp_info:q} \
                2> {log:q}

            esearch -db sra -query "$bs" \
                | efetch -format runinfo \
                >> {output.srr_info:q} \
                2>> {log:q}
        done
        
        cat {output.srr_info:q} \
            | cut -d ',' -f 1 \
            | grep SRR \
            >> {output.srr_list:q} \
            2>> {log:q}
        """

rule get_biosample_reads:
    input:
        fmt_biosample_srr_list,
    output:
        dir_out = temp(directory(dir_biosample_fastq)),
        fn_done = fmt_biosample_fastq_done,
    log:
        "logs/get_biosample_reads/{exp_biosample}.log"
    benchmark:
        "benchmarks/get_biosample_reads/{exp_biosample}.benchmark.txt"
    threads:
        config['get_bioproject_reads']['threads'],
    resources:
        mem_mb=config['get_bioproject_reads']['mem_mb'],
        runtime=config['get_bioproject_reads']['runtime'],
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        prefetch \
            --option-file {input:q} \
            --output-directory {output.dir_out:q} \
            2> {log:q}
        for file in {output.dir_out:q}/*/*.sra; do \
            fasterq-dump \
                --split-files \
                --outdir {output.dir_out:q} \
                --temp {output.dir_out:q}/tmp \
                --threads {threads} \
                --mem {resources.mem_mb}MB \
                "$file" \
                2>> {log:q}
        done
        pigz -p {threads} {output.dir_out:q}/*.fastq
        echo "Done" > {output.fn_done:q}
        """

rule get_srr_info:
    output:
        srr_info = fmt_srr_info,
        srr_list = fmt_srr_list
    log:
        "logs/get_srr_info/{exp_srr}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        srrs = lambda w: DICT_EXP[w.exp_srr]['srrs']
    shell:
        """
        > {output.srr_info:q}
        > {output.srr_list:q}
        > {log:q}
        for bs in {params.srrs}; do
            esearch -db sra -query "$bs" \
                | efetch -format runinfo \
                >> {output.srr_info:q} \
                2>> {log:q}
        done

        cat {output.srr_info:q} \
            | cut -d ',' -f 1 \
            | grep SRR \
            >> {output.srr_list:q} \
            2>> {log:q}        
        """


rule prefetch_srr_reads:
    input:
        fmt_srr_list,
    output:
        dir_out = temp(directory(dir_srr_prefetch)),
    log:
        "logs/prefetch_srr_reads/{exp_srr}.log"
    benchmark:
        "benchmarks/prefetch_srr_reads/{exp_srr}.benchmark.txt"
    threads:
        config['get_bioproject_reads']['threads'],
    resources:
        mem_mb=config['get_bioproject_reads']['mem_mb'],
        runtime=config['get_bioproject_reads']['runtime'],
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        cat {input:q} | parallel -j {threads} prefetch {{}} \
            --output-directory {output.dir_out:q} \
            > {log:q} 2>&1
        """

rule get_srr_reads:
    input:
        dir_srr_prefetch,
    output:
        dir_out = temp(directory(dir_srr_fastq)),
        fn_done = fmt_srr_fastq_done,
    log:
        "logs/get_srr_reads/{exp_srr}.log"
    benchmark:
        "benchmarks/get_srr_reads/{exp_srr}.benchmark.txt"
    threads:
        config['get_bioproject_reads']['threads'],
    resources:
        mem_mb=config['get_bioproject_reads']['mem_mb'],
        runtime=config['get_bioproject_reads']['runtime'],
    conda:
        "../envs/sra_tools.yaml"
    shell:
        """
        for file in {input:q}/*/*.sra; do \
            fasterq-dump \
                --split-files \
                --outdir {output.dir_out:q} \
                --temp {output.dir_out:q}/tmp \
                --threads {threads} \
                --mem {resources.mem_mb}MB \
                "$file" \
                2>> {log:q}
        done
        pigz -p {threads} {output.dir_out:q}/*.fastq
        echo "Done" > {output.fn_done:q}
        """


# TESTING: run download in parallel
# rule get_srr_reads:
#     output:
#         fn_out = temp(fmt_srr_fastq),
#         dir_out = temp(directory(dir_srr_fastq)),
#     log:
#         "logs/get_srr_reads/{exp_srr}/{srr}.log"
#     benchmark:
#         "benchmarks/get_srr_reads/{exp_srr}/{srr}.benchmark.txt"
#     conda:
#         "../envs/sra_tools.yaml"
#     shell:
#         """
#         prefetch \
#             {wildcards.srr} \
#             --output-directory {output.dir_out:q} \
#             2> {log:q}
        
#         file={output.dir_out:q}/{wildcards.srr}/{wildcards.srr}.sra
#         fasterq-dump \
#             --split-files \
#             --outdir {output.dir_out:q} \
#             --temp {output.dir_out:q}/tmp \
#             "$file" \
#             2>> {log:q}

#         pigz -p {output.dir_out:q}/{wildcards.srr}_*.fastq
#         """


# rule srr_download_done:
#     input:
#         get_fns_srr_fastq,
#     output:
#         fmt_srr_fastq_done
#     shell:
#         """
#         echo {input:q} > {output:q}

#         """

rule read_and_count_downloads_done:
    input:
        get_fns_downloads_done,
    output:
        fn_read_and_count_downloads_done,
    shell:
        """
        echo {input:q} > {output:q}
        """

