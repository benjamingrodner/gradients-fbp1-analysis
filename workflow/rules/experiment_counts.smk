rule download_counts:
    output:
        fmt_download_counts,
    conda:
        "../envs/sra_tools.yaml"
    params:
        link = lambda w: DICT_ISO[w.iso]['link_counts']
    shell:
        """
        wget {params.link:q} -O {output:q}
        """

rule merge_download_counts:
    input:
        get_fns_download_counts,
    output:
        fn_download_counts_done,
    shell:
        """
        cat {input:q} > {output:q}
        """

rule get_bioproject_info:
    output:
        bp_info = fmt_bioproject_info,
        srr_info = fmt_bioproject_srr_info,
        srr_list = fmt_bioproject_srr_list,
    log:
        "logs/get_bioproject_info/{iso}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        bioproject = lambda w: DICT_ISO[w.iso]['bioproject']
        elink_target = lambda w: get_elink_target(w.iso)
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
        dir_out = directory(dir_bioproject_fastq),
        fn_done = fmt_bioproject_fasta_done,
    log:
        "logs/get_bioproject_reads/{gene}.log"
    benchmark:
        "benchmarks/get_bioproject_reads/{gene}.benchmark.txt"
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
            --output-directory {output:q} \
            2> {log:q}
        for file in {output:q}/*/*.sra; do \
            fasterq-dump \
                --split-files \
                --outdir {output.dir_out:q} \
                --temp {output:q}/tmp \
                --threads {threads} \
                --mem {resources.mem_mb}MB \
                "$file" \
                2>> {log:q}
        done
        cat "Done" > {output.fn_done:q}
        """

rule merge_bioproject_done:
    input:
        get_fns_bioproject_done,
    output:
        fn_bioproject_merge,
    shell:
        """
        cat {input:q} > {output:q}
        """

rule get_biosample_info:
    output:
        bp_info = fmt_biosample_info,
        srr_info = fmt_biosample_srr_info,
        srr_list = fmt_biosample_srr_list,
    log:
        "logs/get_biosample_info/{iso}.log"
    conda:
        "../envs/sra_tools.yaml"
    params:
        biosamples = lambda w: DICT_ISO[w.iso]['biosamples']
    shell:
        """
        for bs in {params.biosamples}; do
            esearch -db biosample -query "$bs" \
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
        done
        """

rule get_biosample_reads:
    input:
        fmt_biosample_srr_list,
    output:
        dir_out = directory(dir_biosample_fastq),
        fn_done = fn_biosample_fasta_done,
    log:
        "logs/get_biosample_reads/{gene}.log"
    benchmark:
        "benchmarks/get_biosample_reads/{gene}.benchmark.txt"
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
            --output-directory {output:q} \
            2> {log:q}
        for file in {output:q}/*/*.sra; do \
            fasterq-dump \
                --split-files \
                --outdir {output.dir_out:q} \
                --temp {output:q}/tmp \
                --threads {threads} \
                --mem {resources.mem_mb}MB \
                "$file" \
                2>> {log:q}
        done
        cat "Done" > {output.fn_done:q}
        """

rule merge_biosamples_done:
    input:
        get_fns_biosamples_done,
    output:
        fn_biosamples_merge,
    shell:
        """
        cat {input:q} > {output:q}
        """

rule align_reads:
    input:
        assembly = lambda w: get_fn_assembly(w.iso),
        reads = lambda w: get_fn(w.iso)
    output:
        


