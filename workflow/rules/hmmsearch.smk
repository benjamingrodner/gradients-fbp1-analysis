
fmt_seq_parquet = (
    config['dirs_data']['metatranscriptomes'] 
    + '/merge_counts_tax_gene-{batch}.parquet'
)

rule get_fasta_from_parquet:
    input:
        fmt_seq_parquet,
    output:
        fmt_seq,
    params: 
        script="scripts/parquet2fasta.py",
    log:
        "logs/get_fasta_from_parquet/{batch}.log"
    benchmark:
        "benchmarks/get_fasta_from_parquet/{batch}.benchmark.txt"
    conda: 
        "../envs/ibis.yaml"
    shell:
        """
        python {params.script:q} \
            -i {input:q} \
            -o {output:q} \
            -hc contig_name_6tr \
            -h all \
            2> {log:q}
        """
