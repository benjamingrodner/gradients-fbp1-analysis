# Generated with Gemini 3 Flash, operating in the Free tier 27 Apr 2026
# using the following prompt:

# "Please generate a complete Snakemake pipeline for phylogenetic placement based on the following specifications.
# 1. Directory Structure:

# Snakefile in the root.
# Conda environments in ../envs/ (python.yaml, raxml.yaml, hmmer.yaml, seqkit.yaml).
# External Python scripts in scripts/.
# Logs in logs/ and intermediate data in work/.
# Output files written to results/
# 2. Rules to Include:

# merge_seqs_to_place: Uses a Python input function to glob a directory, filter by a target_genes list in the config, and cat them into one FASTA. Redirect stderr to log.
# build_alignment_hmm: Runs hmmbuild using hmmer.yaml. Assign threads,  RAM, and runtime with reference to the config.
# align_seqs_to_place: Runs hmmalign with --trim and --mapali (referencing the original MSA). Use hmmer.yaml, assign threads, and RAM from the config.
# trim_alignment_by_index: Uses seqkit export fa piped into seqkit subseq to hard-trim the alignment based on 1-based start/end indices from the config. Use seqkit.yaml.
# mask_alignment_by_log: Calls a standalone Python script scripts/mask_alignment.py. Use python.yaml.
# append_alignment: Simple cat of the reference alignment and the masked query alignment.
# place_seqs: Runs raxml-ng --epa. Use raxml.yaml, and assign threads,  RAM, and  runtime from the config. Move the resulting .jplace to the final results folder.
# 3. Script Requirements:

# Provide the full code for scripts/mask_alignment.py.
# The script must use argparse, contain a main() function, and use Bio.AlignIO to mask columns with '-' based on indices found in a ClipKit log file.
# Include logic to handle potential 1-based to 0-based index conversion.
# 4. Global Requirements:

# All shell commands must redirect errors using 2> {log}.
# The all rule should target results/final_placement.jplace.
# Include a sample config.yaml structure at the top of the response."


### Edited while troubleshooting

rule merge_env_seqs_to_place:
    input:
        get_target_env_seqs
    output:
        merged = fn_target_env_seqs_merged
    log:
        "logs/merge_env_seqs_to_place.log"
    shell:
        "cat {input} > {output.merged} 2> {log}"

rule build_alignment_hmm:
    input:
        ref_msa = fn_alignment
    output:
        hmm = fn_alignment_hmm
    log:
        "logs/build_alignment_hmm.log"
    benchmark:
        "benchmarks/build_alignment_hmm.benchmark.txt"
    conda:
        "../envs/hmmer.yaml"
    threads: config["build_alignment_hmm"]["threads"]
    resources:
        mem_mb = config["build_alignment_hmm"]["mem_mb"],
        runtime = config["build_alignment_hmm"]["runtime"]
    shell:
        "hmmbuild --cpu {threads} {output.hmm} {input.ref_msa} 2> {log}"

rule align_env_seqs:
    input:
        hmm = fn_alignment_hmm,
        queries = fn_target_env_seqs_merged,
        ref_msa = fn_alignment
    output:
        aln = fn_env_aligned_sto
    log:
        "logs/align_env_seqs.log"
    conda:
        "../envs/hmmer.yaml"
    shell:
        "hmmalign --trim --mapali {input.ref_msa} "
        "-o {output.aln} {input.hmm} {input.queries} 2> {log}"

rule trim_env_alignment_by_index:
    input:
        aln = fn_env_aligned_sto,
        startend = fn_trim_crystal_startend,
    output:
        trimmed = fn_env_aligned_trim,
        fa = fn_env_aligned_fasta,
    log:
        "logs/trim_env_alignment_by_index.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqmagick convert --output-format fasta {input.aln} {output.fa} 2> {log}
        cat {input.startend:q} \
            | (IFS=',' read -r start end; 
                seqkit subseq -r "$start":"$end" {output.fa:q} \
                > {output.trimmed:q} \
                2>> {log:q}
            )
        """

rule mask_env_alignment_by_log:
    input:
        aln = fn_env_aligned_trim,
        log_file = fn_trim_clip + '.log'
    output:
        masked = fn_env_aligned_trim_mask
    log:
        "logs/mask_env_alignment_by_log.log"
    conda:
        "../envs/python.yaml"
    params:
        script = config['dir_scripts'] + "/mask_alignment_from_clipkit_log.py"
    shell:
        "python {params.script} "
        "--input {input.aln:q} "
        "--log {input.log_file:q} "
        "--output {output.masked} 2> {log}"

rule filter_env_seqs_by_len:
    input:
        fn_env_aligned_trim_mask,
    output:
        fn_env_aligned_trim_mask_filt,
    log:
        "logs/filter_env_seqs_by_len.log",
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

rule deduplicate_env_alignment:
    input:
        fasta = fn_env_aligned_trim_mask_filt
    output:
        fasta = fn_env_aligned_trim_mask_filt_dedup,
        mapping = fn_env_aligned_trim_mask_filt_dedup_map
    log:
        err = "logs/deduplicate_env_alignment.log"
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
            --map {output.mapping:q} 2> {log.err:q}
        """

rule place_env_on_tree:
    input:
        msa = fn_env_aligned_trim_mask_filt_dedup,
        tree = fn_full_tree_done
    output:
        fn_place_env_tree_done
    log:
        "logs/place_env_on_tree.log"
    benchmark:
        "benchmarks/place_env_on_tree.benchmark.txt"
    threads: config["build_tree"]["threads"]
    resources:
        mem_mb = config["build_tree"]["mem_mb"],
        runtime = config["build_tree"]["runtime"]
    params:
        model = config["build_tree"]['model'],
    shell:
        """
        DIR_OUT=$( dirname {output:q} )
        CWD=$( pwd )
        DIR_OUT="$CWD"/"$DIR_OUT"

        BN=$( basename {output:q} )
        BN="${{BN%.*}}"

        TREE={input.tree}
        DIR_TREE=$(dirname "$TREE")
        BN_TREE=$(basename "$TREE")
        BN_TREE="${{TREE%.*}}"
        FN_TREE="$DIR_TREE"/RAxML_bestTree."$BN_TREE"
        
        raxmlHPC-PTHREADS-AVX \
            -f v \
            -w "$DIR_OUT" \
            -s {input.msa} \
            -t "$FN_TREE" \
            -m {params.model} \
            -T {threads} \
            -n "$BN" \
            2> {log}

        cat "Done" > {output:q}
        """