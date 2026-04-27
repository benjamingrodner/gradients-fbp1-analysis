# Generated with Gemini 3 Flash, operating in the Free tier
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

import os
import glob

configfile: "config.yaml"

def get_target_seqs(wildcards):
    files = []
    for gene in config["target_genes"]:
        files.extend(glob.glob(f"{config['input_dir']}/*{gene}*.fasta"))
    return files

rule all:
    input:
        "results/final_placement.jplace"

rule merge_seqs_to_place:
    input:
        get_target_seqs
    output:
        merged = "work/merged_queries.fasta"
    log:
        "logs/merge_seqs.log"
    shell:
        "cat {input} > {output.merged} 2> {log}"

rule build_alignment_hmm:
    input:
        ref_msa = config["ref_msa"]
    output:
        hmm = "work/reference.hmm"
    log:
        "logs/hmmbuild.log"
    conda:
        "../envs/hmmer.yaml"
    threads: config["resources"]["threads"]
    resources:
        mem_mb = config["resources"]["mem_mb"],
        runtime = config["resources"]["runtime"]
    shell:
        "hmmbuild --cpu {threads} {output.hmm} {input.ref_msa} 2> {log}"

rule align_seqs_to_place:
    input:
        hmm = "work/reference.hmm",
        queries = "work/merged_queries.fasta",
        ref_msa = config["ref_msa"]
    output:
        aln = "work/queries_aligned.sto"
    log:
        "logs/hmmalign.log"
    conda:
        "../envs/hmmer.yaml"
    threads: config["resources"]["threads"]
    shell:
        "hmmalign --cpu {threads} --trim --mapali {input.ref_msa} "
        "{input.hmm} {input.queries} > {output.aln} 2> {log}"

rule trim_alignment_by_index:
    input:
        aln = "work/queries_aligned.sto"
    output:
        trimmed = "work/queries_trimmed.fasta"
    log:
        "logs/trimming.log"
    conda:
        "../envs/seqkit.yaml"
    params:
        start = config["trim_indices"]["start"],
        end = config["trim_indices"]["end"]
    shell:
        "seqkit export fa {input.aln} 2> {log} | "
        "seqkit subseq -r {params.start}:{params.end} > {output.trimmed} 2>> {log}"

rule mask_alignment_by_log:
    input:
        aln = "work/queries_trimmed.fasta",
        log_file = config["clipkit_log"]
    output:
        masked = "work/queries_masked.fasta"
    log:
        "logs/masking.log"
    conda:
        "../envs/python.yaml"
    shell:
        "python scripts/mask_alignment.py "
        "--input {input.aln} "
        "--log {input.log_file} "
        "--output {output.masked} 2> {log}"

rule append_alignment:
    input:
        ref = config["ref_msa"],
        query = "work/queries_masked.fasta"
    output:
        combined = "work/combined_alignment.fasta"
    log:
        "logs/append.log"
    shell:
        "cat {input.ref} {input.query} > {output.combined} 2> {log}"

rule place_seqs:
    input:
        msa = "work/combined_alignment.fasta",
        tree = config["ref_tree"]
    output:
        jplace = "results/final_placement.jplace"
    log:
        "logs/raxml_epa.log"
    conda:
        "../envs/raxml.yaml"
    threads: config["resources"]["threads"]
    resources:
        mem_mb = config["resources"]["mem_mb"],
        runtime = config["resources"]["runtime"]
    params:
        model = config["model"],
        prefix = "work/epa_run"
    shell:
        "raxml-ng --epa --msa {input.msa} --tree {input.tree} "
        "--model {params.model} --threads {threads} "
        "--prefix {params.prefix} 2> {log} && "
        "mv {params.prefix}.raxml.epa.jplace {output.jplace} 2>> {log}"