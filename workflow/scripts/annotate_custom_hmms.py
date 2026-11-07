#!/usr/bin/env python3

## Prompt to gemini:
# You are an expert Python data engineer with strong experience in ibis, 
# DuckDB, parquet, and parsing HMMER output. Write production-quality Python code to solve the following task.

# Task
# I have:

# a large parquet table with columns including target_id, KO, and optionally E-value and/or score;
# several HMMER hmmsearch --tblout result tables, each associated with a specific gene name.
# For each row in the parquet table:

# if target_id appears in any HMMER result table, consider the corresponding gene name as a candidate replacement for KO;
# replace KO with that gene name only if the HMMER hit is better than the existing row according to the relevant metric:
# lower E-value is better;
# higher score is better.
# The code must use ibis with a DuckDB backend to read and transform the parquet table.

# Requirements
# Read the large parquet table lazily with ibis + DuckDB.
# Parse multiple HMMER --tblout files, each associated with a gene name.
# Normalize the HMMER tables into a consistent schema.
# Combine all HMMER hits into one table.
# If multiple HMMER rows match the same target_id, reduce them to the best hit per target using the correct comparison rule.
# Join the HMMER hits to the parquet table on target_id.
# Update KO only when the HMMER hit is better than the existing row.
# Preserve all original rows and columns from the parquet table.
# Write the updated result to a parquet file.
# Include clear, minimal comments and make the code runnable as a script.
# Important details
# Assume the HMMER --tblout files are tabular but may need parsing because they are whitespace-delimited rather than clean TSV.
# Handle missing values robustly.
# The code should be efficient for large tables and avoid loading the full parquet table into pandas unless absolutely necessary.
# Prefer clean ibis expressions over manual row-by-row Python loops.
# If E-value and score both exist, implement a sensible rule and make it explicit in the code.
# If the update rule is ambiguous, choose a clear precedence and document it briefly in a comment.
# Include command-line arguments for:
# input parquet path,
# output parquet path,
# one or more HMMER --tblout files,
# gene name for each file, or a way to infer it from filename.
# Suggested output structure
# Return:

# A short explanation of the approach.
# A complete Python script.
# Any assumptions made.
# Code quality requirements
# Use functions for parsing and transformation.
# Include type hints where helpful.
# Validate inputs and give informative errors.
# Keep the implementation practical and concise.
# Do not use placeholder pseudocode; produce real code.
# Clarifying assumptions to make if needed
# If the input parquet has both E-value and score, the code should define whether:

# either metric can trigger replacement,
# both must agree, or
# one metric takes precedence.
# Choose a reasonable approach and implement it consistently.
# Use ibis for the main join/update logic, and only use pandas for parsing the small HMMER tables if necessary.

### Code edited by perplexity

### Edits by Ben grodner


import argparse
import re
from pathlib import Path
from typing import List, Optional

import ibis
import pandas as pd


def parse_hmmer_tblout(filepath: Path, gene_name: str) -> pd.DataFrame:
    if not filepath.exists():
        raise FileNotFoundError(f"HMMER file not found: {filepath}")

    df = pd.read_csv(
        filepath,
        sep=r"\s+",
        comment="#",
        header=None,
        engine="python"
    )
    columns = [
        "target_id_6tr_hmm", "target_accession", "query_name", "query_accession",
        "hmmer_evalue", "hmmer_score", "full_bias", "best_domain_evalue",
        "best_domain_score", "best_domain_bias", "exp", "reg", "clu", "ov",
        "env", "dom", "rep", "inc", "description"
    ]
    if not df.shape[1] == len(columns):
        raise ValueError(f"Hmmer table in {filepath} does not have the correct number of columns. Needed: {len(columns)}, Present: {df.shape[1]}")
    
    df.columns = columns
    df['hmmer_gene'] = gene_name

    return df[["target_id_6tr_hmm", "hmmer_gene", "hmmer_evalue", "hmmer_score"]]


def parse_hmmer_inputs(
    tblout_files: List[str],
    join_key: str,
    gene_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    dfs = []

    for idx, file_str in enumerate(tblout_files):
        file_path = Path(file_str)
        gene = gene_names[idx] if gene_names and idx < len(gene_names) else file_path.stem.split(".")[0]

        df = parse_hmmer_tblout(file_path, gene)
        if not df.empty:
            df[join_key] = df["target_id_6tr_hmm"].str.replace(r'_\d+$', '', regex=True)
            dfs.append(df)

    if not dfs:
        return pd.DataFrame(
            columns=[join_key, "target_id_6tr_hmm", "hmmer_gene", "hmmer_evalue", "hmmer_score"]
        )

    return pd.concat(dfs, ignore_index=True)


def get_metric_columns(columns: List[str]):
    evalue_col = "eval_ko" if "eval_ko" in columns else None
    score_col = "score" if "score" in columns else None
    if (evalue_col is None) and (score_col is None):
        raise ValueError("Either 'eval_ko' or 'score' must be a column in the data table, neither is present.")
    return evalue_col, score_col


def get_best_hmmer_df(hmmer_raw: pd.DataFrame, join_key: str) -> pd.DataFrame:
    return (
        hmmer_raw
        .sort_values(
            by=[join_key, 'hmmer_evalue', 'hmmer_score'],
            ascending=[True, True, False]
        )
        .groupby(join_key, as_index=False)
        .head(1)
    )


def build_replacement_expr(joined, columns: List[str], col_6tr: str):
    evalue_col, score_col = get_metric_columns(columns)
    hmmer_present = joined.hmmer_gene.notnull()
    
    new_score = None
    new_eval = None
    if evalue_col and score_col:
        orig_e = joined[evalue_col]
        orig_s = joined[score_col]
        better = (
            orig_e.isnull()
            | (joined.hmmer_score > orig_s)
            | ((joined.hmmer_score == orig_s) & (joined.hmmer_evalue > orig_e))
        )
        new_score = ibis.ifelse(
            hmmer_present & better, joined.hmmer_score, joined[score_col]
        )
        new_eval = ibis.ifelse(
            hmmer_present & better, joined.hmmer_evalue, joined[evalue_col]
        )
    elif evalue_col:
        orig_e = joined[evalue_col]
        better = orig_e.isnull() | (joined.hmmer_evalue < orig_e)
        new_eval = ibis.ifelse(
            hmmer_present & better, joined.hmmer_evalue, joined[evalue_col]
        )
    elif score_col:
        orig_s = joined[score_col]
        better = orig_s.isnull() | (joined.hmmer_score > orig_s)
        new_score = ibis.ifelse(
            hmmer_present & better, joined.hmmer_score, joined[score_col]
        )

    new_ko = ibis.ifelse(hmmer_present & better, joined.hmmer_gene, joined.KO)
    new_6tr = ibis.ifelse(hmmer_present & better, joined["target_id_6tr_hmm"], joined[col_6tr])
    return new_ko, new_score, new_eval, new_6tr


def project_updated_table(joined, columns: List[str], col_6tr: str, new_ko, new_score, new_eval, new_6tr):
    projected = []
    for col in columns:
        if col == "KO":
            projected.append(new_ko.name("KO"))
        elif col == 'score':
            if new_score is None:
                raise ValueError("'score' column present in main table, but replacement expression not created")
            projected.append(new_score.name("score_ko"))
        elif col == 'eval_ko':
            if new_eval is None:
                raise ValueError("'eval_ko' column present in main table, but replacement expression not created")
            projected.append(new_eval.name("eval_ko"))
        elif col == col_6tr:
            projected.append(new_6tr.name(col_6tr))
        else:
            projected.append(joined[col])
    return joined.select(*projected)


def build_updated_table(
    con: ibis.BaseBackend, 
    parquet_path: str, 
    hmmer_df: pd.DataFrame, 
    join_key: str, 
    col_6tr: str
):
    main_tbl = con.read_parquet(parquet_path)

    if hmmer_df.empty:
        print("No hmmer hits, no rows edited.")
        return main_tbl

    if "KO" not in main_tbl.columns:
        raise ValueError("Input parquet must contain a KO column.")
    if join_key not in main_tbl.columns:
        raise ValueError(f"Input parquet must contain a {join_key} column.")

    best_hmmer = get_best_hmmer_df(hmmer_df, join_key)
    best_hmmer = con.create_table("hmmer_raw", best_hmmer, temp=True)

    joined = main_tbl.left_join(best_hmmer, join_key)
    new_ko, new_score, new_eval, new_6tr = build_replacement_expr(
        joined, list(main_tbl.columns), col_6tr
    )
    return project_updated_table(
        joined, list(main_tbl.columns), col_6tr, new_ko, new_score, new_eval, new_6tr
    )


def main():
    parser = argparse.ArgumentParser(description="Annotate custom HMMs into parquet tables.")
    parser.add_argument("--ann-parquet", required=True, help="Path to input parquet annotation file")
    parser.add_argument("--output-parquet", required=True, help="Path to output custom annotation parquet file")
    parser.add_argument("--tbl-files", nargs="+", required=True, help="List of HMMER tblout files")
    parser.add_argument("--regex-pattern", required=True, help="Regex pattern to extract gene names from filenames")
    parser.add_argument("--join-key", required=True, help="Join key column name")
    parser.add_argument("--col-6tr", required=True, help="Target 6tr column name")
    parser.add_argument("--threads", type=int, default=4, help="Threads for DuckDB")
    parser.add_argument("--mem_mb", default=1000, help="Memory limit for DuckDB")

    args = parser.parse_args()

    # Extract gene names based on input regex pattern
    regex = re.compile(args.regex_pattern)
    gene_names = []
    for fn in args.tbl_files:
        match = regex.search(fn)
        if match:
            gene_names.append(match.group('gene'))
        else:
            raise ValueError(f"Regex did not match filename.\nregex: {regex}\nfilename: {fn}")

    hmmer_df = parse_hmmer_inputs(args.tbl_files, args.join_key, gene_names)

    mem = str(args.mem_mb) + 'MB'
    con = ibis.duckdb.connect(memory_limit=mem, threads=args.threads)
    result_tbl = build_updated_table(con, args.ann_parquet, hmmer_df, args.join_key, args.col_6tr)

    result_tbl.to_parquet(args.output_parquet)


if __name__ == "__main__":
    main()