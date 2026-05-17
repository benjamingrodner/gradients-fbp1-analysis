
# adapted from chris' script aggregate_counts.py
# Function written with Gemini 3 flash free tier 5/15/26
### Prompt
# Please write a python function using ibis with duckdb backend to load a table from a file that may be either tab or comma separated and write the table to a parquet and set threads

"""
Script to copy a csv to parquet
"""
import logging
from pathlib import Path

import ibis
import click


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s"
)


def convert_csv_to_parquet(
    input_file_path: str | Path, 
    output_parquet_path: str | Path,
    threads: int | None = None,
) -> None:
    """Loads a CSV or TSV file using Ibis/DuckDB and writes it to a Parquet file.

    Automatically detects whether the delimiter is a comma or a tab.
    """
    # 1. Connect to the ephemeral DuckDB backend
    con = ibis.duckdb.connect()
    
    # 2. Apply performance configurations if provided
    if threads is not None:
        con.raw_sql(f"SET threads = {threads};")
    # 2. Reference the input file.
    # DuckDB's read_csv auto-detects delimiters (like ',' and '\t') by default.
    table = con.read_csv(input_file_path)

    # 3. Write the table expression directly to Parquet
    # Ibis executes this efficiently out-of-core via DuckDB
    table.to_parquet(output_parquet_path)

@click.command()
@click.option("--jobs", "-j", type=int, default=1, show_default=True, help="Number of CSV reader jobs")
@click.argument("input")
@click.argument("output", type=click.Path(exists=False))
def main(jobs, input, output):
    convert_csv_to_parquet(input, output, threads=jobs)


if __name__ == "__main__":
    main()
