"""
Script to aggregate kallisto or salmon counts.

Can read kallisto or salmon TSV files.

Output is in wide data format, with the columns "contig_name" and
"contig_length" plus a count column for each sample, named as `count_<sample>`.
This format assumes that kallisto and salmon output is in stable order. This can
be tested with duckdb, e.g. this should show one value for the checksum.

duckdb -csv \
    -c "SELECT filename, bit_xor(md5_number(target_id)) AS checksum FROM read_csv_auto('kallisto_sample_counts/*/*.abundance.tsv.gz', filename=true) GROUP BY filename;" \
    | awk -F, '{print $2}' | sort -u
"""
import logging
from multiprocessing import Pool
from pathlib import Path

import click
import duckdb
import pandas as pd
from tqdm import tqdm


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(message)s"
)


def get_file_type(file: str) -> str:
    if file.endswith("abundance.tsv.gz") or file.endswith("abundance.tsv"):
        return "kallisto"
    elif file.endswith("quant.sf.gz") or file.endswith("quant.sf"):
        return "salmon"
    else:
        raise ValueError(f"Unknown file type for {file}")


def read_kallisto_tsv(file: str, counts_only=False) -> pd.DataFrame:
    if not Path(file).exists():
        raise OSError(f"File {file} not found")
    df = duckdb.read_csv(file, sep="\t").to_df()
    if counts_only:
        df = df[["est_counts"]]
        df.rename(columns={"est_counts": "count"}, inplace=True)
    else:
        df = df[["target_id", "length", "est_counts",]]
        df.rename(columns={
            "target_id": "contig_name",
            "length": "contig_length",
            "est_counts": "count"
        }, inplace=True)
    return df


def read_salmon_tsv(file: str, counts_only=False) -> pd.DataFrame:
    if not Path(file).exists():
        raise OSError(f"File {file} not found")
    df = duckdb.read_csv(file, sep="\t").to_df()
    if counts_only:
        df = df[["NumReads"]]
        df.rename(columns={"NumReads": "count"}, inplace=True)
    else:
        df = df[["Name", "Length", "NumReads"]]
        df.rename(columns={
            "Name": "contig_name",
            "Length": "contig_length",
            "NumReads": "count"
        }, inplace=True)
    return df


def read_counts_file(work: tuple[int, str]) -> pd.DataFrame:
    i = work[0]
    file = work[1]
    # Sampe should be the name of the directory containing the counts file
    # e.g. kallisto_sample_counts/sample1/abundance.tsv.gz
    # or salmon_sample_counts/sample1/quant.sf.gz
    sample = Path(file).parts[-2]
    if i == 0:
        counts_only = False
    else:
        counts_only = True

    file_type = get_file_type(file)

    if file_type == "kallisto":
        df = read_kallisto_tsv(file, counts_only=counts_only)
    elif file_type == "salmon":
        df = read_salmon_tsv(file, counts_only=counts_only)
    # Should raise if file type not recognized
    df.rename(columns={"count": f"count_{sample}"}, inplace=True)
    return df


def agg(input: list[str], output, jobs: int=1):
    jobs = max(jobs, 1)
    # Assume all files have the same column of contig_name and contig_length,
    # only keep the contig_name and contig_length from the first file.
    logging.info(f"Aggregating {len(input)} files")
    logging.info(f"Using {jobs} jobs")
    for f in input:
        if not Path(f).exists():
            logging.error(f"{f} not found")
            raise OSError(f"{f} not found")
    with Pool(processes=jobs) as pool:
        dfs = []
        for df in tqdm(pool.imap_unordered(read_counts_file, list(enumerate(sorted(input)))), total=len(input)):
            dfs.append(df)
    if dfs:
        logging.info(f"Concatenating CSV file columns and writing to {output}")
        big_df = pd.concat(dfs, axis=1)
        logging.info(f"Final DataFrame shape = {big_df.shape}")
        count_cols = sorted([c for c in big_df.columns if c.startswith("count_")])
        big_df[["contig_name", "contig_length"] + count_cols].to_parquet(output, index=False)
        logging.info("Done")


@click.command()
@click.option("--jobs", "-j", type=int, default=1, show_default=True, help="Number of CSV reader jobs")
@click.argument("input", nargs=-1)
@click.argument("output", type=click.Path(exists=False))
def main(jobs, input, output):
    agg(input, output, jobs=jobs)


if __name__ == "__main__":
    main()
