# written initially by gemini flash 3 free tier 5/18/26
## prompt
# Please write a python script with:

# click arguments: 
# - counts table filename, a parquet file
# - metadata table filename, a csv file
# - config filename, a yaml file
# - experiment, a string
# - min_sum_counts, an integer
# - output directory for deseq stats
# - threads, an integer

# steps:
# - load the yaml, csv, parquet files.
# - make a new table with a subset the columns of the counts table (samples) such that the kept columns names are the index names from the metadata file and the new counts table index is defined by a column from the old table specified by config[experiment]['colname_contigs']
# - subset the rows of the new table such that the mean value of each row is greater than or equal to the min_sum_counts parameter. Print the number of rows kept and the number of rows removed
# - build a DeseqDataset with the the transpose of the new counts table data type int, the metadata table, design as "~condition", and the inference set with the threads argument
# - run deseq2
# - write out the DeseqStats for each comparison between config[experiment]['control_condition'] column and each config[experiment]['test_conditions'] column as a pkl file

# Note: use ibis with threads set by argument threads to subset the parquet file before loading to pandas
###

# edited for troubleshooting

import os
import re
import sys
import ibis
import yaml
import click
import pickle
import pandas as pd
from pydeseq2.ds import DeseqStats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference


@click.command()
@click.argument("counts_file", type=click.Path(exists=True))
@click.argument("metadata_file", type=click.Path(exists=True))
@click.argument("config_file", type=click.Path(exists=True))
@click.argument("experiment", type=str)
@click.argument("colname_contigs", type=str)
@click.argument("min_sum_counts", type=float)
@click.argument("output_dir", type=click.Path())
@click.argument("threads", type=int)
def main(
    counts_file,
    metadata_file,
    config_file,
    experiment,
    colname_contigs,
    min_sum_counts,
    output_dir,
    threads,
):
    """Run PyDESeq2 analysis """

    # -------------------------------------------------------------------------
    # 1. Load the configuration and metadata
    # -------------------------------------------------------------------------
    click.echo("Loading configuration and metadata...")
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    if experiment not in config:
        click.echo(
            f"Error: Experiment '{experiment}' not found in config.", err=True
        )
        sys.exit(1)

    exp_config = config[experiment]
    dict_ctrl_tests = exp_config['deseq_ctrl_tests']

    metadata_df = pd.read_csv(metadata_file, index_col=0)
    samples = metadata_df.index.tolist()

    # -------------------------------------------------------------------------
    # 2. Use Ibis to read and subset the Parquet file
    # -------------------------------------------------------------------------
    click.echo(
        f"Initializing Ibis (DuckDB backend) with {threads} threads..."
    )
    # Connect using DuckDB backend and set thread configuration
    con = ibis.duckdb.connect()
    con.raw_sql(f"SET threads TO {threads};")

    # Lazily reference the parquet file
    counts_table = con.read_parquet(counts_file)

    # Validate schema columns before executing
    table_columns = counts_table.columns
    if colname_contigs not in table_columns:
        click.echo(
            f"Error: Index column '{colname_contigs}' not found in counts table.",
            err=True,
        )
        sys.exit(1)

    missing_samples = [s for s in samples if s not in table_columns]
    if missing_samples:
        click.echo(
            f"Error: {len(missing_samples)} samples from metadata are missing "
            f"in counts columns. Examples: {missing_samples[:5]}",
            err=True,
        )
        sys.exit(1)

    # Get total row count before filtering
    initial_row_count = counts_table.count().execute()

    # Subset columns lazily
    columns_to_keep = [colname_contigs] + samples
    subset_table = counts_table.select(columns_to_keep)

    # Lazily calculate the row mean across sample columns
    # Ibis enables math operations across an array of columns natively
    row_sum = sum(subset_table[sample] for sample in samples) # / len(samples)

    # Apply row filter and execute the expression directly to a Pandas DataFrame
    click.echo("Filtering and loading data into memory...")
    filtered_table = subset_table.filter(row_sum >= min_sum_counts)
    new_counts = filtered_table.execute()

    # Set index now that data is in Pandas
    new_counts = new_counts.set_index(colname_contigs)

    # Report filtering metrics
    rows_kept = len(new_counts)
    rows_removed = initial_row_count - rows_kept
    click.echo(f"Rows kept (sum >= {min_sum_counts}): {rows_kept}")
    click.echo(f"Rows removed: {rows_removed}")

    if rows_kept == 0:
        click.echo(
            "Error: No rows left after filtering. Adjust min_sum_counts.",
            err=True,
        )
        sys.exit(1)


    # -------------------------------------------------------------------------
    # 4. Run DESeq2 and Write Outputs
    # -------------------------------------------------------------------------

    os.makedirs(output_dir, exist_ok=True)

    for control_condition, tests in dict_ctrl_tests.items():
        click.echo(f"Building DeseqDataset for control {control_condition}...")
        # Transpose matrix for PyDESeq2 compatibility (samples as rows, genes as columns)
        conds = [control_condition] + tests
        cols = metadata_df[metadata_df['condition'].isin(conds)].index
        counts_transposed = new_counts[cols].T.astype(int)
        metadata_sort = metadata_df.loc[counts_transposed.index]

        inference = DefaultInference(n_cpus=threads)
        design = f"C(condition, contr.treatment(base='{control_condition}'))"
        dds = DeseqDataSet(
            counts=counts_transposed,
            metadata=metadata_sort,
            design=f'~ {design}',
            inference=inference,
        )
        click.echo("Running DESeq2...")
        dds.deseq2()
        for test in tests:
            click.echo(f"Extracting stats for contrast: {test} vs {control_condition}")
            stat_res = DeseqStats(
                dds,
                contrast=["condition", test, control_condition],
                n_cpus=threads,
            )
            stat_res.summary()
            stat_res.lfc_shrink(coeff=f'{design}[T.{test}]'); # Account for high variance, low count, and low sample number
            out_filename = f"{experiment}_{test}_vs_{control_condition}_stats.pkl"
            out_path = os.path.join(output_dir, out_filename)
            out_path = re.sub(r'\s', '_', out_path)

            with open(out_path, "wb") as f:
                pickle.dump(stat_res, f)

            click.echo(f"Saved stats to {out_path}")

    click.echo("Analysis complete!")


if __name__ == "__main__":
    main()