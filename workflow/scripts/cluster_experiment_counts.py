# framework written by Gemini 3 flash free tier 5/17/26
## Prompt
# Please write a python script taking click arguments of filenames for a counts parquet file, a config yaml, a cluster tsv table, a sample info text file, a string indicating the experiment, and an output parquet. Load the parquets and the csv using ibis with duckb backend. Load the yaml config as well. 

import os
import click
import ibis
import yaml
from ibis import _, selectors

def aggregate_mapping(t, mapping, col):
    # mamke groups for contigs that we clustered, 
    # Leave other contigs unclustered
    category_col = t[col].cases(*mapping.items(), else_ = t[col]).name("cluster")
    return t.group_by(category_col).agg(
        selectors.across(selectors.numeric(), _.sum())
    )
def load_config(config_path: str) -> dict:
    """Loads and returns the YAML configuration."""
    try:
        with open(config_path, "r") as f:
            return yaml.safe_load(f)
    except Exception as e:
        raise click.ClickException(f"Error reading config YAML: {e}")


@click.command()
@click.option(
    "--counts",
    "-c",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the input counts Parquet file.",
)
@click.option(
    "--clusters",
    "-l",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    required=False,
    default=None,
    help="Paths to the input cluster TSV files.",
)
@click.option(
    "--colname_contig",
    "-n",
    type=str,
    required=True,
    help="Name of the column with contig names to group into clusters.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Path where the output Parquet file will be saved.",
)
@click.option(
    "--mem",
    "-m",
    type=str,
    required=True,
    help="Memory to use in e.g. 50G",
)
@click.option(
    "--threads",
    "-t",
    type=int,
    required=True,
    help="Number of threads to use",
)


def main(
    counts: str, clusters: str, colname_contig: str, output: str,
    mem: str, threads: int
):
    """Process counts and cluster data using Ibis and DuckDB."""

    # 2. Initialize the Ibis DuckDB backend
    # Using an ephemeral in-memory DuckDB connection
    con = ibis.duckdb.connect(memory_limit=mem, threads=threads)

    click.echo("Loading data files into Ibis...")
    try:
        # 3. Load the Parquet file
        counts_table = con.read_parquet(counts, table_name="counts")

        # 4. Load the TSV files and add to mapping
        if clusters is not None:
            dict_contig_clust = {}
            for fn_cl in clusters:
                clusters_table = con.read_csv(
                    fn_cl, table_name="clusters", sep="\t", 
                    column_names=["cluster","contig"], header=False
                ).to_pandas()
                for clust, cont in clusters_table.values:
                    dict_contig_clust[cont] = clust

    except Exception as e:
        raise click.ClickException(f"Error loading data tables: {e}")

    click.echo("Data successfully loaded.")
    # Get mapping

    # Aggregate counts into clusters
    if clusters is not None:
        out_table = aggregate_mapping(counts_table, dict_contig_clust, colname_contig)
    else:
        out_table = counts_table
        click.echo(f"No clusters passed, copying over the input counts table")

    # Write table
    click.echo(f"Writing output to: {output}")
    try:
        # Ensure the output directory exists
        out_dir = os.path.dirname(output)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        con.to_parquet(out_table, output)
        click.echo("Processing complete!")

    except Exception as e:
        raise click.ClickException(f"Error writing output Parquet: {e}")


if __name__ == "__main__":
    main()