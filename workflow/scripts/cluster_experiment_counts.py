# framework written by Gemini 3 flash free tier 5/17/26
## Prompt
# Please write a python script taking click arguments of filenames for a counts parquet file, a config yaml, a cluster tsv table, a sample info text file, a string indicating the experiment, and an output parquet. Load the parquets and the csv using ibis with duckb backend. Load the yaml config as well. 

import re
import os
import click
import ibis
import yaml
import pandas as pd
from ibis import _, selectors
from collections import defaultdict

def aggregate_mapping(t, mapping, col):
    # mamke groups for contigs that we clustered, 
    # Leave other contigs unclustered
    if not all([c in t[col].to_pandas().values for c in mapping.keys()]):
        raise ValueError(f"One or more clusters do not align with the counts table contig names.")
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
    '--exp-info-path', '-e', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the experiment info YAML file.'
)
@click.option(
    '--experiment', '-x', 
    required=True, 
    type=str, 
    help='Experiment name/string identifier.'
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
    counts: str, clusters: str, colname_contig: str, exp_info_path: str, 
    experiment: str, output: str,
    mem: str, threads: int
):
    """Process counts and cluster data using Ibis and DuckDB."""

    with open(exp_info_path, 'r') as f:
        exp_info = yaml.safe_load(f)

    dict_exp_info = exp_info[experiment]
    re_subs = dict_exp_info.get('re_subs_fasta2counts')
    # prefix = dict_exp_info['prefix']

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
            dict_gene_contigs = defaultdict(list)
            for fn_cl in clusters:
                clusters_table = con.read_csv(
                    fn_cl, table_name="clusters", sep="\t", 
                    column_names=["cluster","contig"], header=False
                ).to_pandas()
                # bn_cl = os.path.basename(fn_cl)
                for clust, cont in clusters_table.values:
                    # remove prefix and 6tr frame
                    cs = []
                    for c in [clust, cont]:
                        c_ = re.sub(r'_\d+$','',c)
                        # c_ = re.sub(f'{prefix}_','', c_)
                        if re_subs is not None:
                            for subs in re_subs:
                                c_ = re.sub(subs[0],subs[1],c_)
                        cs.append(c_)
                    
                    dict_contig_clust[cs[1]] = cs[0]
                    dict_gene_contigs[fn_cl].append(cs[1])

    except Exception as e:
        raise click.ClickException(f"Error loading data tables: {e}")

    click.echo("Data successfully loaded.")
    # Get mapping

    # Aggregate counts into clusters
    if clusters is not None:
        # aggregate contig clusters
        out_table = aggregate_mapping(counts_table, dict_contig_clust, colname_contig)
        out_table = out_table.to_pandas()
        # aggregate all contigs for each gene
        for gene, contigs in dict_gene_contigs.items():
            df_gene = counts_table.filter(
                counts_table[colname_contig].isin(contigs)
            ).agg(
                selectors.across(selectors.numeric(), _.sum())
            ).to_pandas()
            df_gene['cluster'] = gene
            out_table = pd.concat([out_table, df_gene], axis=0, ignore_index=True)
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

        out_table.to_parquet(output)
        click.echo("Processing complete!")

    except Exception as e:
        raise click.ClickException(f"Error writing output Parquet: {e}")


if __name__ == "__main__":
    main()