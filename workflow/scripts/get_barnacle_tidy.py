#!/usr/bin/env python3

# Written by gemini 3.5 flash in free tier 6 Jul 2026
### Prompt
# Please write a python script that does the following:
# - Loads a parquet file using ibis with duckdb
# - Loads a config yaml 
# - Loads a text file as a list of strings
# - For values in the parquet column 'target_id' that match values in the text file, set the value in the parquet column ''KO" as "FBP1" 
# - Subset the parquet file to only those rows where the "phylum" column is one of ['Bacillariophyta','Haptophyta','Chlorophyta'] or the "class" column is one of ['Pelagophyceae','Dinophyceae'].
# - Group the subset rows by the phylum column for ['Bacillariophyta','Haptophyta','Chlorophyta'] rows and by the class column for ['Pelagophyceae','Dinophyceae']. The within those groups, group the rows by the "KO" column, summing the values in columns that are not ['taxid','domain','kingdom','phylum','class','order','family','genus','species','KO','score','target_id_right ','contig_name_6tr','target_id'].
# - Using those groupings make a dataframe with columns "taxname" for the first grouping and "KO" for the second grouping and then all the summed columns. 
# - Now convert this dataframe into "tidy" format with columns "taxname", "KO", "sample" (where a value in this column is a summed column name), "counts" (where a value in this column is the summed value for that column name).
# - Write this dataframe to a csv
# Use click arguments for inputs and add options for setting the number of threads and memory in ibis
# Make it so that the specific strings I have provided are taken from the config file
###

# Edited by Ben Grodner

import sys
import click
import ibis
import pandas as pd
import yaml
import re

@click.command()
@click.option(
    '--parquet-file', '-p', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the input Parquet file.'
)
@click.option(
    '--config-file', '-c', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the YAML configuration file.'
)
@click.option(
    '--target-gene-ids', '-t', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the text file containing target IDs (one per line).'
)
@click.option(
    '--metadata-file', '-m', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the text file containing metadata for counts columns.'
)
@click.option(
    '--output-file', '-o', 
    required=True, 
    type=click.Path(writable=True), 
    help='Path where the final tidy CSV will be saved.'
)
@click.option(
    '--threads', type=int, default=None,
    help='Number of threads for DuckDB execution. Defaults to system maximum.'
)
@click.option(
    '--memory', type=str, default=None,
    help='Memory limit for DuckDB (e.g., "8GB", "500MB"). Defaults to system maximum.'
)
def main(parquet_file, config_file, target_gene_ids, metadata_file, output_file, threads, memory):
    """
    Process taxonomic data using Ibis (DuckDB) and Pandas.
    Filters, updates, aggregates, and transforms data into a tidy format.
    """
    click.echo("🚀 Starting the pipeline...")

    # 1. Loading Setup inputs
    try:
        click.echo(" -> Loading configuration and input text files...")
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)
        
        with open(target_gene_ids, "r") as f:
            target_ids = [line.strip() for line in f if line.strip()]
            
    except Exception as e:
        click.echo(f"❌ Error loading input files: {e}", err=True)
        sys.exit(1)

    # 2. Setup Ibis Backend with DuckDB and resource limits
    click.echo(" -> Connecting to DuckDB backend via Ibis...")
    
    # Passing specific configurations directly to DuckDB backend setup
    con = ibis.duckdb.connect(
        threads=threads,
        memory_limit=memory
    )

    # Print confirmed limits if user provided them
    if threads:
        click.echo(f"   - Execution threads limited to: {threads}")
    if memory:
        click.echo(f"   - Execution memory limit set to: {memory}")

    try:
        t = con.read_parquet(parquet_file)
    except Exception as e:
        click.echo(f"❌ Error reading Parquet file: {e}", err=True)
        sys.exit(1)

    # 3. Data Transformation 

    # target gene naming
    tg_map = [(t[config['contig_colname']].isin(target_ids), config['target_gene_name'])]
    tg_mutate = {
        config['ko_colname']: ibis.cases( # Add ko name for target gene
            *tg_map, else_=t[config['ko_colname']]
        )
    }
    tn_map = [] # taxonomic level and taxa to include
    for tax, inlist in config['dict_tax_include'].items():
        tn_map.append((t[tax].isin(inlist), t[tax]))

    tn_mutate = {
        config['new_tax_colname']: ibis.cases(*tn_map, else_=None)
    }

    click.echo(" -> Updating 'KO' column and applying taxonomic filters...")
    t_subset = t.mutate(
        **tg_mutate
    ).mutate( # new column for taxon grouping
        **tn_mutate
    ).drop_null([config['new_tax_colname']]) # Filter for target taxa only

    # target_phyla = config['target_phyla']
    # target_classes = config['target_classes']

    # t_subset = t.filter(t.phylum.isin(target_phyla) | t.class_.isin(target_classes))

    # 4. Grouping Logic & Dynamic Summation
    click.echo(" -> Structuring grouping levels and computing sums...")

    # t_subset = t.mutate(
    #     taxname=ibis.cases(*tn_map, else_=None)
    # ).drop_null(subset=['taxname'])

    sum_cols = [col for col in t_subset.columns if col not in config['sum_exclude_cols']]
    
    if not sum_cols:
        click.echo("❌ Error: No valid numeric columns found to sum after exclusions.", err=True)
        sys.exit(1)

    # Build aggregations dictionary dynamically
    aggs = {col: t_subset[col].sum() for col in sum_cols}
    grouped_expr = t_subset.group_by(
        [config['new_tax_colname'], config['ko_colname']]
    ).aggregate(**aggs)

    # 5. Execution and Tidy Transformation (Pandas)
    click.echo(" -> Fetching results into memory and converting to tidy format...")
    try:
        wide_df = grouped_expr.execute()
    except Exception as e:
        click.echo(f"❌ Error during DuckDB query execution: {e}", err=True)
        sys.exit(1)

    tidy_df = pd.melt(
        wide_df,
        id_vars=[config['new_tax_colname'], config['ko_colname']],
        value_vars=sum_cols,
        var_name='kallisto_fn',
        value_name=config['counts_colname']
    )

    # Get sample name and rep from metadata
    meta = pd.read_csv(metadata_file)
    meta_sub = meta[config['metadata_columns']]

    tidy_df = tidy_df.merge(meta_sub, 
        how='left', 
        left_on=config['values_merge_colname'],
        right_on=config['meta_merge_colname'],
    )
    
    samples, filts, reps = [],[],[]
    for sr in tidy_df[config['sample_filt_rep_colname']]:
        match = re.search(config['sample_filt_rep_regex'], sr)

        samples.append(match.group('sample'))
        filts.append(match.group('filt'))
        reps.append(match.group('rep'))
    
    tidy_df[config['sample_colname']] = samples
    tidy_df[config['filt_colname']] = filts
    tidy_df[config['rep_colname']] = reps

    # Sum between filters
    tidy_df[config['trl_colname']] = (
        tidy_df[config['counts_colname']] * tidy_df[config['norm_colname']]
    )
    tidy_df_grp = tidy_df.groupby([
        config['new_tax_colname'], config['ko_colname'], 
        config['sample_colname'], config['rep_colname'],
    ])[config['trl_colname']].sum().reset_index()

    # 6. Exporting Results
    try:
        click.echo(f" -> Saving final output to: {output_file}")
        tidy_df_grp.to_csv(output_file, index=False)
        click.echo("🎉 Pipeline finished successfully!")
    except Exception as e:
        click.echo(f"❌ Error saving output CSV: {e}", err=True)
        sys.exit(1)

if __name__ == '__main__':
    main()