# Framework written by gemini flash 3 free tier 5/18/26

## prompt
# Please write a python script with:

# click arguments: 
# - counts table filename, a parquet file
# - metadata table filename, a csv file
# - pydeseq2 stats result directory
# - config filename, a yaml file
# - experiment info filename, a yaml file
# - experiment, a string
# - output filename
# - threads, an integer

# steps:
# - load the yaml, csv
# - lazy load the parquet file using ibis with the specified number of threads
# - load DeseqStats objects from the .pkl files in the stats result directory
# - placeholder code for the operation of the script
##

# Plot code added and edited for troubleshooting

#!/usr/bin/env python3
import os
import glob
import ibis
import yaml
import click
import pickle
import numpy as np
import pandas as pd
import matplotlib as plt
from collections import defaultdict

@click.command()
@click.option(
    '--counts-path', '-c', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the counts table Parquet file.'
)
@click.option(
    "--clusters-format", "-l",
    type=str,
    required=True,
    help="Format string to the input cluster TSV files.",
)
@click.option(
    '--metadata-path', '-m', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the metadata table CSV file.'
)
@click.option(
    '--stats-dir', '-s', 
    required=True, 
    type=click.Path(exists=True, file_okay=False), 
    help='Directory containing pydeseq2 DeseqStats .pkl files.'
)
@click.option(
    '--config-path', '-g', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the config YAML file.'
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
    '--output-path', '-o', 
    required=True, 
    type=click.Path(writable=True), 
    help='Path for the output file.'
)
@click.option(
    '--threads', '-t', 
    required=True, 
    type=int, 
    default=4, 
    show_default=True, 
    help='Number of threads for Ibis/DuckDB processing.'
)
def main(counts_path, clusters_format, metadata_path, stats_dir, config_path, exp_info_path, experiment, output_path, threads):
    """
    Bioinformatics pipeline script to process counts, metadata, and PyDeSeq2 stats results.
    """
    click.echo(f"Starting pipeline for experiment: {experiment}")

    # Configure DuckDB backend via Ibis to use the specified threads
    con = ibis.duckdb.connect(threads=threads)
    ibis.set_backend(con)
    
    # -------------------------------------------------------------------------
    # Step 1: Load stuff
    # -------------------------------------------------------------------------
    click.echo("Loading files...")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    with open(exp_info_path, 'r') as f:
        exp_info = yaml.safe_load(f)
        
    metadata_df = pd.read_csv(metadata_path)
    
    # Lazy load the parquet file
    counts_table = ibis.read_parquet(counts_path)
    
    # -------------------------------------------------------------------------
    # Step 4: Plot
    # -------------------------------------------------------------------------
    
    exp_config = exp_info[experiment]
    colname_contigs = exp_config["colname_contigs"]
    dict_ctrl_tests = exp_config['deseq_ctrl_tests']

    plotg = config['experiment_fold_change_plots']
    plot_e = exp_config['plot']

    for gene in config['target_genes_for_exp_placement']:
        fig, ax = plt.subplots()
        # Load clusters for target gene
        dict_clust_contig = {}
        fn_cl = clusters_format.format(exp=experiment, gene=gene)
        clusters_table = ibis.read_csv(
            fn_cl, table_name="clusters", sep="\t", 
            column_names=["cluster","contig"]
        ).to_pandas()
        for clust, cont in clusters_table.values:
            dict_clust_contig[clust].appen(cont)

        # Get order of clusts based on sum of relative abundances
        colns_sam = metadata_df.index.to_list()
        df_counts = counts_table.select(colns_sam).filter(
            counts_table[colname_contigs].isin(dict_clust_contig)
        ).to_pandas()
        df_counts['relabund_sum'] = df_counts.div(
            df_counts.sum(axis=0), axis=1
            ).sum(axis=1)
        df_counts = df_counts.sort_values(by='relabund_sum', ascending=False)
        ncontigs = df_counts.shape[0]

        # Get each stat test
        nctrls = len(dict_ctrl_tests)
        for i, (ctrl, tests) in enumerate(dict_ctrl_tests.items()):
            ntests = len(tests)

            # Get the ctrl columns
            colns_ctrl = metadata_df[metadata_df['condition'] == ctrl].index.to_list()
            means_ctrl = df_counts[colns_ctrl].mean(axis=1)
            # Get all scatter plot values
            df_logfc = np.log2(df_counts.div(means_ctrl, axis=0))

            # Plot ctrl scatter
            xsc = [r for c in df_logfc[colns_ctrl].values for r in c]
            jit = plotg['xjit']
            spread_clust = plotg['spread_clust']
            spread_ctrl = plotg['spread_ctrl']
            xjits = np.linspace(-jit, jit, len(colns_ctrl))
            ysc = [
                (k - spread_clust/2 + i*ntests + spread_ctrl + xjits[k])
                for c in df_logfc[colns_ctrl].values
                for k, _ in enumerate(c)
            ]
            ax.scatter(
                xsc, ysc, 
                facecolor=plotg['dotcolor'],
                edgecolor='none', 
                s=plotg['dotsize'], 
                alpha=plotg['dotalpha']
            )

            for j, test in enumerate(tests):
                # bar plot
                colns_tst = metadata_df[metadata_df['condition'] == test].index.to_list()
                yb = np.log2(df_counts[colns_tst].mean(axis=1) / means_ctrl)
                xb = (np.arange(ncontigs) + (j + 1)*spread_clust/ntests 
                      - spread_clust/2 + i*ntests + spread_ctrl)
                ax.barh(xb, yb, plotg['barwidth'], color=plotg['barcolor'])

                # scatter plot
                xsc = [r for c in df_logfc[colns_tst].values for r in c]
                xjits = np.linspace(-jit, jit, len(colns_tst))
                ysc = [
                    (k + (j + 1)*spread_clust/ntests - spread_clust/2 
                     + i*ntests + spread_ctrl + xjits[k])
                    for c in df_logfc[colns_tst].values
                    for k, _ in enumerate(c)
                ]
                ax.scatter(
                    xsc, ysc, 
                    facecolor=plotg['dotcolor'],
                    edgecolor='none', 
                    s=plotg['dotsize'], 
                    alpha=plotg['dotalpha']
                )

                # Load deseq stats 
                pkl_path = f'{stats_dir}/{experiment}_{test}_vs_{ctrl}_stats.pkl'
                try:
                    with open(pkl_path, 'rb') as f:
                        # Expecting a pydeseq2.ds.DeseqStats object inside
                        deseq_object = pickle.load(f)
                except Exception as e:
                    click.echo(f"Failed to load {pkl_path}: {e}", err=True)

        # line at 0 for control
        ax.axvline(x=0, color='k', linestyle='-', linewidth=1)

        # adjust the plot
        spines_invisible = ['left','right']
        for s in spines_invisible:
            ax.spines[s].set_visible(False)
        ax.tick_params(axis='both', labelsize=plotg['ft1']) 
        ax.tick_params(axis='y', length=0, width=0, which='both')
        ax.tick_params(axis='x', direction='in')


    

    click.echo(f"Output saved to {output_path}")

if __name__ == '__main__':
    main()

