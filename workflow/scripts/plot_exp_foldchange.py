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
import re
import glob
import ibis
import yaml
import math
import click
import pickle
import numpy as np
import pandas as pd
from ibis import _, selectors
import matplotlib.pyplot as plt
from collections import defaultdict

def save_fig(bn, exts=['png','pdf'], dpi=500):
    fns_out = [f'{bn}.{ext}' for ext in exts]
    for fn_out in fns_out:
        plt.savefig(fn_out, dpi=dpi, bbox_inches='tight')


@click.command()
@click.option(
    '--counts-path', '-c', 
    required=True, 
    type=click.Path(exists=True, dir_okay=False), 
    help='Path to the counts table Parquet file.'
)
@click.option(
    "--colname_contigs",
    "-n",
    type=str,
    required=True,
    help="Name of the column with contig names to group into clusters.",
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
    '--output-dir', '-o', 
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
def main(counts_path, colname_contigs, clusters_format, metadata_path, stats_dir, config_path, exp_info_path, experiment, output_dir, threads):
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
        
    metadata_df = pd.read_csv(metadata_path, index_col=0)
    
    # Lazy load the parquet file
    counts_table = ibis.read_parquet(counts_path)
    
    # -------------------------------------------------------------------------
    # Step 4: Plot
    # -------------------------------------------------------------------------
    
    exp_config = exp_info[experiment]
    dict_ctrl_tests = exp_config['deseq_ctrl_tests']
    prefix = exp_config['prefix']

    plotg = config['experiment_fold_change_plots']
    plot_e = exp_config['plot_deseq']

    for gene in config['target_genes_for_exp_placement']:
        fig, ax = plt.subplots(figsize=plot_e['figsize'])
        # Load clusters for target gene
        dict_clust_contig = defaultdict(list)
        fn_cl = clusters_format.format(exp=experiment, gene=gene)
        clusters_table = ibis.read_csv(
            fn_cl, table_name="clusters", sep="\t", 
            column_names=["cluster","contig"], header=False
        ).to_pandas()
        for clust, cont in clusters_table.values:
            clust_ = re.sub(prefix + '_','',clust)
            clust_ = re.sub(r'_\d+$','',clust_)
            dict_clust_contig[clust_].append(cont)

        # Get order of clusts based on sum of relative abundances
        colns_sam = metadata_df.index.to_list()
        df_sums = counts_table.aggregate(
            selectors.across(colns_sam, _.sum())
        ).to_pandas()
        df_counts = counts_table.select([colname_contigs] + colns_sam).filter(
            counts_table[colname_contigs].isin(dict_clust_contig)
        ).to_pandas()
        df_pct = df_counts[colns_sam].div(
            df_sums.values, axis=1
        ) * 100
        df_pct['relabund_sum'] = df_pct.sum(axis=1)
        df_pct[colname_contigs] = df_counts[colname_contigs]
        df_pct = df_pct.sort_values(by='relabund_sum', ascending=False)
        ncontigs = df_pct.shape[0]

        # Get each control case
        nctrls = len(dict_ctrl_tests)
        for i, (ctrl, tests) in enumerate(dict_ctrl_tests.items()):
            ntests = len(tests)

            # Get the ctrl columns
            colns_ctrl = metadata_df[
                metadata_df['condition'] == ctrl
            ].index.to_list()
            # means_ctrl = df_counts[colns_ctrl].mean(axis=1)

            # # Get all scatter plot values
            # df_scat = np.log2(df_counts[colns_sam].div(means_ctrl, axis=0))
            df_scat = df_pct

            # bar plot ctrl
            spread_clust = plotg['spread_clust']
            spread_ctrl = plotg['spread_ctrl']
            bar_width = spread_clust*(1 - plotg['bar_gap_frac'])/(ntests + 1)
            yb = df_scat[colns_ctrl].mean(axis=1)
            # yb = np.log2(df_counts[colns_tst].mean(axis=1) / means_ctrl)
            xb = (np.arange(ncontigs) - spread_clust/2 
                  + i*nctrls + spread_ctrl).tolist()
            ax.barh(xb, yb, bar_width, color=plotg['barcolor'])
            yticks = xb
            yticklabels = [ctrl]*len(xb)

            # Plot ctrl scatter
            xsc = [r for c in df_scat[colns_ctrl].values for r in c]
            jit = plotg['xjit_bar_frac']*bar_width
            yjits = np.linspace(-jit, jit, len(colns_ctrl))
            ysc = [
                (k - spread_clust/2 + i*nctrls + spread_ctrl + yjits[l])
                for k, c in enumerate(df_scat[colns_ctrl].values)
                for l, _ in enumerate(c)
            ]
            ax.scatter(
                xsc, ysc, 
                facecolor=plotg['dotcolor'],
                edgecolor='none', 
                s=plotg['dotsize'], 
                alpha=plotg['dotalpha']
            )

            # Get non detected ctrl
            dict_xb_nd = {}
            for x_, ys in zip(xb, df_scat[colns_ctrl].values):
                nd = 0
                for y_ in ys:
                    if y_ == 0:
                        nd += 1
                if nd > 0:
                    dict_xb_nd[x_] = nd


            # Get each test case
            dict_xb_pval = {}
            for j, test in enumerate(tests):

                # bar plot
                colns_tst = metadata_df[
                    metadata_df['condition'] == test
                ].index.to_list()
                yb = df_scat[colns_tst].mean(axis=1)
                # yb = np.log2(df_counts[colns_tst].mean(axis=1) / means_ctrl)
                xb = (np.arange(ncontigs) + (j + 1)*spread_clust/ntests 
                      - spread_clust/2 + i*nctrls + spread_ctrl).tolist()
                ax.barh(xb, yb, bar_width, color=plotg['barcolor'])
                yticks += xb
                yticklabels += [test]*len(xb)

                # scatter plot
                xsc = [r for c in df_scat[colns_tst].values for r in c]
                yjits = np.linspace(-jit, jit, len(colns_tst))
                ysc = [
                    (k + (j + 1)*spread_clust/ntests - spread_clust/2 
                     + i*nctrls + spread_ctrl + yjits[l])
                    for k, c in enumerate(df_scat[colns_tst].values)
                    for l, _ in enumerate(c)
                ]
                ax.scatter(
                    xsc, ysc, 
                    facecolor=plotg['dotcolor'],
                    edgecolor='none', 
                    s=plotg['dotsize'], 
                    alpha=plotg['dotalpha']
                )

                # Get non detected
                for x_, ys in zip(xb, df_scat[colns_tst].values):
                    nd = 0
                    for y_ in ys:
                        if y_ == 0:
                            nd += 1
                    if nd > 0:
                        dict_xb_nd[x_] = nd

                # Deseq stats 
                test_, ctrl_ = [re.sub(r'\s', '_', tc) for tc in [test, ctrl]]
                pkl_path = f'{stats_dir}/{experiment}_{test_}_vs_{ctrl_}_stats.pkl'
                try:
                    with open(pkl_path, 'rb') as f:
                        # Expecting a pydeseq2.ds.DeseqStats object inside
                        ds = pickle.load(f)
                except Exception as e:
                    click.echo(f"Failed to load {pkl_path}: {e}", err=True)
                # get_p_vals
                for x_, cl in zip(xb, df_scat[colname_contigs].values):
                    pval = ds.results_df['pvalue'].get(cl,1)
                    dict_xb_pval[x_] = pval
        
        # Print significant p values to the figure
        xlims = ax.get_xlim()
        xticks = ax.get_xticks()
        xtickrange = xticks[1] - xticks[0]
        for x_, pval in dict_xb_pval.items():
            if (pval is not None) and (pval < 0.05):
                pv = round(pval,4)
                pv = f'={pv}' if pv > 0 else f'<0.00005'
                ax.text(
                    xlims[1] + xtickrange*plotg['pval_adj'], x_, f'p{pv}',
                    fontsize=plotg['ft0'] - 1,
                    va='baseline',
                    ha='right'
                )

        # Print nds to the figure
        if dict_xb_nd:
            for y_, nd in dict_xb_nd.items():
                ax.text(
                    xlims[0] + xtickrange*plotg['nd_shift'],  y_, f'n.d.({nd})',
                    fontsize=plotg['ft0'] - 1,
                    va='baseline',
                    ha='left'
                )

        # # line at 0 for control
        # ax.axvline(x=0, color='k', linestyle='-', linewidth=1)

        # adjust the plot
        spines_invisible = ['right','top']
        for s in spines_invisible:
            ax.spines[s].set_visible(False)
        ax.set_yticks(yticks, labels=yticklabels)
        ax.tick_params(axis='y', labelsize=plotg['ft1']) 
        ax.tick_params(axis='x', labelsize=plotg['ft0'], direction='in') 
        # ax.tick_params(axis='y', length=0, width=0, which='both')
        ax.tick_params(axis='y', pad=plotg['ft1']*1.5)

        # # Set up second axis
        # ax2 = ax.twiny()
        # spines_invisible = ['left','right']
        # for s in spines_invisible:
        #     ax2.spines[s].set_visible(False)
        # ax2.tick_params(axis='x', labelsize=plotg['ft1'], direction='in') 
        # ## Get tick translation
        # # max and min plotted vals
        # yvals = counts_table[colns_sam].values
        # yvals_gt0 = yvals[yvals > 0]
        # mxval = np.max(yvals_gt0)
        # mnval = np.min(yvals_gt0)
        # # max and min exponents used in sci notation
        # expmn = math.floor(np.log10(mnval))
        # expmx = math.floor(np.log10(mxval))
        # tickexps = np.arange(expmn, expmx + 1)
        # # values to plot as ticks for each exponent
        # tvals = np.array([1,5]) if expmx - expmn > 1 else np.array([1,2,4,8])
        # y2ticklabels = []
        # for exp in tickexps:
        #     newticks = (tvals/10**(-exp)).tolist()
        #     y2ticklabels += newticks
        # # Convert to the original ax scale
        # y2ticks = np.log2(y2ticklabels/mean_ctrl)
        # ax2.set_xticks(y2ticks, labels=y2ticklabels, rotation=45, ha='right',rotation_mode="anchor");
        # # # Adjust horizontal position
        # # Now limit to the original axis
        # ax2.set_xlim(ax.get_xlim())
        # # Flip around
        ax.invert_yaxis()
        # ax.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)
        # ax.xaxis.set_ticks_position('top')
        # ax2.xaxis.set_ticks_position('bottom')
        # ax2.xaxis.set_label_position('bottom')

        # Save plot
        os.makedirs(output_dir, exist_ok=True)
        bn = f'{output_dir}/horizontal_bar_{gene}-{experiment}'
        save_fig(bn)
        plt.close()

        click.echo(f"Output saved to {bn}.pdf")
    return

if __name__ == '__main__':
    main()

