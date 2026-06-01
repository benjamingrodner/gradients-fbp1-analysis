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

def save_fig(fig, bn, exts=['png','pdf'], dpi=500):
    fns_out = [f'{bn}.{ext}' for ext in exts]
    for fn_out in fns_out:
        fig.savefig(fn_out, dpi=dpi, bbox_inches='tight')


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
    # prefix = exp_config['prefix']
    re_subs = exp_config.get('re_subs_fasta2counts')

    plotg = config['experiment_fold_change_plots']
    plot_e = exp_config['plot_deseq']

    for gene in config['target_genes_for_exp_placement']:
        # Load clusters for target gene
        dict_clust_contig = defaultdict(list)
        fn_cl = clusters_format.format(exp=experiment, gene=gene)
        clusters_table = ibis.read_csv(
            fn_cl, table_name="clusters", sep="\t", 
            column_names=["cluster","contig"], header=False
        ).to_pandas()
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
            dict_clust_contig[cs[0]].append(cs[1])

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
        df_pct.index = df_counts[colname_contigs]
        df_pct = df_pct.sort_values(by='relabund_sum', ascending=False)
        ncontigs = df_pct.shape[0]
        nctrls = len(dict_ctrl_tests)

        # # info for stacked plot
        # dict_ctrl_stack = {
        #     "grouping": {},
        #     "scat": [[],[]],
        #     "pvals": {},
        #     "nds": {},
        #     "yticks": [[],[]],
        # }

        # iterate through contigs
        for h, df_scat in df_pct.iterrows():
            # clust = df_scat[colname_contigs]
            clust = h
            for i, (ctrl, tests) in enumerate(dict_ctrl_tests.items()):
                ntests = len(tests)
                
                # Percent of transcripts plot
                fig, ax = plt.subplots(figsize=plot_e['figsize'])
                # Log fold change plot
                fig1, ax1 = plt.subplots(figsize=plot_e['figsize'])
                

                # Get the ctrl columns
                colns_ctrl = metadata_df[
                    metadata_df['condition'] == ctrl
                ].index.to_list()
                # if h == 0:
                #     for coln in colns_ctrl:
                #         dict_ctrl_stack[ctrl]['grouping'][coln] = ctrl

                # means_ctrl = df_counts[colns_ctrl].mean(axis=1)

                # # Get all scatter plot values
                # df_scat = np.log2(df_counts[colns_sam].div(means_ctrl, axis=0))
                # df_scat = df_pct

                # bar plot ctrl
                # spread_clust = plotg['spread_clust']
                # spread_ctrl = plotg['spread_ctrl']
                # bar_width = spread_clust*(1 - plotg['bar_gap_frac'])/(ntests + 1)
                bar_width = plotg['bar_width']
                yb = df_scat.loc[colns_ctrl].mean()
                # yb = np.log2(df_counts[colns_tst].mean(axis=1) / means_ctrl)
                # xb = (np.arange(ncontigs) - spread_clust/2 
                #     + i*nctrls + spread_ctrl).tolist()
                x_ = 0
                xb = [x_]
                ax.barh(xb, yb, bar_width, color=plotg['barcolor'])
                yticks = xb
                yticklabels = [ctrl]*len(xb)

                # Plot ctrl scatter
                xsc = df_scat.loc[colns_ctrl].values.tolist()
                # xsc = [r for c in df_scat.loc[colns_ctrl].values for r in c]
                jit = plotg['xjit_bar_frac']*bar_width / 2
                yjits = np.linspace(-jit, jit, len(colns_ctrl))
                ysc = [yjits[l] for l in range(len(xsc))]
                # ysc = [
                #     (k - spread_clust/2 + i*nctrls + spread_ctrl + yjits[l])
                #     for k, c in enumerate(df_scat.loc[colns_ctrl].values)
                #     for l, _ in enumerate(c)
                # ]
                ax.scatter(
                    xsc, ysc, 
                    facecolor=plotg['dotcolor'],
                    edgecolor='none', 
                    s=plotg['dotsize'], 
                    alpha=plotg['dotalpha']
                )

                # Get non detected ctrl
                dict_xb_nd = {}
                # for ys in df_scat.loc[colns_ctrl].values:
                nd = 0
                for y_ in df_scat.loc[colns_ctrl].values:
                    if y_ == 0:
                        nd += 1
                if nd > 0:
                    dict_xb_nd[x_] = nd


                # Get each test case
                dict_xb_pval = {}
                for j, test in enumerate(tests):
                    x_ = j+1
                    # bar plot
                    colns_tst = metadata_df[
                        metadata_df['condition'] == test
                    ].index.to_list()
                    # if h == 0:
                    #     for coln in colns_tst:
                    #         dict_ctrl_stack[ctrl]['grouping'][coln] = test
                    yb = df_scat.loc[colns_tst].mean()
                    xb = [x_]
                    # yb = np.log2(df_counts[colns_tst].mean(axis=1) / means_ctrl)
                    # xb = (np.arange(ncontigs) + (j + 1)*spread_clust/ntests 
                    #     - spread_clust/2 + i*nctrls + spread_ctrl).tolist()
                    ax.barh(xb, yb, bar_width, color=plotg['barcolor'])
                    yticks += xb
                    yticklabels += [test]*len(xb)


                    # scatter plot
                    # xsc = [r for c in df_scat[colns_tst].values for r in c]
                    xsc = df_scat.loc[colns_tst].values.tolist()
                    yjits = np.linspace(-jit, jit, len(colns_tst))
                    ysc = [x_ + yjits[l] for l in range(len(xsc))]
                    # ysc = [
                    #     (k + (j + 1)*spread_clust/ntests - spread_clust/2 
                    #     + i*nctrls + spread_ctrl + yjits[l])
                    #     for k, c in enumerate(df_scat[colns_tst].values)
                    #     for l, _ in enumerate(c)
                    # ]
                    ax.scatter(
                        xsc, ysc, 
                        facecolor=plotg['dotcolor'],
                        edgecolor='none', 
                        s=plotg['dotsize'], 
                        alpha=plotg['dotalpha']
                    )

                    # Get non detected
                    # for ys in zip(xb, df_scat.loc[colns_tst].values):
                    nd = 0
                    for y_ in df_scat.loc[colns_tst].values:
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

                    # Shrunk lfc bar plot
                    if not ds.shrunk_LFCs:
                        raise ValueError(f"LFCs are not shrunk for comparison {test} vs {ctrl} in file {pkl_path}")
                    else:
                        ybl = ds.results_df['log2FoldChange'].get(clust,1)
                        yblerr = ds.results_df['lfcSE'].get(clust,1)
                        ax1.barh(
                            xb, ybl, bar_width, 
                            xerr=yblerr, 
                            error_kw={
                                'elinewidth': plotg['elinewidth'], 
                                'capsize': plotg['ecapsize']
                            }, 
                            color=plotg['barcolor']
                        )

                    # get_p_vals
                    # for x_, cl in zip(xb, df_scat[colname_contigs].values):
                    pval = ds.results_df['pvalue'].get(clust,1)
                    dict_xb_pval[x_] = pval

                # Print significant p values to the figures
                for ax_ in [ax, ax1]:
                    xlims = ax_.get_xlim()
                    xticks = ax_.get_xticks()
                    xtickrange = xticks[1] - xticks[0]
                    for y_, pval in dict_xb_pval.items():
                        if (pval is not None) and (pval < 0.05):
                            pv = round(pval,4)
                            pv = f'={pv}' if pv > 0 else f'<0.00005'
                            txt = f'p{pv}'
                            # txt = f'p{pv}' if ax_ == ax1 else '*'
                            ax_.text(
                                xlims[1] + xtickrange*plotg['pval_adj'], y_, txt,
                                fontsize=plotg['ft0'] - 1,
                                va='center',
                                ha='right'
                            )

                    # Print nds to the figure
                    if dict_xb_nd:
                        for y_, nd in dict_xb_nd.items():
                            ax_.text(
                                xlims[0] + xtickrange*plotg['nd_shift'],  y_, f'n.d.({nd})',
                                fontsize=plotg['ft0'] - 1,
                                va='center',
                                ha='left'
                            )


                    # adjust the plots
                    ax_.set_yticks(yticks, labels=yticklabels)
                    ax_.tick_params(axis='y', labelsize=plotg['ft1']) 
                    ax_.tick_params(axis='x', labelsize=plotg['ft0'], direction='in') 
                    ax_.tick_params(axis='y', length=0, width=0, which='both')
                    ax_.tick_params(axis='y', pad=plotg['ft1']*1.5)
                    # Flip around
                    ax_.invert_yaxis()

                # line at 0 for lfc control
                lw = ax.spines['bottom'].get_linewidth()
                ax1.axvline(x=0, color='k', linestyle='-', linewidth=lw)
                xlims = ax1.get_xlim()
                xticks = ax1.get_xticks()
                xtickrange = xticks[1] - xticks[0]
                adj = 0.01*xtickrange
                if xlims[0] == 0:
                    ax1.set_xlim(xlims[0]-adj, xlims[1])
                elif xlims[1] == 0:
                    ax1.set_xlim(xlims[0], xlims[1]+adj)
                
                # # labels
                # ax.set_xlabel('Estimated percent of transcripts', fontsize=plotg['ft1'])
                # ax1.set_xlabel(f'Log2(Fold Change vs {ctrl})', fontsize=plotg['ft1'])

                # remove spines
                spines_invisible = {ax:['right','top'], ax1:['left','right','top']}
                for ax_, si in spines_invisible.items():
                    for s in si:
                        ax_.spines[s].set_visible(False)
  
                # Save plot
                for tp, fig_ in zip(['pct','lfc'],[fig, fig1]):
                    od = f'{output_dir}/{tp}'
                    os.makedirs(od, exist_ok=True)
                    ctrl_ = re.sub('_','',ctrl)
                    bn = f'{od}/horizontal_bar_{tp}-{gene}-{experiment}_{clust}_{ctrl_}'
                    plt.figure(fig_)
                    save_fig(fig_, bn)
                    plt.close(fig_)

                    click.echo(f"Output saved to {bn}.pdf")

        # Plot stacked bar
        for i, (ctrl, tests) in enumerate(dict_ctrl_tests.items()):
            # Get plot info
            colns_ctrl = metadata_df[
                metadata_df['condition'] == ctrl
            ].index.to_list()
            # Get values
            gene_clust_name = fn_cl
            # df_countsg = counts_table.select([colname_contigs] + colns_sam).filter(
            #     counts_table[colname_contigs] == gene_clust_name
            # ).to_pandas()
            df_countsg = counts_table.select([colname_contigs] + colns_sam).filter(
                counts_table[colname_contigs] == gene_clust_name
            ).to_pandas()
            df_countsg.index = df_countsg[colname_contigs]
            df_pctg = df_countsg[colns_sam].div(
                df_sums.values, axis=1
            ) * 100
            df_scat = df_pctg
            # for grouping
            dict_coln_xb = {c: 0 for c in colns_ctrl}
            # scatter
            coln_order = colns_ctrl
            ysc2 = np.linspace(-jit, jit, len(colns_ctrl)).tolist()
            # bar
            yticks = [0]
            yticklabels = [ctrl]
            dict_xb_nd = {}
            dict_xb_pval = {}
            for j, test in enumerate(tests):
                x_ = j+1
                colns_tst = metadata_df[
                    metadata_df['condition'] == test
                ].index.to_list()
                for c in colns_tst:
                    dict_coln_xb[c] = x_
                yticks.append(x_)
                yticklabels.append(test)
                coln_order += colns_tst
                yjits = np.linspace(-jit, jit, len(colns_tst))
                ysc2 += [x_ + l for l in yjits]
                # Get non detected
                nd = 0
                for y_ in df_scat.loc[:,colns_tst].values[0]:
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
                pval = ds.results_df['pvalue'].get(gene_clust_name,1)
                dict_xb_pval[x_] = pval

            fig2, ax2 = plt.subplots(figsize=plot_e['figsize'])
            
            # Bar
            df_stack = df_pct.T.groupby(dict_coln_xb).mean()
            cmap = plot_e.get('stacked_cmap')
            cmap = cmap if cmap is not None else 'tab20'
            tmp = df_stack.plot(ax=ax2, kind='barh', 
                                stacked=True, cmap=cmap,
                                legend=True
                                )
            
            # scatter
            xsc2 = df_scat[coln_order].values
            ax2.scatter(xsc2, ysc2, color='k', s=plotg['dotsize'])
            
            fig_, ax_ = fig2, ax2
            # Plot pvals
            xlims = ax_.get_xlim()
            xticks = ax_.get_xticks()
            xtickrange = xticks[1] - xticks[0]
            for y_, pval in dict_xb_pval.items():
                if (pval is not None) and (pval < 0.05):
                    pv = round(pval,4)
                    pv = f'={pv}' if pv > 0 else f'<0.00005'
                    txt = f'p{pv}'
                    # txt = f'p{pv}' if ax_ == ax1 else '*'
                    ax_.text(
                        xlims[1] + xtickrange*plotg['pval_adj'], y_, txt,
                        fontsize=plotg['ft0'] - 1,
                        va='center',
                        ha='right'
                    )

            # Print nds to the figure
            if dict_xb_nd:
                for y_, nd in dict_xb_nd.items():
                    ax_.text(
                        xlims[0] + xtickrange*plotg['nd_shift'],  y_, f'n.d.({nd})',
                        fontsize=plotg['ft0'] - 1,
                        va='center',
                        ha='left'
                    )

            # adjust the plots
            ax_.set_yticks(yticks, labels=yticklabels)
            ax_.tick_params(axis='y', labelsize=plotg['ft1']) 
            ax_.tick_params(axis='x', labelsize=plotg['ft0'], direction='in') 
            ax_.tick_params(axis='y', length=0, width=0, which='both')
            ax_.tick_params(axis='y', pad=plotg['ft1']*1.5)
            # Flip around
            ax_.invert_yaxis()

            # remove spines
            for s in ['top','right']:
                ax_.spines[s].set_visible(False)

            # save plot
            od = f'{output_dir}/stacked'
            os.makedirs(od, exist_ok=True)
            ctrl_ = re.sub('_','',ctrl)
            bn = f'{od}/horizontal_bar_stacked-{gene}-{experiment}_{ctrl_}'
            plt.figure(fig_)
            save_fig(fig_, bn)
            plt.close(fig_)



    return

if __name__ == '__main__':
    main()

