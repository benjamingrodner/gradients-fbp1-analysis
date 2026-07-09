#!/usr/bin/env python3

import argparse
import glob
from collections import defaultdict
from pathlib import Path
import yaml
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from ete4 import NCBITaxa
from matplotlib.gridspec import GridSpec


def die(msg, code=2):
    raise SystemExit(f"ERROR: {msg}")


def save_fig(output_prefix, exts=("png", "pdf"), dpi=500):
    for ext in exts:
        plt.savefig(f"{output_prefix}.{ext}", dpi=dpi, bbox_inches="tight")


def plot_tree_and_heatmap(
    t,
    leaf_order,
    y_coords,
    x_coords,
    reordered_df,
    treelabels,
    output_prefix,
    cfg,
    dpi=500,
):
    fig = plt.figure(figsize=cfg['figsize'])
    gs = GridSpec(1, 2, width_ratios=cfg['width_ratios'], wspace=cfg['wspace'])

    ax_tree = fig.add_subplot(gs[0])
    ax_heatmap = fig.add_subplot(gs[1])

    for node in t.traverse("postorder"):
        if not node.is_root:
            parent = node.up
            ax_tree.plot(
                [x_coords[parent], x_coords[node]],
                [y_coords[node], y_coords[node]],
                color="black",
            )
            sciname = node.props["sci_name"]
            if sciname in treelabels:
                ax_tree.text(
                    x_coords[node],
                    y_coords[node],
                    sciname,
                    fontsize=cfg['ft1'],
                    ha="right",
                    va="bottom",
                )
            ax_tree.plot(
                [x_coords[parent], x_coords[parent]],
                [y_coords[parent], y_coords[node]],
                color="black",
            )

    ax_tree.set_ylim(-0.5, len(leaf_order) - 0.5)
    ax_tree.axis("off")
    ax_tree.invert_yaxis()

    mx = np.max(reordered_df.values) + 1
    bounds = np.arange(mx)
    current_cmap = plt.cm.get_cmap(cfg['cmap']).copy()
    colors = [current_cmap(b / mx) for b in bounds]
    colors[0] = (0.5, 0.5, 0.5)
    custom_cmap = mcolors.ListedColormap(colors)
    custom_norm = mcolors.BoundaryNorm(bounds, custom_cmap.N)
    midpoints = [(bounds[i] + bounds[i + 1]) / 2.0 for i in range(len(bounds) - 1)]

    sns.heatmap(
        reordered_df,
        cmap=custom_cmap,
        norm=custom_norm,
        ax=ax_heatmap,
        cbar_kws={"ticks": midpoints, "shrink":1.25},
        xticklabels=reordered_df.columns
    )
    for i in range(reordered_df.shape[1] + 1):
        ax_heatmap.axvline(i, color="white", lw=cfg['heatmap_whitespace_width'])

    cbar = ax_heatmap.collections[0].colorbar
    cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(midpoints))
    cbar.set_ticklabels(bounds[:-1])
    cbar.ax.yaxis.set_minor_locator(ticker.NullLocator())
    cbar.ax.tick_params(labelsize=cfg['ft0'])
    cbar.set_label(cfg['cbar_lab'], fontsize=cfg['ft1'])
    
    ax_heatmap.get_yaxis().set_visible(False)
    ax_heatmap.tick_params(axis='both', labelsize=cfg['ft1'])
    for label in ax_heatmap.get_xticklabels():
        # label.set_ha('right')
        # label.set_va('top')
        label.set_rotation(45)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{output_prefix}.{ext}", dpi=dpi, bbox_inches="tight")
    plt.show()


def build_tree_layout(t, prefix):
    leaf_order = []
    y_coords = {}
    y_counter = 0

    for node in t.traverse("postorder"):
        if node.name.startswith(prefix):
            leaf_order.append(node.name)
            y_coords[node] = y_counter
            y_counter += 1
        else:
            children_y = [y_coords[child] for child in node.children if child in y_coords]
            if children_y:
                y_coords[node] = sum(children_y) / len(children_y)

    x_coords = {}
    max_depth = 1.0
    for node in t.traverse(strategy="postorder"):
        if node.is_leaf:
            x_coords[node] = max_depth
        else:
            child_distances = [
                x_coords[child] - (child.dist if child.dist is not None else 1.0)
                for child in node.children
            ]
            x_coords[node] = min(child_distances)

    if t not in x_coords:
        x_coords[t] = min(
            [x_coords[child] - (child.dist if child.dist is not None else 1.0) for child in t.children]
        )

    return leaf_order, y_coords, x_coords

def build_presence_df(leaf_order, feature_to_values, prefix):
    data = {}
    entry_ids = [int(e.lstrip(prefix)) for e in leaf_order]

    for feature_name, value_map in feature_to_values.items():
        data[feature_name] = [value_map.get(entry_id, 0) for entry_id in entry_ids]

    df = pd.DataFrame(data, index=leaf_order)
    return df.reindex(leaf_order)
# def build_presence_df(leaf_order, dict_entry_ncopies, prefix, gene):
#     vals = [dict_entry_ncopies[int(e.lstrip(prefix))] for e in leaf_order]
#     df = pd.DataFrame({gene: vals}, index=leaf_order)
#     return df.reindex(leaf_order)


def require_file(path, label):
    path = Path(path)
    if not path.exists():
        die(f"{label} not found: {path}")
    if not path.is_file():
        die(f"{label} is not a file: {path}")
    return path


def load_cluster_map(cluster_glob):
    dict_gene_cl = {}
    fns_clusters = sorted(glob.glob(str(cluster_glob)))
    if not fns_clusters:
        die(f"no cluster files found matching pattern: {cluster_glob}")

    for fn in fns_clusters:
        fn = require_file(fn, "cluster file")
        with open(fn, "r") as f:
            for lineno, line in enumerate(f, start=1):
                parts = line.split()
                if len(parts) != 2:
                    die(
                        f"malformed cluster file {fn} line {lineno}: "
                        f"expected 2 fields, found {len(parts)}"
                    )
                cl, gene = parts
                dict_gene_cl[gene] = cl

    if not dict_gene_cl:
        die(f"cluster files matching {cluster_glob} contained no usable mappings")

    return dict_gene_cl


def load_contig_names(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"contig list not found: {path}")
    if not path.is_file():
        raise IsADirectoryError(f"contig list is not a file: {path}")

    contigs = []
    seen = set()

    with open(path, "r") as f:
        for lineno, raw in enumerate(f, start=1):
            line = raw.strip()
            if not line:
                continue
            if " " in line or "\t" in line:
                raise ValueError(
                    f"malformed contig name at line {lineno} in {path}: "
                    f"expected one name per line"
                )
            if line in seen:
                raise ValueError(
                    f"duplicate contig name at line {lineno} in {path}: {line}"
                )
            seen.add(line)
            contigs.append(line)

    if not contigs:
        raise ValueError(f"contig list is empty: {path}")

    return contigs


def get_dict_entry_count(contig_name_file, cluster_glob, dict_contig_entry, meta1):
    contig_names = load_contig_names(contig_name_file)
    dict_contig_cl = load_cluster_map(cluster_glob)
    # map isolate to list of clusters
    dict_entry_clusts = defaultdict(list)
    skipped = []
    # TODO: make this cleaner...Don't rely on .get() 
    for contig in contig_names:
        cl = dict_contig_cl.get(contig)
        entry = dict_contig_entry.get(contig)
        if (cl is not None) & (entry is not None):
            dict_entry_clusts[entry].append(cl)
        else:
            skipped.append(contig)
    if len(dict_entry_clusts) == 0:
        raise ValueError(f"No contigs from hmmsearch in clustering \n first ten: {skipped[:10]}")
    # only unique set of clusters
    for entry, clusts in dict_entry_clusts.items():
        dict_entry_clusts[entry] = list(set(clusts))
    # map isolate to number of copies
    dict_tid_entries = defaultdict(list)
    dict_entry_ncopies = {}
    for _, row in meta1.iterrows():
        entry = row["entry_id"]
        dict_tid_entries[row["tax_id"]].append(entry)
        clusts = dict_entry_clusts.get(entry)
        dict_entry_ncopies[entry] = len(clusts) if clusts is not None else 0
    return dict_entry_ncopies


def load_config(path):
    path = require_file(path, "config")
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        die("config must be a YAML mapping")

    required = {"dtp", "ttnames", "treelabels"}
    missing = required - set(cfg)
    if missing:
        die(f"config missing keys: {', '.join(sorted(missing))}")

    if not isinstance(cfg["ttnames"], list) or not all(isinstance(x, str) for x in cfg["ttnames"]):
        die("config key 'ttnames' must be a list of strings")

    if not isinstance(cfg["treelabels"], list) or not all(isinstance(x, str) for x in cfg["treelabels"]):
        die("config key 'treelabels' must be a list of strings")

    if not isinstance(cfg["dtp"], str) or not cfg["dtp"]:
        die("config key 'dtp' must be a non-empty string")

    return cfg


def validate_metadata(df):
    required = {"entry_id", "tax_id", "data_type"}
    missing = required - set(df.columns)
    if missing:
        die(f"metadata missing columns: {', '.join(sorted(missing))}")


def load_metadata(path):
    path = require_file(path, "metadata")
    meta1 = pd.read_csv(path)
    validate_metadata(meta1)

    if meta1.empty:
        die(f"metadata file is empty: {path}")

    if meta1["entry_id"].isna().any():
        die("metadata contains missing entry_id values")
    if meta1["tax_id"].isna().any():
        die("metadata contains missing tax_id values")
    if meta1["data_type"].isna().any():
        die("metadata contains missing data_type values")

    dict_entry_tid = dict(zip(meta1["entry_id"].values, meta1["tax_id"].values))
    dict_entry_dtp = dict(zip(meta1["entry_id"].values, meta1["data_type"].values))
    return meta1, dict_entry_tid, dict_entry_dtp


def validate_entry_id_map(df):
    required = {"aa_id", "entry_id", "source_defline"}
    missing = required - set(df.columns)
    if missing:
        die(f"entry-id map missing columns: {', '.join(sorted(missing))}")


def load_entry_id_map(path):
    path = require_file(path, "entry-id map")
    df_id = pd.read_csv(path, sep="\t", header=None)
    if df_id.shape[1] != 3:
        die(f"entry-id map must have exactly 3 columns, found {df_id.shape[1]} in {path}")
    df_id.columns = ["aa_id", "entry_id", "source_defline"]
    validate_entry_id_map(df_id)
    return dict(zip(df_id["aa_id"].values, df_id["entry_id"].values))


def parse_args():
    p = argparse.ArgumentParser(description="Build FBP1 presence table")
    p.add_argument("--entry-id-map", required=True, type=Path)
    p.add_argument("--genes", required=True, nargs="+", help="Gene names to include")
    p.add_argument("--contig-name-files", required=True, nargs="+", type=Path, help="Same order as '--genes'")  
    p.add_argument("--cluster-globs", required=True, nargs="+", type=Path, help="Same order as '--genes'")  
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--output-prefix", required=True, type=Path)
    p.add_argument(
        "--config",
        required=True,
        type=Path,
        help="YAML config containing dtp, ttnames, and treelabels.",
    )
    p.add_argument("--dpi", type=int, default=500)
    return p.parse_args()


def main():
    args = parse_args()
    ncbi = NCBITaxa()

    ## File Loading

    dict_contig_entry = load_entry_id_map(args.entry_id_map)
    meta1, dict_entry_tid, dict_entry_dtp = load_metadata(args.metadata)

    cfg = load_config(args.config)['prevalence']


    ## Useful dicts

    # Major taxonomic groups to plot
    ttnames = cfg["ttnames"]
    ttn_trans = ncbi.get_name_translator(ttnames)
    ttids = [ttn_trans[n][0] for n in ttnames]
   
    # Map different data types to taxids
    dict_dtp_tids = defaultdict(list)
    for entry, tid in dict_entry_tid.items():
        lin = ncbi.get_lineage(tid)
        if any([t in lin for t in ttids]):
            typ = dict_entry_dtp[entry]
            dict_dtp_tids[typ].append(tid)

    # map taxid to list of isolates
    dict_tid_entries = defaultdict(list)
    for _, row in meta1.iterrows():
        entry = row["entry_id"]
        dict_tid_entries[row["tax_id"]].append(entry)

    # map gene to isolate to number of copies
    feature_to_values = {
        gene: get_dict_entry_count(contig_name_file, cluster_glob, dict_contig_entry, meta1)
        for gene, contig_name_file, cluster_glob 
        in zip(args.genes, args.contig_name_files, args.cluster_globs)
    }

    ## Tree building

    dtp = cfg["dtp"]
    prefix = "marf_"
    t = ncbi.get_topology(dict_dtp_tids[dtp])

    for n in t.traverse():
        tid = n.name
        entries = dict_tid_entries.get(int(tid))
        if entries is not None:
            n.name = prefix + str(entries[0])
            sciname = n.props["sci_name"]
            if len(entries) > 1:
                for m in entries[1:]:
                    n_new = n.add_sister(name=prefix + str(m))
                    n_new.add_props(sci_name=sciname)

    ## Plotting

    leaf_order, y_coords, x_coords = build_tree_layout(t, prefix)
    reordered_df = build_presence_df(leaf_order, feature_to_values, prefix)

    treelabels = cfg["treelabels"]

    plot_tree_and_heatmap(
        t=t,
        leaf_order=leaf_order,
        y_coords=y_coords,
        x_coords=x_coords,
        reordered_df=reordered_df,
        treelabels=treelabels,
        output_prefix=args.output_prefix,
        cfg=cfg,
        dpi=args.dpi,
    )

if __name__ == "__main__":
    main()



