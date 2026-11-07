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


def require_file(path, label):
    path = Path(path)
    if not path.exists():
        die(f"{label} not found: {path}")
    if not path.is_file():
        die(f"{label} is not a file: {path}")
    return path


def save_fig(output_prefix, exts=("png", "pdf"), dpi=500):
    for ext in exts:
        plt.savefig(f"{output_prefix}.{ext}", dpi=dpi, bbox_inches="tight")


def plot_tree_and_heatmap(
    t,
    leaf_order,
    y_coords,
    x_coords,
    reordered_df,
    output_prefix,
    cfg,
    dict_node_marker,
    dict_node_markercolor,
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
                lw=cfg["tree_lw"],
            )
            sciname = node.props["sci_name"]
            bl = sciname in cfg["treelabels"]
            if cfg['label_leaves']:
                bl = node.is_leaf
            if bl:
                ax_tree.text(
                    x_coords[node],
                    y_coords[node],
                    sciname,
                    fontsize=cfg['ft1'],
                    ha="right",
                    va="bottom",
                )
            marker = dict_node_marker.get(node.name)
            if marker is not None:
                ax_tree.plot(
                    x_coords[node],
                    y_coords[node],
                    marker=marker,
                    markersize=cfg['markersize'],
                    color=dict_node_markercolor[node.name],
                )
            ax_tree.plot(
                [x_coords[parent], x_coords[parent]],
                [y_coords[parent], y_coords[node]],
                color="black",
                lw=cfg['tree_lw'],
            )

    ax_tree.set_ylim(-0.5, len(leaf_order) - 0.5)
    ax_tree.axis("off")
    ax_tree.invert_yaxis()

    if cfg['max_count'] is not None:
        reordered_df = reordered_df.clip(upper=cfg['max_count']+1)

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
        cbar_kws={"ticks": midpoints, 
                  "shrink": cfg['cbar_shrink'],
                  "aspect": cfg['cbar_aspect'],
                  "fraction": cfg['cbar_fraction']
                  },
        xticklabels=reordered_df.columns
    )
    for i in range(reordered_df.shape[1] + 1):
        ax_heatmap.axvline(i, color=colors[0], lw=cfg['heatmap_whitespace_width'])

    cbar = ax_heatmap.collections[0].colorbar
    cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(midpoints))
    tlab = [str(l) for l in bounds[:-1]]
    tlab[-1] = f"≥{tlab[-1]}"
    cbar.set_ticklabels(tlab)
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
    dout = output_prefix.parent
    if not dout.exists():
        dout.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        plt.savefig(f"{output_prefix}.{ext}", dpi=cfg['dpi'], bbox_inches="tight")
    plt.show()
    plt.close()
    return


def get_marker_dicts(t, df, prefix, cfg):
    dict_node_marker, dict_node_markercolor = {}, {}
    for node in t.traverse():
        nm = node.name
        if nm.startswith(prefix):
            entry = nm.lstrip(prefix)
            dtp = df.loc[df["entry_id"] == entry, "data_type"].unique()[0]
            source = df.loc[df["entry_id"] == entry, "source"].unique()[0]
            dict_node_marker[nm] = cfg["dict_dtp_marker"].get(dtp)
            dict_node_markercolor[nm] = cfg["dict_source_markercolor"].get(source)
    return dict_node_marker, dict_node_markercolor


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
    entry_ids = [e.lstrip(prefix) for e in leaf_order]

    for feature_name, value_map in feature_to_values.items():
        data[feature_name] = [value_map.get(entry_id, 0) for entry_id in entry_ids]

    df = pd.DataFrame(data, index=leaf_order)
    return df.reindex(leaf_order)


def build_tree(dict_tid_entries, dict_tidnew_tidold, ncbi, prefix):
    t = ncbi.get_topology(list(dict_tid_entries.keys()))
    for n in t.traverse():
        tid = int(n.name)
        # tidold = dict_tidnew_tidold.get(int(tid))
        # if tidold is not None:
        #     tid = tidold
        entries = dict_tid_entries.get(int(tid))
        nodes = []
        if entries is not None:
            nodes.append(n)
            n.name = prefix + str(entries[0])
            sciname = n.props["sci_name"]
            if len(entries) > 1:
                for m in entries[1:]:
                    n_new = n.add_sister(name=prefix + str(m))
                    n_new.add_props(sci_name=sciname)
                    nodes.append(n_new)
    return t


def get_dict_tid_entries(df, feature_to_values, cfg, ncbi, ttids):
    bool_source = np.ones(df.shape[0]).astype(bool)
    if not cfg["source"] == "all":
        bool_source = df["source"] == cfg["source"]

    # counts of target gene
    dict_entry_count = feature_to_values[cfg["count_min_gene"]]
    df["count"] = df["entry_id"].map(dict_entry_count).fillna(0)

    bool_dtp = np.ones(df.shape[0]).astype(bool)
    # custom dtype and count minimum
    if not cfg["dtp"] == 'all':
        bool_dtp = df["data_type"] == cfg["dtp"]
        bool_count = df["count"] >= cfg["count_min"]
    # if merging then accept all sags and genomes, but only TSAs with detection
    else:
        bool_count = np.ones(df.shape[0]).astype(bool)
        bool_tsa = df["data_type"] == "TSA"
        bool_count[bool_tsa] = df[bool_tsa]["count"] >= 1

    arr_sub = (
        df.loc[bool_source & bool_count & bool_dtp, ["entry_id", "taxid"]]
        .drop_duplicates()
        .values
    )
    tids = []
    entries = []
    for entry, tid in arr_sub:
        lin = ncbi.get_lineage(tid)
        if any([t in lin for t in ttids]):
            entries.append(entry)
            tids.append(int(tid))

    # map taxid to list of isolates
    dict_tid_entries = defaultdict(list)
    for entry, tid in zip(entries, tids):
        dict_tid_entries[tid].append(entry)
    return dict_tid_entries


def get_dict_entry_count(df):
    df_sub = df[['entry_id','cluster']].drop_duplicates()
    return df_sub.groupby('entry_id').size().to_dict()


def load_config(path):
    path = require_file(path, "config")
    with open(path, "r") as f:
        cfg = yaml.safe_load(f)

    if not isinstance(cfg, dict):
        die("config must be a YAML mapping")

    if 'prevalence' not in cfg:
        die(f"config missing keys: prevalence")
    cfg = cfg['prevalence']

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


def validate_table(df, path):
    required = {"contig", "cluster", "entry_id", "gene", "taxid", "data_type", "source"}
    missing = required - set(df.columns)
    if missing:
        die(f"metadata missing columns: {', '.join(sorted(missing))}")
    if df.empty:
        die(f"Prevalence table file is empty: {path}")

    for c in required:
        if df[c].isna().any():
            die(f"Prevalence table contains missing {c} values")
    return


def load_table(path):
    path = require_file(path, "Prevalence table")
    meta1 = pd.read_csv(path)
    validate_table(meta1, path)
    meta1["entry_id"] = meta1["entry_id"].astype(str)
    meta1["taxid"] = meta1["taxid"].astype(int)
    return meta1

def parse_args():
    p = argparse.ArgumentParser(description="Build FBP1 presence table")
    p.add_argument("--prevalence-table", required=True, type=Path)
    p.add_argument("--genes", required=True, nargs="+", help="Gene names to include")
    p.add_argument("--output-prefix", required=True, type=Path)
    p.add_argument(
        "--config",
        required=True,
        type=Path,
        help="YAML config containing dtp, ttnames, and treelabels.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    ncbi = NCBITaxa()

    ## File Loading
    df = load_table(args.prevalence_table)
    cfg = load_config(args.config)

    ## Useful dicts
    # map gene to isolate to number of copies
    feature_to_values = {
        gene: get_dict_entry_count(df[df['gene'] == gene])
        for gene in args.genes
    }

    # Major taxonomic groups to plot
    ttnames = cfg["ttnames"]
    ttn_trans = ncbi.get_name_translator(ttnames)
    ttids = [ttn_trans[n][0] for n in ttnames]

    # Get sublist of taxa and their entries
    dict_tid_entries = get_dict_tid_entries(df, feature_to_values, cfg, ncbi, ttids)

    # get mapping of changed ids
    tids_old = df["taxid"].unique().astype(int)
    nms = [ncbi.get_taxid_translator([tid])[tid] for tid in tids_old]
    dict_nm = ncbi.get_name_translator(nms)
    tids_new = [dict_nm[nm][0] for nm in nms]
    dict_tidnew_tidold = dict(zip(tids_new, tids_old))

    ## Tree building
    prefix = "marf_"
    t = build_tree(dict_tid_entries, dict_tidnew_tidold, ncbi, prefix)
    dict_node_marker, dict_node_markercolor = get_marker_dicts(t, df, prefix, cfg)
    ## Plotting
    leaf_order, y_coords, x_coords = build_tree_layout(t, prefix)
    reordered_df = build_presence_df(leaf_order, feature_to_values, prefix)

    plot_tree_and_heatmap(
        t=t,
        leaf_order=leaf_order,
        y_coords=y_coords,
        x_coords=x_coords,
        reordered_df=reordered_df,
        output_prefix=args.output_prefix,
        cfg=cfg,
        dict_node_marker=dict_node_marker,
        dict_node_markercolor=dict_node_markercolor,
    )

if __name__ == "__main__":
    main()
