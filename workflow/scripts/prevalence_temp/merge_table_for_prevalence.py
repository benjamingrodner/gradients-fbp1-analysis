import pandas as pd
from pathlib import Path
import glob
import re


def load_cluster_map(cluster_glob):
    gene_cl = []
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
                gene_cl.append(parts)

    if not gene_cl:
        die(f"cluster files matching {cluster_glob} contained no usable mappings")

    return pd.DataFrame(gene_cl, columns=["contig","cluster"])


def validate_metadata(df):
    required = {"entry_id", "tax_id", "data_type"}
    missing = required - set(df.columns)
    if missing:
        die(f"metadata missing columns: {', '.join(sorted(missing))}")


def load_metadata(path):
    path = require_file(path, "metadata")
    meta1 = pd.read_csv(path)
    validate_metadata(meta1)
    meta1["entry_id"] = meta1["entry_id"].astype(str)

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


def die(msg, code=2):
    raise SystemExit(f"ERROR: {msg}")


def require_file(path, label):
    path = Path(path)
    if not path.exists():
        die(f"{label} not found: {path}")
    if not path.is_file():
        die(f"{label} is not a file: {path}")
    return path


def load_entry_id_map(path):
    path = require_file(path, "entry-id map")
    df_id = pd.read_csv(path, sep="\t", header=None)
    if df_id.shape[1] != 3:
        die(
            f"entry-id map must have exactly 3 columns, found {df_id.shape[1]} in {path}"
        )
    df_id.columns = ["aa_id", "entry_id", "source_defline"]
    validate_entry_id_map(df_id)
    df_id["entry_id"] = df_id["entry_id"].astype(str)
    return dict(zip(df_id["aa_id"].values, df_id["entry_id"].values))


def merge_write_table(fn_entry_id_map, fns_cluster, regex_clust_gene, fn_meta, fn_out):
    # Load files
    dict_contig_entry = load_entry_id_map(fn_entry_id_map)
    _, dict_entry_tid, dict_entry_dtp = load_metadata(fn_meta)
    df_clust = []
    for fn in fns_cluster:
        df_cl = load_cluster_map(fn)
        gene = re.search(regex_clust_gene, fn).group('gene')
        df_cl['gene'] = gene
        df_clust.append(df_cl)
    df_clust = pd.concat(df_clust, ignore_index=True)

    # merge
    df_merge['entry_id'] = df_clust['contig'].map(dict_contig_entry)
    df_merge["taxid"] = df_clust["entry_id"].map(dict_entry_tid)
    df_merge["data_type"] = df_clust["entry_id"].map(dict_entry_dtp)

    # write
    df_merge.to_csv(fn_out, index=False)

    return


def main():
    # Merge marferret table
    fn_entry_id_map = "/Users/benjamingrodner/work/armbrust/data/prevalence/marferret/marfmmdb_contig_source-fbp1_isips_fre_aldo_fldA.grep.tsv"

    fns_cluster = [
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/fbp1_hmmhitseqs.mmseqs2_cluster.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/aldo_fre_fldA_hmmhitseqs.mmseqs2_cluster.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/isip2_isip3_isip1_hmmhitseqs.mmseqs2_cluster.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/isip2_isip3_isip1_hmmhitseqs.mmseqs2_cluster.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/isip2_isip3_isip1_hmmhitseqs.mmseqs2_cluster.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/aldo_fre_fldA_hmmhitseqs.mmseqs2_cluster.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/clustering/mmseqs2/aldo_fre_fldA_hmmhitseqs.mmseqs2_cluster.tsv",
    ]

    fn_meta = "/Users/benjamingrodner/work/armbrust/data/prevalence/marferret/MarFERReT.v1.1.1.metadata.csv"
    

    return

if __name__ == "__main__":

    main()
