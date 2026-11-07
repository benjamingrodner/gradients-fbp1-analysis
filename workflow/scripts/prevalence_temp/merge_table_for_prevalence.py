import pandas as pd
from pathlib import Path
from ete4 import NCBITaxa
import glob
import re


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

    # get mapping of changed ids
    ncbi = NCBITaxa()
    tids_old = meta1["tax_id"].unique().astype(int)
    nms = [ncbi.get_taxid_translator([tid])[tid] for tid in tids_old]
    dict_nm = ncbi.get_name_translator(nms)
    tids_new = [dict_nm[nm][0] for nm in nms]
    dict_tidold_tidnew = dict(zip(tids_old, tids_new))

    meta1['tax_id'] = meta1['tax_id'].map(dict_tidold_tidnew)
    meta1["entry_id"] = meta1["entry_id"].astype(str)

    dict_entry_tid = dict(zip(meta1["entry_id"].values, meta1["tax_id"].values))
    dict_entry_dtp = dict(zip(meta1["entry_id"].values, meta1["data_type"].values))
    return meta1, dict_entry_tid, dict_entry_dtp


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


def merge_write_table(fn_entry_id_map, genes, fns_names, fns_cluster, fn_meta, source, fn_out):
    # Load files
    dict_contig_entry = load_entry_id_map(fn_entry_id_map)
    df = []
    for gene, fn_n, fn_c in zip(genes, fns_names, fns_cluster):
        dict_contig_cl = load_cluster_map(fn_c)
        contig_names = load_contig_names(fn_n)
        # Biuld df
        skipped = []
        dfg = []
        for contig in contig_names:
            cl = dict_contig_cl.get(contig)
            entry = dict_contig_entry.get(contig)
            if (cl is not None) & (entry is not None):
                dfg.append([contig, cl, entry, gene])
            else:
                skipped.append(contig)
        if len(dfg) == 0:
            raise ValueError(
                f"No contigs from hmmsearch in clustering \n first ten: {skipped[:10]}"
            )
        dfg = pd.DataFrame(dfg, columns=['contig','cluster','entry_id','gene'])
        df.append(dfg)
    df_merge = pd.concat(df, ignore_index=True)

    # merge
    _, dict_entry_tid, dict_entry_dtp = load_metadata(fn_meta)
    df_merge["taxid"] = df_merge["entry_id"].map(dict_entry_tid)
    df_merge["data_type"] = df_merge["entry_id"].map(dict_entry_dtp)
    df_merge["source"] = source

    # write
    df_merge.to_csv(fn_out, index=False)

    return


def main():
    # Merge marferret table
    fn_entry_id_map = "/Users/benjamingrodner/work/armbrust/data/prevalence/marferret/marfmmdb_contig_source-fbp1_isips_fre_aldo_fldA.grep.tsv"

    genes = [
        "FBP1",
        "FRE",
        "ISIP1",
        "ISIP2",
        "ISIP3",
        "fldA",
        "ALDO",
    ]

    fns_names = [
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/fbp1_hmmsearch_T30.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/fre_hmmsearch_T30.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/isip1_hmmsearch_T30.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/isip2_hmmsearch_T30.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/isip3_hmmsearch_T30.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/flda_hmmsearch_T30.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/hmmsearch_names/aldo_hmmsearch_T30.names",
    ]

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

    fn_out_marf = "/Users/benjamingrodner/work/armbrust/data/prevalence/marferret/prevalence_merge_table.csv"

    source = 'marferret'

    merge_write_table(
        fn_entry_id_map, genes, fns_names, fns_cluster, fn_meta, source, fn_out_marf
    )

    # shiri isolate table
    fn_entry_id_map = "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/entry_id_map_shiri.tsv"

    genes = [
        "FBP1",
        "FRE",
        "ISIP1",
        "ISIP2",
        "ISIP3",
        "fldA",
        "ALDO",
    ]

    fns_names = [
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/FBP1.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/FRE.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ISIP1.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ISIP2.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ISIP3.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/fldA.names",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ALDO.names",
    ]

    fns_cluster = [
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/FBP1_clusters.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/FRE_clusters.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ISIP1_clusters.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ISIP2_clusters.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ISIP3_clusters.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/fldA_clusters.tsv",
        "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch/ALDO_clusters.tsv",
    ]

    fn_meta = "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/metadata_custom_shiri.csv"

    fn_out_shiri = "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/prevalence_merge_table.csv"

    source = "shiri_isolates"

    merge_write_table(
        fn_entry_id_map, genes, fns_names, fns_cluster, fn_meta, source, fn_out_shiri
    )

    # combine tables
    df = []
    for fn in [fn_out_marf, fn_out_shiri]:
        df.append(pd.read_csv(fn))

    fn_out = "/Users/benjamingrodner/work/armbrust/data/prevalence/prevalence_merge_table.csv"
    pd.concat(df, axis=0, ignore_index=True).to_csv(fn_out, index=False)

    return

if __name__ == "__main__":

    main()
