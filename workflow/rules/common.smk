import glob


def get_fns_seqs_to_search(wildcards):
    fns = [fmt_seqs_fasta.format(batch=b) for b in config['batches']]
    fns += glob.glob(config['dirs_data']['isolates'] + '/*')
    fns.append(config['path_database'])
    return fns

def aggregate_env_clusters(wildcards):
    chunk_ids = []
    for batch in BATCHES:
        # 1. Access the checkpoint output
        checkpoint_output = checkpoints.group_env_hitnames.get(
            batch=batch
        ).output[0]
        # 2. Get filenames and strip extensions/prefixes
        # Suppose files are: chunk_aa.txt, chunk_ab.txt
        # We want chunk_id to be just: aa, ab
        for f in os.listdir(checkpoint_output):
            if f.startswith(tg_prefix):
                # Remove 'chunk_' prefix and '.txt' suffix
                cid = f.replace(tg_prefix, "").replace(ext_env_hitnames, "")
                chunk_ids.append(cid)
    # Get unique since in checkpoint we get taxgene repeats from multiple batches, but later we merged batches
    chunk_ids = list(set(chunk_ids))
    # 3. Reconstruct the paths for the next rule
    out = [
        fmt_clustered_env_hitseqs.format(taxgene=cid) 
        for cid in chunk_ids
    ]
    return out


GENES = list(config['dict_gene_hmmprofile'].keys())
BATCHES = list(config['batches'])

fmt_seqs_parquet = (
    config['dirs_data']['metatranscriptomes'] 
    + "/seqs-{batch}.parquet"
)
# Hmmsearch
dir_hmmsearch = config['dir_out'] + "/hmmsearch"
fmt_seqs_fasta = dir_hmmsearch + "/fastas_to_search/{batch}.fasta"
fn_seqs_to_search = dir_hmmsearch + "/fastas_to_search/to_search.fasta"
bn_hmmsearch = (
    dir_hmmsearch 
    + "/{gene}_hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
)
fmt_table_hmm = bn_hmmsearch + ".tbl"
fmt_table_hmm_domain = bn_hmmsearch + ".domain"
fmt_stdout_hmm = bn_hmmsearch + ".stdout"
fmt_hmm_hitnames = bn_hmmsearch + ".names"
fn_hmm_hitnames_all = (
    dir_hmmsearch + "/hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
    + '_all.names'
)
fn_hmms_best_hit = dir_hmmsearch + '/hmms_best_hit.tsv'


# Cluster db
dir_clust_db = config['dir_out'] + '/cluster/db'
fn_db_seqs_to_cluster = dir_clust_db + '/hmmhitseqs.fasta'
ident = re.sub('0.','',str(config['cluster_db_seqs']['min_seq_id']))
cov = re.sub('0.','',str(config['cluster_db_seqs']['coverage']))
mode = config['cluster_db_seqs']['cov_mode']
bn_clust = f'mmseqs2_i{ident}_c{cov}_mode{mode}'
fn_db_rep_seqs = dir_clust_db + '/' + bn_clust + '/db_clust_rep_seq.fasta'
sub = re.sub('0.','',str(config['subset_db_clusts']['pct']))
fn_db_rep_seqs_sub = dir_clust_db + '/' + bn_clust + f'/db_clust_rep_seq_sub{sub}pct.fasta'
fn_db_clusters = dir_clust_db + '/' + bn_clust + '/db_clust_cluster.tsv'
bn_db_seqs = f'{bn_clust}_sub{sub}pct'

# Alignment db
dir_aln = config['dir_out'] + f'/alignment/db/{bn_db_seqs}'
fmt_crystal_seqs = config['dir_rcsb'] + '/{rcsb_id}.fasta'
fn_db_crystal_seqs = dir_aln + '/db_and_crystal_struct_seqs.fasta'
fn_alignment = dir_aln + '/db_and_crystal_struct_seqs.aln'
fn_trim_crystal = fn_alignment + '.trim_crystal'
fn_trim_crystal_startend = fn_alignment + '.trim_crystal_startend'
fn_trim_clip = fn_trim_crystal + '.clipkit'
frac_range = str(config['filter_alignment']['frac_thresh'])
fn_trim_clip_filt = fn_trim_clip + '.frac_range' + frac_range

# Build DB Tree
dir_tree = config['dir_out'] + f'/tree/db/{bn_db_seqs}'
fn_ml_tree = dir_tree + '/ML/RAxML_bestTree.db_and_crystal_struct_seqs.tree'
fn_bootstraps_done = dir_tree + '/bootstrap/RAxML_bootstrap.db_and_crystal_struct_seqs.done'
fn_ml_bootstrap = dir_tree + '/ML-bootstrap/RAxML_bipartitions.db_and_crystal_struct_seqs.ml_bootstrap'
fn_mrc = dir_tree + '/MRC/RAxML_result.db_and_crystal_struct_seqs.mrc'
fn_full_tree_done = dir_tree + '/full_analysis/db_and_crystal_struct_seqs.done'

# Cluster env
dir_env_clust = config['dir_out'] + '/cluster/env'
fmt_bigtable = (
    config['dirs_data']['metatranscriptomes'] 
    + "/merge_counts_tax_gene-{batch}.parquet"
)
tg_prefix = 'taxgene_'
ext_env_hitnames = '.txt'
fmt_env_hitnames_dir_tmp = dir_env_clust + '/{batch}/hitnames_tmp'
fmt_env_hitnames_tmp = (
    dir_env_clust + '/{batch}/hitnames_tmp/' 
    + tg_prefix + '{taxgene}' + ext_env_hitnames
)
fmt_env_hitnames = (
    dir_env_clust + '/{batch}/hitnames/' 
    + tg_prefix + '{taxgene}' + ext_env_hitnames
)
# fmt_env_hitseqs_dir = (
#     dir_env_clust + '/hitseqs/{batch}'
# )
fmt_env_hitseqs = (
    dir_env_clust + '/{batch}/hitseqs/' + tg_prefix + '{taxgene}.faa'
)
fmt_grouped_env_hitseqs = (
    dir_env_clust + '/hitseqs_batch_grouped/' 
    + tg_prefix + '{taxgene}.faa'
)
ident = re.sub('0.','',str(config['cluster_env_hitseqs']['min_seq_id']))
cov = re.sub('0.','',str(config['cluster_env_hitseqs']['coverage']))
mode = config['cluster_env_hitseqs']['cov_mode']
bn_env_clust = f'mmseqs2_i{ident}_c{cov}_mode{mode}'
fmt_clustered_env_hitseqs = (
    dir_env_clust + '/' + bn_env_clust + '/taxgene_groups/{taxgene}/' 
    + tg_prefix + '{taxgene}_rep_seq.fasta'
)
fn_env_seqs_clust_cat = (
    dir_env_clust + '/' + bn_env_clust + '/env_seqs_clust_cat.faa'
)