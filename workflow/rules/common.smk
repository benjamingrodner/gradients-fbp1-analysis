import glob


def get_fns_seqs_to_search(wildcards):
    fns = [fmt_seqs_fasta.format(batch=b) for b in config['batches']]
    fns += glob.glob(config['dirs_data']['isolates'] + '/*')
    fns.append(config['path_database'])
    return fns


GENES = list(config['dict_gene_hmmprofile'].keys())

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

# Build Tree
dir_tree = config['dir_out'] + '/tree'
fn_db_seqs_to_cluster = dir_tree + '/db_cluster/hmmhitseqs.fasta'
ident = re.sub('0.','',str(config['cluster_db_seqs']['min_seq_id']))
cov = re.sub('0.','',str(config['cluster_db_seqs']['coverage']))
mode = config['cluster_db_seqs']['cov_mode']
bn_clust = dir_tree + f'/db_cluster/mmseqs2_i{ident}_c{cov}_mode{mode}'
fn_db_rep_seqs = bn_clust + '_rep_seq.fasta'
fn_db_clusters = bn_clust + '_cluster.tsv'
fmt_crystal_seqs = config['dir_rcsb'] + '/{rcsb_id}.fasta'
fn_db_crystal_seqs = dir_tree + '/alignment/db_and_crystal_struct_seqs.fasta'
fn_db_crystal_seqs = dir_tree + '/alignment/db_and_crystal_struct_seqs.fasta'
fn_alignment = dir_tree + '/alignment/db_and_crystal_struct_seqs.aln'
fn_trim_crystal = fn_alignment + '.trim_crystal'
fn_trim_crystal_startend = fn_alignment + '.trim_crystal_startend'
fn_trim_clip = fn_trim_crystal + '.clipkit'
frac_range = str(config['filter_alignment']['frac_thresh'])
fn_trim_clip_filt = fn_trim_clip + '.frac_range' + frac_range
fn_ml_tree = dir_tree + '/ML/RAxML_bestTree.db_and_crystal_struct_seqs.tree'
fn_bootstraps = dir_tree + '/bootstrap/RAxML_bootstrap.db_and_crystal_struct_seqs.bootstrap'
fn_ml_bootstrap = dir_tree + '/ML-bootstrap/RAxML_bipartitions.db_and_crystal_struct_seqs.ml_bootstrap'
fn_mrc = dir_tree + '/MRC/RAxML_result.db_and_crystal_struct_seqs.mrc'