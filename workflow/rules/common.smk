import glob
import re

def get_target_env_seqs(wildcards):
    files = []
    for gene in config["target_genes"]:
        files.extend(glob.glob(
            fmt_clustered_env_hitseqs.format(taxgene='*')
        ))
    return files


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
fn_db_clusters = dir_clust_db + '/' + bn_clust + '/db_clust_cluster.tsv'
sub = re.sub('0.','',str(config['subset_db_clusts']['pct']))
fn_db_rep_seqs_sub = dir_clust_db + '/' + bn_clust + f'/db_clust_rep_seq_sub{sub}pct.fasta'
bn_db_seqs = f'{bn_clust}_sub{sub}pct'

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

# Alignment
dir_aln = config['dir_out'] + f'/alignment/db_{bn_db_seqs}-env_{bn_env_clust}'
fmt_crystal_seqs = config['dir_rcsb'] + '/{rcsb_id}.fasta'
fn_db_crystal_seqs = dir_aln + '/db-env-crystal-manual.fasta'
fn_alignment = re.sub('.fasta','.aln',fn_db_crystal_seqs)
fn_trim_crystal = fn_alignment + '.trim_crystal'
fn_trim_crystal_startend = fn_alignment + '.trim_crystal_startend'
fn_trim_clip = fn_trim_crystal + '.clipkit'
frac_range = str(config['filter_alignment']['frac_thresh'])
fn_trim_clip_filt = fn_trim_clip + '.len_filt' + frac_range
fn_trim_clip_filt_dedup = fn_trim_clip_filt + '.dedup'
fn_trim_clip_filt_dedup_map = fn_trim_clip_filt + '.json'

# Build DB Tree
dir_tree = config['dir_out'] + f'/tree/db/{bn_db_seqs}'
bn_tree = 'db-crystal-manual'
fn_alignment_noenv = f'{dir_tree}/{bn_tree}.aln.trim_crystal.clipkit.len_filt{frac_range}'
fn_alignment_noenv_clip = f'{dir_tree}/{bn_tree}.aln.trim_crystal.clipkit.len_filt{frac_range}.clipkit'
fn_alignment_noenv_cliplog = f'{fn_alignment_noenv_clip}.log'
dir_fasttree = f'{dir_tree}/fasttree'
dir_ft_boot = f'{dir_fasttree}/bootstraps'
ext_bootstrap = '{rep}.fa'  # extension must be {rep}.fa 
fmt_bootstrap_resample = dir_ft_boot + f'/resamples/rep_{ext_bootstrap}'
bn_bootstrap_resample = re.sub(ext_bootstrap,'',fmt_bootstrap_resample)
fmt_bootstrap_fasttree = dir_ft_boot + '/trees/rep_{rep}.tree'
fn_bootstrap_fasttree_merged = dir_ft_boot + '/fasttree_boostraps_merged.tree'
fn_fasttree = f'{dir_fasttree}/fasttree_full.tree'
fn_headers_crystal_manual = f'{dir_tree}/crystal-manual.headers'
dir_roguenarok = f'{dir_tree}/roguenarok'
fn_roguenarok = f'{dir_roguenarok}/RogueNaRok_droppedRogues.{bn_tree}'  # Must be RogueNaRok_droppedRogues.<basename>
fn_rogues_to_drop = f'{dir_roguenarok}/rogues_to_drop.{bn_tree}'  # Must be RogueNaRok_droppedRogues.<basename>
fn_alignment_noenv_clip_drop = fn_alignment_noenv_clip + '.drop_rogues'
fn_full_tree_done = f'{dir_tree}/raxml_tree_done.{bn_tree}'


# Env tree placement
dir_env_tree = f'{dir_tree}/env_placement/{bn_env_clust}'
bn_env_target = ''
for gene in config['target_genes']:
    bn_env_target += gene + '_'
bn_env_target += 'env_clust'
bn_env_tree = f'{bn_env_target}-db-crystal-manual'
fn_env_aligned_mask = f'{dir_env_tree}/{bn_env_tree}_masked.aln'
fn_place_env_tree_done = (
    f'{dir_env_tree}/tree_done.{bn_env_tree}'
)
