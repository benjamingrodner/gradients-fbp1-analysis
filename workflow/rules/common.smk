from pathlib import Path
import glob
import yaml
import re
import os

def get_dirs_exp_deseq_plot(wildcards):
    drs = []
    for exp in DICT_EXP.keys():
        d = checkpoints.group_exp_hitnames.get(
            exp=exp
        ).output[0]
        fns = glob.glob(f'{d}/*')
        if len(fns) > 0:
            drs.append(dir_exp_deseq_plot.format(exp=exp))
    return drs

def get_fn_exp_counts_merge(exp):
    m = DICT_EXP[exp]['method_counts']
    if 'salmon' in m:
        return fmt_salmon_counts_merge.format(exp_quant=exp)
    elif m in ['from_author','download']:
        return fmt_fromauthor_and_downloaded_counts_merge.format(exp_auth=exp)
    else:
        raise ValueError(f"method_counts is not defined for experiment {exp}")


def get_fn_sample_info(exp):
    m = DICT_EXP[exp]['method_counts']
    if m == 'salmon-bioproject':
        return fmt_bioproject_info.format(exp_bioproject=exp)
    elif m == 'salmon-biosamples':
        return fmt_biosample_info.format(exp_biosample=exp)
    elif m == 'salmon-srr':
        return fmt_srr_info.format(exp_srr=exp)
    else:
        return 'None'

def get_colname_contigs(exp):
    fn = DICT_EXP[exp].get('colname_contigs')
    if fn is None:
        fn = "contig_name"
    return fn

def get_script_merge_counts_auth(exp):
    fn = DICT_EXP[exp].get('script_merge_counts')
    if fn is None:
        fn = "workflow/scripts/cp_counts_author_to_parquet.py"
    return fn

def get_colname_merge_counts_auth(exp):
    fn = DICT_EXP[exp].get('script_merge_counts')
    colname = ""
    if fn is None:
        colname = DICT_EXP[exp]['colname_contigs']
    return colname


def get_glob_counts_auth(exp):
    fn = DICT_EXP[exp].get('glob_counts')
    if fn is None:
        d = dir_download_counts.format(exp_download_counts=exp)
        fn = f'{d}/*'
    return fn

def expand_exp_quant(fmt):
    fns = []
    for exp, dinfo in DICT_EXP.items():
        m = dinfo['method_counts']
        if 'salmon' in m:
            fns.append(fmt.format(exp_quant=exp))
    return fns
def expand_exp_auth(fmt):
    fns = []
    for exp, dinfo in DICT_EXP.items():
        m = dinfo['method_counts']
        if m in ['from_author','download']:
            fns.append(fmt.format(exp_auth=exp))
    return fns

def aggregate_exp_salmon_counts(wildcards):
    return expand_biosamples(wildcards, fmt_quant, get_reads=False)

def aggregate_fastqc_trimmed(wildcards):
    return expand_biosamples(wildcards, fmt_fastqc_trimmed)

def aggregate_fastqc_raw(wildcards):
    return expand_biosamples(wildcards, fmt_fastqc_raw)

def expand_biosamples(wildcards, fmt, get_reads=True):
    exp = wildcards.exp_quant
    # _ = checkpoints.gzip_biosample_fastqs.get(**wildcards)
    d = checkpoints.merge_biosample_fastqs.get(exp_quant=exp).output[0]
    fns_d = glob.glob(f'{d}/*')
    regex = fmt_biosample_fastq_merged.format(
        exp_quant=exp, sample=r"(?P<sample>\w+)", read=r"(?P<read>\w+)"
    )
    fns = []
    for fn in fns_d:
        match = re.search(regex, fn)
        try:
            sample = match.group('sample')
        except:
            raise ValueError(f"Could not match sample with regex\n{regex}\non file\n{fn}")
        read = ""
        if get_reads:
            try:
                read = match.group('read')
            except:
                raise ValueError(f"Could not match read with regex\n{regex}\non file\n{fn}")
        fn_ = fmt.format(exp_quant=exp, sample=sample, read=read)     
        fns.append(fn_)
    return list(set(fns))


def get_adapter_file(exp):
    fn = DICT_EXP[exp].get('adapter_file')
    if fn is None:
        fn = "None"
    return fn

def get_trim_opts(exp):
    try:
        fn = DICT_EXP[exp].get('trimmomatic_opts')
    except: 
        raise ValueError(f"trimmomatic_opts is not defined for experiment {exp}")
    return fn

def check_merge_biosample(exp):
    m = DICT_EXP[exp]['method_counts']
    merge = "yes"
    if m == 'salmon-srr':
        merge = "no"
    return merge

def get_srr_info(exp):
    m = DICT_EXP[exp]['method_counts']
    if m == 'salmon-biosamples':
        return fmt_biosample_srr_info.format(exp_biosample=exp)
    elif m == 'salmon-bioproject':
        return fmt_bioproject_srr_info.format(exp_bioproject=exp)
    elif m == 'salmon-srr':
        return fmt_srr_info.format(exp_srr=exp)
    else:
        raise ValueError(f"Experiment {exp} is not defined for quantification with method {m}")
        

def get_dir_fastq(exp):
    m = DICT_EXP[exp]['method_counts']
    if m == 'salmon-biosamples':
        return dir_biosample_fastq.format(exp_biosample=exp)
    elif m == 'salmon-bioproject':
        return dir_bioproject_fastq.format(exp_bioproject=exp)
    elif m == 'salmon-srr':
        return dir_srr_fastq.format(exp_srr=exp)
    else:
        raise ValueError(f"Experiment {exp} is not defined for salmon quant with method {m}")
    

def get_fns_downloads_done(wildcards):
    fns = []
    for exp, dinfo in DICT_EXP.items():
        m = dinfo['method_counts']
        if m == 'salmon-biosamples':
            fn = fmt_biosample_fastq_done.format(exp_biosample=exp)
        elif m == 'salmon-bioproject':
            fn = fmt_bioproject_fastq_done.format(exp_bioproject=exp)
        elif m == 'salmon-srr':
            fn = fmt_srr_fastq_done.format(exp_srr=exp)
        elif m == 'download':
            fn = fmt_download_counts_done.format(exp_download_counts=exp)
        elif 'from_author' in m:
            continue
        else:
            raise ValueError(f"Method {m} is not available for getting the experiment counts")
        fns.append(fn)
    return fns


def get_elink_target(exp):
    elink = DICT_EXP[exp].get('elink_target')
    if elink is None:
        elink = 'biosample'
    return elink

def get_exp_assm_fn_or_link(exp):
    method = DICT_EXP[exp]['method_assembly']
    if method == 'download':
        return DICT_EXP[exp]['link_assembly']
    elif method == 'from_author':
        return DICT_EXP[exp]['fn_assembly']
    else:
        raise ValueError(f"Method {method} is not available for getting the experiment assembly")

def get_fn_exp_seqs_clust_target(wildcards):
    fn_sub = config['fn_manual_subset_clusters']
    fn = fn_exp_seqs_clust_target
    if fn_sub is not None:
        fn = fn_exp_seqs_clust_target_sub
    return fn



# def get_exp_rep_seqs(wildcards):
#     fns = []
#     for exp in DICT_EXP.keys():
#         d = checkpoints.group_exp_hitnames.get(
#             exp=exp
#         ).output[0]
#         for f in os.listdir(d):
#             if f.startswith(best_hit_prefix):
#                 gene = f.replace(best_hit_prefix, "").replace(best_hit_ext,"")
#                 fn = fmt_exp_rep_seqs.format(exp=exp, gene=gene)
#                 fns.append(fn)
#     return fns
def aggregate_exp_rep_seqs_target(wildcards):
    fns = []
    for exp in DICT_EXP.keys():
        d = checkpoints.group_exp_hitnames.get(
            exp=exp
        ).output[0]
        for f in os.listdir(d):
            if f.startswith(best_hit_prefix):
                gene = f.replace(best_hit_prefix, "").replace(best_hit_ext,"")
                if gene in TGENES:
                    fn = fmt_exp_rep_seqs.format(exp=exp, gene=gene)
                    fns.append(fn)
    return fns

def get_exp_clusters(exp):
    fns = []
    d = checkpoints.group_exp_hitnames.get(
        exp=exp
    ).output[0]
    for f in os.listdir(d):
        if f.startswith(best_hit_prefix):
            gene = f.replace(best_hit_prefix, "").replace(best_hit_ext,"")
            fn = fmt_exp_clusters.format(exp=exp, gene=gene)
            fns.append(fn)
    return fns

def get_fn_env_seqs_clust_target(wildcards):
    fn_sub = config['fn_manual_subset_clusters']
    fn = fn_env_seqs_clust_target
    if fn_sub is not None:
        fn = fn_env_seqs_clust_target_sub
    return fn


def get_env_rep_seqs(wildcards):
    return glob.glob(
            fmt_clustered_env_hitseqs.format(taxgene='*')
        )

def get_target_exp_rep_seqs(wildcards):
    fns = []
    for exp in DICT_EXP.keys():
        d = checkpoints.group_exp_hitnames.get(
            exp=exp
        ).output[0]
        for f in os.listdir(d):
            if f.startswith(best_hit_prefix):
                gene = f.replace(best_hit_prefix, "").replace(best_hit_ext,"")
                if gene in config['target_genes_for_exp_placement']:
                    fn = fmt_exp_rep_seqs.format(exp=exp, gene=gene)
                    fns.append(fn)
    return fns

def get_target_env_rep_seqs(wildcards):
    fns = glob.glob(
        fmt_clustered_env_hitseqs.format(taxgene='*')
    )
    out = []
    for fn in fns:
        if any([g in fn for g in config['target_genes_for_env_placement']]):
            out.append(fn)
    return out

def get_regex_seqname(exp):
    r = DICT_EXP[exp].get('regex_replace_seqname')
    if r is None:
        r = ['^','']
    return r

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

def aggregate_env_clusters_target(wildcards):
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
    out = []
    for cid in chunk_ids:
        if any([gene in cid for gene in TGENES]):
            out.append(fmt_clustered_env_hitseqs.format(taxgene=cid))
    return out

def get_fns_hmm_hitnames_nontarget(wildcards):
    fns = []
    for gene in OGENES:
        fns.append(fmt_hmm_hitnames.format(gene=gene))
    return fns

def get_fns_hmm_hitnames_target(wildcards):
    fns = []
    for gene in TGENES:
        fns.append(fmt_hmm_hitnames.format(gene=gene))
    return fns

# Global 
GENES = list(config['dict_gene_hmmprofile'].keys())
TGENES = list(set(
    config['target_genes_for_exp_placement'] 
    + config['target_genes_for_env_placement']
))
OGENES = config['outgroup_genes']
BATCHES = list(config['batches'])
with open(config['fn_experiment_info'], 'r') as f:
    DICT_EXP = yaml.safe_load(f)

READS=['1','2']
# Data format
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
bn_hmm_genes = ''
for gene in GENES:
    bn_hmm_genes += f'{gene}_'
bn_hmm_genes = bn_hmm_genes.rstrip('_')

bn_hmm_target = ''
for gene in TGENES:
    bn_hmm_target += f'{gene}_'
bn_hmm_target = bn_hmm_target.rstrip('_')

bn_hmm_outgroup = ''
for gene in OGENES:
    bn_hmm_outgroup += f'{gene}_'
bn_hmm_outgroup = bn_hmm_outgroup.rstrip('_')

fn_hmm_hitnames_all = (
    dir_hmmsearch + f"/{bn_hmm_genes}_hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
    + '_all.names'
)
fn_hmms_best_hit = dir_hmmsearch + f'/{bn_hmm_genes}_hmms_best_hit.tsv'
fn_hmm_hitnames_target = (
    dir_hmmsearch + f"/{bn_hmm_target}_hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
    + '_all.names'
)
fn_hmm_hitnames_outgroup = (
    dir_hmmsearch + f"/{bn_hmm_outgroup}_hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
    + '_all.names'
)
fn_outgroup_db_seqs = (
    dir_hmmsearch + f"/{bn_hmm_outgroup}_hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
    + '_db.faa'
)
fn_outgroup_db_seqs_sub = (
    dir_hmmsearch + f"/{bn_hmm_outgroup}_hmmsearch_T" + str(config['hmmsearch']['thresh_score'])
    + f"_db_subn{config['subset_db_outgroups']['number']}.faa"
)

# Cluster db
dir_clust_db = config['dir_out'] + '/cluster/db'
fn_db_seqs_to_cluster = f'{dir_clust_db}/{bn_hmm_target}_hmmhitseqs.fasta'
ident = re.sub('0.','',str(config['cluster_db_seqs']['min_seq_id']))
cov = re.sub('0.','',str(config['cluster_db_seqs']['coverage']))
mode = config['cluster_db_seqs']['cov_mode']
bn_clust = f'mmseqs2_i{ident}_c{cov}_mode{mode}'
dir_clust_db_genes = f'{dir_clust_db}/{bn_hmm_target}_{bn_clust}'
fn_db_rep_seqs = f'{dir_clust_db_genes}/db_clust_rep_seq.fasta'
fn_db_clusters = f'{dir_clust_db_genes}/db_clust_cluster.tsv'
sub = config['subset_db_clusts']['pct']
fn_db_rep_seqs_sub = f'{dir_clust_db_genes}/db_clust_rep_seq_sub{sub}pct.fasta'
bn_clust_sub = f'{bn_clust}_sub{sub}pct'

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
fn_env_seqs_clust_target = (
    f'{dir_env_clust}/{bn_env_clust}/{bn_hmm_target}_env_seqs_clust_cat.faa'
)
fn_env_seqs_clust_target_sub = (
    f'{dir_env_clust}/{bn_env_clust}/{bn_hmm_target}_env_seqs_clust_cat-manual_sub.faa'
)

# Experiment downloads
dir_exp_data = config['dirs_data']['experiments']
fmt_exp_assembly = dir_exp_data + '/{exp}/assembly_og.fasta'
fmt_exp_assembly_rename = dir_exp_data + '/{exp}/assembly_rename.fasta'
fmt_exp_assembly_6tr = dir_exp_data + '/{exp}/assembly.faa'
# fn_download_assemblies_done = dir_exp_data + '/fns_downloaded.txt'

# Hmmsearch experiments
dir_exp = config['dir_out'] + '/experiments'
dir_exp_hmm = dir_exp + '/{exp}/hmmsearch'
bn_hmm_exp_gene = (
    dir_exp_hmm + '/{gene}/hmmsearch_T' 
    + str(config['hmmsearch']['thresh_score'])
)
fmt_table_hmm_exp = bn_hmm_exp_gene + '.tbl'
fmt_table_hmm_domain_exp = bn_hmm_exp_gene + '.domtab'
fmt_stdout_hmm_exp = bn_hmm_exp_gene + '.out'
fmt_hmm_hitnames_exp = bn_hmm_exp_gene + '.hitnames'
fn_hmm_hitnames_exp_all = dir_exp_hmm + '/headers_merged.txt'
fmt_hmms_best_hit_exp = dir_exp_hmm + '/best_hit_genes.tsv'
dir_exp_hitnames_grouped_gene = dir_exp_hmm + '/best_hit_genes_grouped'
best_hit_prefix = 'hitnames_best_hit_'
best_hit_ext = '.txt'
fmt_exp_hitnames_grouped_gene = (
    dir_exp_hitnames_grouped_gene + '/' + best_hit_prefix + '{gene}' + best_hit_ext
)

# Cluster experiments
ident = re.sub('0.','',str(config['cluster_exp_hitseqs']['min_seq_id']))
cov = re.sub('0.','',str(config['cluster_exp_hitseqs']['coverage']))
mode = config['cluster_exp_hitseqs']['cov_mode']
bn_exp_clust = f'mmseqs2_i{ident}_c{cov}_mode{mode}'
dir_exp_clust = dir_exp + '/{exp}/cluster/' + bn_exp_clust + '/{gene}'
fmt_exp_seqs_to_cluster = f'{dir_exp_clust}/seqs_to_cluster.fasta'
fmt_exp_rep_seqs = f'{dir_exp_clust}/exp_clust_rep_seq.fasta'
fmt_exp_clusters = f'{dir_exp_clust}/exp_clust_cluster.tsv'
fn_exp_seqs_clust_target = (
    f'{dir_exp}/merge_clusters/{bn_exp_clust}'
    + f'/{bn_hmm_target}_all_exp_clust_rep_seqs.fasta'
)
fn_exp_seqs_clust_target_sub = (
    f'{dir_exp}/merge_clusters/{bn_exp_clust}'
    + f'/{bn_hmm_target}_all_exp_clust_rep_seqs-manual_sub.fasta'
)


# Alignment
dir_aln = (
    config['dir_out'] 
    + f'/alignment/target_{bn_hmm_target}-db_{bn_clust_sub}-env_{bn_env_clust}-outgroup_{bn_hmm_outgroup}'
    )
fmt_crystal_seqs = config['dir_rcsb'] + '/{rcsb_id}.fasta'
fn_db_crystal_seqs = dir_aln + '/db-env-exp-crystal-manual.fasta'
fn_db_crystal_seqs_dedup = dir_aln + '/db-env-exp-crystal-manual.dedup.fasta'
fn_db_crystal_seqs_dedup_removed = dir_aln + '/db-env-exp-crystal-manual.removed.fasta'
fn_alignment = re.sub('.fasta','.aln',fn_db_crystal_seqs)
fn_trim_crystal = fn_alignment + '.trim_crystal'
fn_trim_crystal_startend = fn_alignment + '.trim_crystal_startend'
fn_trim_clip = fn_trim_crystal + '.clipkit'
frac_range = str(config['filter_alignment']['frac_thresh'])
fn_trim_clip_filt = fn_trim_clip + '.len_filt' + frac_range
fn_trim_clip_filt_dedup = fn_trim_clip_filt + '.dedup'
fn_trim_clip_filt_dedup_map = fn_trim_clip_filt + '.json'
fn_trim_clip_dedup = fn_trim_clip + '.dedup'
fn_trim_clip_dedup_map = fn_trim_clip_dedup + '.json'

# Build DB Tree
bn_backbone = f'backbone_{bn_hmm_genes}_{bn_clust_sub}'
dir_tree = config['dir_out'] + f'/tree/backbone/{bn_backbone}'
bn_tree = 'db-crystal-manual'
fn_alignment_noenv = f'{dir_tree}/{bn_tree}.aln.trim_crystal.clipkit'
fn_alignment_noenv_clip = f'{fn_alignment_noenv}.clipkit'
fn_alignment_noenv_cliplog = f'{fn_alignment_noenv_clip}.log'
fn_alignment_noenv_clip_filt = f'{fn_alignment_noenv_clip}.len_filt{frac_range}'
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
fn_alignment_noenv_clip_filt_drop = fn_alignment_noenv_clip_filt + '.drop_rogues'
dir_raxml = f'{dir_tree}/raxml'
fn_full_tree_done = f'{dir_raxml}/raxml_tree_done.{bn_tree}'
fn_tree_support = f'{dir_raxml}/{bn_tree}.raxml.support'
nml = re.sub(r'[\{\},]','',config['build_tree']['n_trees_extra'])
nboot = config['build_tree']['n_bootstraps_extra']
seed = config['build_tree']['seed_extra']
bn_extra_ml_and_bootstraps = f'extra_ml_{nml}_and_bootstraps_{nboot}_seed_{seed}'
fn_extra_ml_and_bootstraps_done = f'{dir_tree}/extra_raxml_done.txt'
fn_extra_tree_support = f'{dir_raxml}/{bn_extra_ml_and_bootstraps}.raxml.supportTBE'

# Env tree placement
bn_env_target = ''
for gene in config['target_genes_for_env_placement']:
    bn_env_target += gene + '_'
bn_env_target = bn_env_target.rstrip('_')
dir_env_tree = config['dir_out'] + f'/tree/env_placement/{bn_backbone}/{bn_env_clust}/{bn_env_target}'
bn_env_tree = f'env_{bn_env_target}-db_{bn_hmm_genes}-crystal-manual'
fn_env_aligned_mask = f'{dir_env_tree}/{bn_env_tree}_masked.aln'
fn_env_aligned_mask_dedup = f'{fn_env_aligned_mask}.dedup'
fn_env_aligned_mask_dedup_map = f'{fn_env_aligned_mask_dedup}.map'
fn_env_aligned_mask_dedup_tfilt = f'{fn_env_aligned_mask_dedup}.target_gene_filt'
fn_env_aligned_mask_dedup_tfilt_filt = f'{fn_env_aligned_mask_dedup_tfilt}.short_long_filt'
dir_env_tree_raxml = f'{dir_env_tree}/raxml'
fn_place_env_tree_done = (
    f'{dir_env_tree}/tree_done.{bn_env_tree}'
)
fn_place_env_tree = (
    f'{dir_env_tree_raxml}/RAxML_labelledTree.{bn_env_tree}'
)

# Experiment tree placement
bn_exp_target = ''
for gene in config['target_genes_for_exp_placement']:
    bn_exp_target += gene + '_'
bn_exp_target = bn_exp_target.rstrip('_')
dir_exp_tree = config['dir_out'] + f'/tree/exp_placement/{bn_backbone}/{bn_exp_clust}/{bn_exp_target}'
bn_exp_tree = f'exp_{bn_exp_target}-db_{bn_hmm_genes}-crystal-manual'
fn_exp_aligned_mask_dedup_tfilt = f'{dir_exp_tree}/{bn_exp_tree}_masked.dedup.target_gene_filt'
fn_exp_aligned_mask_dedup_tfilt_filt = f'{fn_exp_aligned_mask_dedup_tfilt}.short_long_filt'
dir_exp_tree_raxml = f'{dir_exp_tree}/raxml'
fn_place_exp_tree_done = (
    f'{dir_exp_tree}/tree_done.{bn_exp_tree}'
)
fn_place_exp_tree = (
    f'{dir_exp_tree_raxml}/RAxML_labelledTree.{bn_exp_tree}'
)

# Target placement for env + exp
bn_target_pl = ''
for gene in TGENES:
    bn_target_pl += gene + '_'
bn_target_pl = bn_target_pl.rstrip('_')
dir_target_pl_tree = (
    config['dir_out'] 
    + f'/tree/target_placement/{bn_backbone}/exp_{bn_exp_clust}-env_{bn_env_clust}/{bn_target_pl}'
)
bn_target_pl_tree = f'expenv_{bn_target_pl}-db_{bn_hmm_genes}-crystal-manual'
fn_target_aligned_mask_dedup_tfilt = f'{dir_target_pl_tree}/{bn_target_pl_tree}_masked.dedup.target_gene_filt'
fn_target_aligned_mask_dedup_tfilt_filt = f'{fn_target_aligned_mask_dedup_tfilt}.short_long_filt'
dir_target_pl_tree_raxml = f'{dir_target_pl_tree}/raxml'
fn_place_target_tree_done = (
    f'{dir_target_pl_tree}/tree_done.{bn_target_pl_tree}'
)
fn_place_target_tree = (
    f'{dir_target_pl_tree_raxml}/RAxML_labelledTree.{bn_target_pl_tree}'
)
fn_place_target_jplace = (
    f'{dir_target_pl_tree_raxml}/RAxML_portableTree.{bn_target_pl_tree}.jplace'
)
fn_place_target_jplace_support = (
    f'{dir_target_pl_tree_raxml}/RAxML_portableTree.{bn_target_pl_tree}.support.jplace'
)

# TODO: clear out old separate placement of experiment and environmental seqs
# Tree annotation
dir_annot = f'{dir_env_tree}/annotation'
# fn_table_annotate_env = f'{dir_annot}/annotations_env.csv'
# fn_table_annotate_db = f'{dir_annot}/annotations_db.csv'
# fn_table_annotate_crystal = f'{dir_annot}/annotations_crystal.csv'
# fn_table_annotate_manual = f'{dir_annot}/annotations_manual.csv'
# fn_table_annotations = f'{dir_annot}/annotations_merge.csv'
# fn_taxon_colorstrip = f'{dir_annot}/Taxon_colorstrip.txt'
# fn_domain_colorstrip = f'{dir_annot}/Domain_colorstrip.txt'
# fn_substrate_treecolors = f'{dir_annot}/Substrate_treecolors.txt'
# fn_source_treecolors = f'{dir_annot}/Source_treecolors.txt'
# fn_fbp1_colorstrip = f'{dir_annot}/Gene_colorstrip.txt'
# fn_crystal_symbol = f'{dir_annot}/crystal_symbol.txt'

# Experiment tree annotation
dir_annot_exp = f'{dir_exp_tree}/annotation'
# fn_table_annotate_exp = f'{dir_annot_exp}/annotations_exp.csv' 
# fn_taxon_colorstrip_exp = f'{dir_annot_exp}/Taxon_colorstrip.txt'
# fn_domain_colorstrip_exp = f'{dir_annot_exp}/Domain_colorstrip.txt'
# fn_substrate_treecolors_exp = f'{dir_annot_exp}/Substrate_treecolors.txt'
# fn_source_treecolors_exp = f'{dir_annot_exp}/Source_treecolors.txt'
# fn_fbp1_colorstrip_exp = f'{dir_annot_exp}/Gene_colorstrip.txt'

# exp + env annotations
dir_annot_target_pl = f'{dir_target_pl_tree}/annotation'
dir_annot_tbl = f'{dir_annot_target_pl}/tables'
fn_table_annotate_env = f'{dir_annot_tbl}/annotations_env.csv'
fn_table_annotate_db = f'{dir_annot_tbl}/annotations_db.csv'
fn_table_annotate_crystal = f'{dir_annot_tbl}/annotations_crystal.csv'
fn_table_annotate_manual = f'{dir_annot_tbl}/annotations_manual.csv'
fn_table_annotate_exp = f'{dir_annot_tbl}/annotations_exp.csv' 
fn_table_annotations = f'{dir_annot_tbl}/annotations_merge.csv'
dir_annot_itol = f'{dir_annot_target_pl}/itol'
fn_target_pl_annot_done = f'{dir_annot_target_pl}/annotation_done.txt'
dir_target_pl_tree_gappa = f'{dir_target_pl_tree}/gappa'

# experiment downloads
dir_download_counts = dir_exp_data + '/{exp_download_counts}/counts_download'
fmt_download_counts_done = dir_exp_data + '/{exp_download_counts}/download_counts_done.txt'

dir_exp_data_bpj = dir_exp_data + '/{exp_bioproject}/bioproject'
fmt_bioproject_info = f'{dir_exp_data_bpj}/biosample_info.txt'
fmt_bioproject_srr_info = f'{dir_exp_data_bpj}/srr_info.txt'
fmt_bioproject_srr_list = f'{dir_exp_data_bpj}/srr_list.txt'
dir_bioproject_fastq = f'{dir_exp_data_bpj}/reads'
fmt_bioproject_fastq_done = f'{dir_exp_data_bpj}/download_done.txt'

dir_exp_data_bs = dir_exp_data + '/{exp_biosample}/biosamples'
fmt_biosample_info = f'{dir_exp_data_bs}/biosample_info.txt'
fmt_biosample_srr_info = f'{dir_exp_data_bs}/srr_info.txt'
fmt_biosample_srr_list = f'{dir_exp_data_bs}/srr_list.txt'
dir_biosample_fastq = f'{dir_exp_data_bs}/reads'
fmt_biosample_fastq_done = f'{dir_exp_data_bs}/download_done.txt'

dir_exp_data_srr = dir_exp_data + '/{exp_srr}/srrs'
fmt_srr_info = f'{dir_exp_data_srr}/srr_info.txt'
fmt_srr_list = f'{dir_exp_data_srr}/srr_list.txt'
dir_srr_prefetch = f'{dir_exp_data_srr}/prefetch_reads'
dir_srr_fastq = f'{dir_exp_data_srr}/reads'
fmt_srr_fastq_done = f'{dir_exp_data_srr}/download_done.txt'

fn_read_and_count_downloads_done = f'{dir_exp_data}/download_reads_and_counts_done.txt'

# Experiment read prep
dir_read_prep = dir_exp + '/{exp_quant}/read_prep'
dir_biosamples_fastq_merged = f'{dir_read_prep}/merge_biosample'
fmt_biosample_fastq_merged_r1 = dir_biosamples_fastq_merged + '/{sample}_1.fastq.gz'
fmt_biosample_fastq_merged_r2 = dir_biosamples_fastq_merged + '/{sample}_2.fastq.gz'
fmt_biosample_fastq_merged = dir_biosamples_fastq_merged + '/{sample}_{read}.fastq.gz'
fmt_gzip_biosample_fastq_done = f'{dir_read_prep}/biosample_merge_and_gzip_done.txt'
# fn_agg_biosample_dirs = f'{dir_exp}/biosample_merge_done.txt'
dir_trimmed = f'{dir_read_prep}/trimmed'
fmt_fastq_trim = dir_trimmed + '/{sample}_{read}.fastq.gz'
fmt_fastq_trim_r1, fmt_fastq_trim_r2 = [
    fmt_fastq_trim.format(exp_quant='{exp_quant}', sample='{sample}', read=r) 
    for r in ['1','2']
]
fmt_fastq_trim_r1_unpaired, fmt_fastq_trim_r2_unpaired = [
    re.sub('.fastq.gz','_unpaired.fastq.gz', fn) 
    for fn in [fmt_fastq_trim_r1, fmt_fastq_trim_r2]
]
# (fmt_fastq_trim, fmt_fastq_trim_r1, fmt_fastq_trim_r1_unpaired, 
# fmt_fastq_trim_r2, fmt_fastq_trim_r2_unpaired) = [
#     fn + '.gz' for fn in [
#         fmt_fastq_trim_unz, fmt_fastq_trim_r1_unz, fmt_fastq_trim_r1_unpaired_unz,
#         fmt_fastq_trim_r2_unz, fmt_fastq_trim_r2_unpaired_unz
#     ]
# ]
dir_fastqc =  f'{dir_read_prep}/fastqc'
fmt_fastqc_raw = dir_fastqc + '/raw/{sample}_{read}_fastqc.html'
fmt_fastqc_trimmed = dir_fastqc + '/trimmed/{sample}_{read}_fastqc.html'
fmt_multiqc_raw = f'{dir_read_prep}/multiqc/raw/multiqc_report.html'
fmt_multiqc_trimmed = f'{dir_read_prep}/multiqc/trimmed/multiqc_report.html'
fn_read_prep_done = f'{dir_exp}/read_prep_done.txt'

# Experiment quant
fmt_exp_assembly_quant = dir_exp_data + '/{exp_quant}/assembly.fasta'
# fmt_exp_assembly_quant = fmt_exp_assembly_rename.format(exp='{exp_quant}')
dir_quant = dir_exp + '/{exp_quant}/sample_quant'
dir_salmon_idx = f'{dir_quant}/salmon_index'
fmt_quant = dir_quant + '/{sample}/quant.sf.gz'
ext_counts_agg = 'counts_agg.parquet'
fmt_salmon_counts_merge = f'{dir_quant}/salmon_{ext_counts_agg}'
fmt_fromauthor_and_downloaded_counts_merge = dir_exp + '/{exp_auth}/sample_quant/author_' + ext_counts_agg
fn_quant_done = f'{dir_exp}/salmon_and_author_quant_done.txt'

# deseq
fmt_exp_counts_merge = dir_exp + '/{exp}/sample_quant/' + ext_counts_agg
dir_exp_deseq = dir_exp + '/{exp}/deseq'
dir_exp_counts_clust = f'{dir_exp_deseq}/counts_clust-genes_{bn_hmm_target}-{bn_exp_clust}'
fmt_exp_counts_clust = f'{dir_exp_counts_clust}/counts.parquet'
fmt_exp_meta = f'{dir_exp_deseq}/metadata.csv'
dir_exp_deseq = f'{dir_exp_counts_clust}/stats'
dir_exp_deseq_plot = f'{dir_exp_counts_clust}/horizontal_barplots'
fn_deseq_done = f'{dir_exp}/deseq_done.txt'


