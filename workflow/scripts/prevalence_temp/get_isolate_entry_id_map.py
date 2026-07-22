import pandas as pd
from glob import glob
import os
import re
from collections import defaultdict

dir_hmmsearch = (
    "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch"
)

fn_glob = f'{dir_hmmsearch}/*/*.names'
fns = glob(fn_glob)

entry_ids = []
aa_ids = []
dict_gene_names = defaultdict(list)
for fn in fns:
    with open(fn, 'r') as f:
        lines = f.read().splitlines()
    aa_ids += lines

    d, bn = os.path.split(fn)
    d, gene = os.path.split(d)
    entry_id = re.search(rf"(.+)_{gene}.names", bn).group(1)
    entry_id = entry_id.rstrip("_RF")
    entry_ids += [entry_id]*len(lines)

    dict_gene_names[gene] += lines

df = pd.DataFrame(zip(aa_ids, entry_ids), columns=['aa_ids','entry_ids'])
df['source_defline'] = ''

fn_out = "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/entry_id_map_shiri.csv"
df.to_csv(fn_out, header=False, index=False, sep='\t')


dict_gene_clusts = defaultdict(list)
fn_clust_glob = f"{dir_hmmsearch}/*/*.tsv"
fns = glob(fn_clust_glob)
for fn in fns:
    d, bn = os.path.split(fn)
    d, gene = os.path.split(d)
    df = pd.read_csv(fn, header=None, sep='\t')
    dict_gene_clusts[gene].append(df)

dir_out = "/Users/benjamingrodner/work/armbrust/data/prevalence/isolates/hmmsearch"
for gene, dfs in dict_gene_clusts.items():
    df_out = pd.concat(dfs, axis=0)
    fn_clust_out = f"{dir_out}/{gene}_clusters.tsv"
    df_out.to_csv(fn_clust_out, header=False, index=False, sep='\t')

    fn_names_out = f"{dir_out}/{gene}.names"
    names = dict_gene_names[gene]
    with open(fn_names_out, 'w') as f:
        f.write("\n".join(names) + "\n")
