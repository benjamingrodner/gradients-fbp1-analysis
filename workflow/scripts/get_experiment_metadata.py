# framework written by Gemini 3 flash free tier 5/17/26
## Prompt
# Please write a python script taking click arguments of filenames for a counts parquet file, a config yaml, a cluster tsv table, a sample info text file, a string indicating the experiment, and an output parquet. Load the parquets and the csv using ibis with duckb backend. Load the yaml config as well. 

import os
import click
import pandas as pd
import yaml
from collections import defaultdict
import re


def parse_sample_info(fn_sample_info, cfg, experiment):
    dict_config = cfg[experiment]
    if fn_sample_info is None: # Did not download reads from NCBI
        dict_meta = dict_config.get('sample_metadata')
        cols = dict_config.get('colnames_metadata')
        if (dict_meta is None) | (cols is None):
            raise ValueError(f"Metadata not defined for experiment {experiment}, with fn_sample_info {fn_sample_info}")
        return pd.DataFrame.from_dict(dict_meta, orient='index', columns=cols)

    else: 
        dict_meta = defaultdict(dict)
        regex = dict_config['regex_condition']  # condition string matching 
        subs = dict_config.get('subs_condition')  # list of substitutions for conditions
        if dict_config['method_counts'] == 'salmon-srr':  # Get info from the srr_info table
            df_info = pd.read_csv(fn_sample_info)
            col = dict_config['srr_info_col']
            for sam, info in df_info[['Run',col]].values:
                match = re.search(regex, info)
                if match is not None:
                    cond = match.group('condition')
                    if subs is not None:
                        for s in subs:
                            cond = re.sub(s[0],s[1],cond)
                    sam_col = 'count_' + sam
                    dict_meta['condition'][sam_col] = cond
        else: # Get info from the Biosample_info table
            with open(fn_sample_info, 'r') as f:
                sam, cond = None, None
                for line in f: # Parse each line for matches
                    m_sam = re.search(r'(?<=BioSample: )SAMN\d+',line)
                    match = re.search(regex, line)
                    if m_sam is not None:
                        sam = m_sam.group(0)
                        sam_col = 'count_' + sam
                    if match is not None:
                        cond = match.group('condition')
                    if (sam is not None) & (cond is not None): # if you collect two matches, write them out
                        dict_meta['condition'][sam_col] = cond
                        sam, cond = None, None
        return pd.DataFrame(dict_meta)


def load_config(config_path: str) -> dict:
    """Loads and returns the YAML configuration."""
    try:
        with open(config_path, "r") as f:
            return yaml.safe_load(f)
    except Exception as e:
        raise click.ClickException(f"Error reading config YAML: {e}")


@click.command()
@click.option(
    "--config",
    "-g",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the configuration YAML file.",
)
@click.option(
    "--sample_info",
    "-s",
    type=click.Path(exists=True, dir_okay=False),
    required=False,
    default=None,
    help="Path to the sample information text file.",
)
@click.option(
    "--experiment",
    "-e",
    type=str,
    required=True,
    help="String identifier indicating the experiment.",
)
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Path where the output Parquet file will be saved.",
)
def main(
    config: str, sample_info: str, experiment: str, output: str
):

    # 1. Load the YAML configuration
    click.echo(f"Loading configuration from: {config}")
    cfg = load_config(config)


    out_table = parse_sample_info(sample_info, cfg, experiment)

    click.echo(f"Writing output to: {output}")
    try:
        # Ensure the output directory exists
        out_dir = os.path.dirname(output)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        out_table.to_csv(output)
        click.echo("Processing complete!")

    except Exception as e:
        raise click.ClickException(f"Error writing output Parquet: {e}")


if __name__ == "__main__":
    main()