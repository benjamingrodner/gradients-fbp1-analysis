# Framework written with gemini 3 in the free tier 26 Oct 2026
## prompt:
# Please write a python script that loads a config file from yaml, iterates through a list of parquet files stored in the yaml, and reads each parquet with ibis. Specify inputs with argparse including config path, output path, as well as memory and number of cores for ibis
# Edited by Ben grodner


import argparse
from pathlib import Path
import ibis
import pandas as pd
import yaml
import re
import numpy as np


def load_config(config_path: str | Path) -> dict:
  """Load configuration parameters from a YAML file."""
  path = Path(config_path)
  if not path.is_file():
    raise FileNotFoundError(f"Configuration file not found: {path}")

  with open(path, "r", encoding="utf-8") as f:
    return yaml.safe_load(f)

def load_meta(meta_path: str | Path) -> dict:
  """Load configuration parameters from a YAML file."""
  path = Path(meta_path)
  if not path.is_file():
    raise FileNotFoundError(f"Metadata file not found: {path}")

  return pd.read_csv(path)


def parse_args() -> argparse.Namespace:
  """Parse command-line arguments."""
  parser = argparse.ArgumentParser(
      description=(
          "Process a list of Parquet files defined in a YAML config using Ibis."
      )
  )
  parser.add_argument(
      "-c",
      "--config",
      type=Path,
      default=Path("config.yaml"),
      help="Path to the YAML configuration file (default: config.yaml).",
  )
  parser.add_argument(
      "-o",
      "--output",
      type=Path,
      required=True,
      help="Path to save the output summary file (e.g., CSV).",
  )
  parser.add_argument(
      "--cores",
      type=int,
      default=None,
      help="Number of CPU cores/threads to allocate for DuckDB execution.",
  )
  parser.add_argument(
      "--memory",
      type=str,
      default=None,
      help="Memory limit for DuckDB (e.g., '4GB', '512MB').",
  )
  return parser.parse_args()


def main() -> None:
  args = parse_args()

  # Load configuration
  config = load_config(args.config)
  fn_meta = config['metadata_metat']
  meta = load_meta(fn_meta)
  batches = config.get("batches", [])
  if not batches:
    print("No batches specified")
    return
  fmt_parq = config['dirs_data']['fmt_metatranscriptome_annotations']
#   parquet_files = [fmt_parq.format(batch=b) for b in batches]
      

  # Build connection configuration dictionary dynamically
  connect_kwargs = {}
  if args.cores is not None:
    connect_kwargs["threads"] = args.cores
  if args.memory is not None:
    connect_kwargs["memory_limit"] = args.memory
  
  # Initialize the Ibis backend (DuckDB)
  con = ibis.duckdb.connect(**connect_kwargs)
  results = []

  for batch in batches:
    file_path_str = fmt_parq.format(batch=batch)
    file_path = Path(file_path_str)

    if not file_path.is_file():
      print(f"[Warning] File not found or invalid path: {file_path}")
      continue

    print(f"\n[Processing] {file_path}")

    try:
      # Read the Parquet file lazily into an Ibis table expression
      table = con.read_parquet(str(file_path))

      # Extract metadata
      fns_kallisto = meta['kallisto_fn'].values
      cols_sams = [c for c in table.columns if c in fns_kallisto]
      sid = meta.loc[meta['kallisto_fn'].isin(cols_sams),'SampleID']
      sams_w_size = np.unique([re.sub(r'\.\w+$', '', s) for s in sid])
      sams = np.unique([re.sub(r'\.0\.2um|\.3um', '', s) for s in sams_w_size])


      results.append({
          "batch": str(batch),
          "n_w_filt_rep": len(cols_sams),
          "n_w_filt": len(sams_w_size),
          "n_sams": len(sams),
      })

    except Exception as e:
      print(f"[Error] Failed to process {file_path}: {e}")

  # Save results to the specified output file
  if results:
    df_out = pd.DataFrame(results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(args.output, index=False)
    print(f"\n[Success] Summary results written to {args.output}")
  else:
    print("\n[Warning] No results generated to write out.")


if __name__ == "__main__":
  main()