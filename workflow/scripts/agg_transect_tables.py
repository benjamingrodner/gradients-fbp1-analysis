# Framework of produced by Gemini 3.5 Flash free tier 12 Jun 2026
### Prompt
# Please write a python script with click arguments:
# - directory of input csv files
# - input config yaml
# - output file for processed csv file
# - output path with basename but no extension for plots
# That does the following:
# - If the output csv does not exist: For each input csv takes the values in a certain column, e.g. "Fe", and makes a new "tidy" table where the column name is entered as a value in a new column named "target" and the column value as a value in the new column named "concentration". Appends the tables together and writes the table to the output csv filename
# - Using the appended dataframe produces a plot where concentration values are grouped by lattitude, cruise, and target, averaged over a range of pre-defined latitude bins (from the config file), and plotted with dots for the means with means within the same cruise connected by lines (points ordered by ascending latitude value) and error bars for standard deviations. 
###

# Then edited by Ben Grodner

import os
import click
import yaml
import pandas as pd


def load_config(config_path):
    """Loads the YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def process_csv_files(output_csv, cfg, tzcf_csv):
    """
    Reads all CSVs in input_dir, extracts the target_col, 
    reshapes into a tidy format, and appends them together.
    """
    # csv_files = glob.glob(os.path.join(input_dir, "*.csv"))
    # if not csv_files:
    #     raise FileNotFoundError(f"No CSV files found in directory: {input_dir}")
    
    all_dfs = []
    
    for _, dict_cfg in cfg.items():
        file_path = dict_cfg['fn_data']
        df = pd.read_csv(file_path)

        dr, bn = os.path.split(file_path)

        dict_col = dict_cfg['columns']
        # if bn in cfg:
        #     dict_cfg = cfg[bn]['columns']
        # else:
        #     raise ValueError(f"File {bn}, which is present in the input dir {dr}, has no defined params in the config file")
        
        required_cols = list(dict_col.values())

        # Verify columns ex
        
        for col in required_cols:
            if col not in df.columns:
                raise ValueError(f"Missing required column '{col}' in {file_path}")
        
        # Keep only necessary columns for the tidy conversion

        df_filtered = df[required_cols].copy()
        
        # Rename target column to concentration, and create the 'target' label column
        dict_rename = {col: ncol for ncol, col in dict_col.items()}
        df_filtered = df_filtered.rename(columns=dict_rename)
        df_tidy = []
        for t, units in dict_cfg['targets'].items():
            df_filtered['target'] = t
            df_filtered['value'] = df[t]
            df_filtered['units'] = units
            df_tidy.append(df_filtered)


        # Set the cruise column to be consistent with previous naming
        if 'cruise' not in df_filtered.columns:
            df_tidy['cruise'] = dict_cfg['cruise']
        else:
            dict_crs = dict_cfg['cruise_rename']
            cruise_rename = {val: nval for nval, val in dict_crs.items()}
            df_tidy['cruise'] = df_tidy['cruise'].map(cruise_rename)


        # Get dist from tzcf column
        df_tzcf = pd.read_csv(tzcf_csv)
        map_tzcf = dict(zip(df_tzcf['cruise'], df_tzcf['lat_tzcf']))
        df_tidy['lat_tzcf'] = df_tidy['cruise'].map(map_tzcf)
        df_tidy['dist_tzcf'] = df_tidy['lat'] - df_tidy['lat_tzcf']

        
        all_dfs.append(df_filtered)
        
    if not all_dfs:
        raise ValueError("No data was successfully processed from the input CSVs.")
        
    # Combine and save
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df.to_csv(output_csv, index=False)
    click.echo(f"Successfully created processed data file: {output_csv}")
    return

@click.command()
@click.option('--config', '-c', type=click.Path(exists=True, dir_okay=False), required=True,
              help='Path to the YAML configuration file.')
@click.option('--tzcf-csv', '-t', type=click.Path(dir_okay=False), required=True,
              help='Path to a CSV with "cruise" and "lat_tzcf" columns.')
@click.option('--output-csv', '-o', type=click.Path(dir_okay=False), required=True,
              help='Output path for the processed tidy CSV file.')
@click.option('--plot-base', '-p', type=click.Path(dir_okay=False), required=True,
              help='Output path base name (no extension) for the generated plots.')
def main(config, tzcf_csv, output_csv, plot_base):
    """
    Processes oceanographic CSV data, bins concentrations by latitude, 
    and generates an aggregated line/error bar plot.
    """
    # 1. Load Configuration
    cfg = load_config(config)['agg_transect_tables']
    process_csv_files(output_csv, cfg, tzcf_csv)

if __name__ == '__main__':
    main()





