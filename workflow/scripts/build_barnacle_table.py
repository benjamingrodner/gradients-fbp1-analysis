# skeleton produced using gemini 3.5 flash in free tier 5/31/26
### prompt
# Please write a python script with click arguments:
# - a filename to a csv table with metadata
# - a filename to another csv table with normalization factors
# - a string that is a key to decide which function to use to parse the table
# - an output filename

# Steps:
# - load the csvs to dataframes using pandas
# - use the string to select the correct function
# - execute the function on the dataframes which produces a  new dataframe
# - write the new dataframe to the output filename
###

import sys
import click
import pandas as pd


# -----------------------------------------------------------------------------
# Parsing Functions
# -----------------------------------------------------------------------------
def parse_method_a(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    """Example parsing function A: Multiplies metadata by normalization factors."""
    click.echo("Executing Parsing Method A...")
    # Dynamic example logic: standardizing data
    # (Replace this with your actual business logic)
    return df_meta.multiply(df_norm, fill_value=1)


def parse_method_b(df_meta: pd.DataFrame, df_norm: pd.DataFrame) -> pd.DataFrame:
    """Example parsing function B: Joins tables or performs alternative logic."""
    click.echo("Executing Parsing Method B...")
    # Dynamic example logic: adding dataframes
    # (Replace this with your actual business logic)
    return df_meta.add(df_norm, fill_value=0)


# -----------------------------------------------------------------------------
# Function Registry
# -----------------------------------------------------------------------------
# This maps your string keys to the actual Python functions
FUNCTION_REGISTRY = {
    "method_a": parse_method_a,
    "method_b": parse_method_b,
}

# -----------------------------------------------------------------------------
# Click CLI Command Configuration
# -----------------------------------------------------------------------------
@click.command()
@click.option(
    "--meta-file",
    "-m",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the CSV file containing metadata.",
)
@click.option(
    "--norm-file",
    "-n",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Path to the CSV file containing normalization factors.",
)
@click.option(
    "--method-key",
    "-k",
    type=click.Choice(list(FUNCTION_REGISTRY.keys()), case_sensitive=False),
    required=True,
    help="Key deciding which parsing function to execute.",
)
@click.option(
    "--output-file",
    "-o",
    type=click.Path(dir_okay=False, writable=True),
    required=True,
    help="Path where the output CSV will be saved.",
)
def main(meta_file, norm_file, method_key, output_file):
    """
    Load metadata and normalization CSVs, apply a selected parsing function
    based on a method key, and save the resulting DataFrame to an output file.
    """
    try:
        # 1. Load the CSVs into pandas DataFrames
        click.echo(f"Loading files:\n - Metadata: {meta_file}\n - Norm Factors: {norm_file}")
        df_meta = pd.read_csv(meta_file)
        df_norm = pd.read_csv(norm_file)

        # 2. Select the correct function using the string key
        # (click.Choice already validates that the key exists in our dict)
        parse_function = FUNCTION_REGISTRY[method_key.lower()]

        # 3. Execute the function
        result_df = parse_function(df_meta, df_norm)

        # 4. Write the new dataframe to the output filename
        click.echo(f"Saving results to: {output_file}")
        result_df.to_csv(output_file, index=False)
        
        click.echo("Processing complete successfully!")

    except Exception as e:
        click.echo(f"Error occurred during execution: {e}", err=True)
        sys.exit(1)


if __name__ == "__main