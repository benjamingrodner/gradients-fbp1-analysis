# Generated with Gemini 3 Flash, operating in the Free tier 
# 30 Apr 2026
# using the following prompt:

# "Write a Python script that generates iTOL (Interactive Tree Of Life) annotation files—specifically DATASET_COLORSTRIP and TREE_COLORS—using the ete4 library for tree traversal.
# Inputs:
# -t (Tree): A Newick format tree file.
# -a (Annotations): A CSV file where leaf names are in a column named Sequence_ID.
# -y (YAML): A configuration file containing a dictionary dict_cstrip_val_color. This dictionary maps column names (e.g., 'Taxon', 'Gene') to sub-dictionaries of {value: hex_color}.
# -d (Directory): The output directory for generated files.
# --default: A command-line argument for a default hex color (defaulting to #ffffff).
# Logic Requirements:

# Target Columns: The script should generate separate files for the columns Taxon and Gene.
# Colorstrips: For each leaf in the tree, find its value in the CSV and assign the color from the YAML. If the value is missing from the YAML, use the default color.
# Tree Colors (Ranges): * Iterate through the tree using preorder traversal.
# For each node, check if all descendant leaves share the same value for the target column.
# If they do, write a range entry for that node in the TREE_COLORS file and skip its descendants to avoid redundant labeling.
# Use the provided default color if the uniform value is not found in the YAML.
# Formatting: Use standard iTOL headers (SEPARATOR COMMA, DATASET_LABEL, etc.). Output filenames should follow the pattern: {tree_basename}.{type}_{column}.txt."

# Edited for debugging and 

import os
import argparse
import pandas as pd
import yaml
from ete4 import Tree
import re

def create_itol_strip(tree, leaf_map, color_map, column_name, output_dir, base_name, default_color="#ffffff"):
    """Generates a DATASET_COLORSTRIP file."""
    output_fn = os.path.join(output_dir, f"{column_name}_colorstrip.txt")
    
    header = f"DATASET_COLORSTRIP\nSEPARATOR COMMA\nDATASET_LABEL,{column_name}\nCOLOR,#000000\nDATA\n"
    
    with open(output_fn, 'w') as f:
        f.write(header)
        for node in tree.leaves():
            n = re.sub('QUERY___','',node.name)
            val = leaf_map.get(n)
            # Use color_map if exists, otherwise use default_color
            color = color_map.get(val, default_color)
            label = val if val else "Unknown"
            f.write(f"{node.name},{color},{label}\n")

def create_itol_treecolors(tree, leaf_map, color_map, column_name, output_dir, base_name, default_color="#ffffff"):
    """Generates a TREE_COLORS file for background ranges."""
    output_fn = os.path.join(output_dir, f"{column_name}_treecolors.txt")
    
    header = f"TREE_COLORS\nSEPARATOR COMMA\nDATASET_LABEL,{column_name}_range\nDATA\n"
    
    nodes_labeled = set()
    
    with open(output_fn, 'w') as f:
        f.write(header)
        for node in tree.traverse("preorder"):
            if node.name is not None:
                if node.name in nodes_labeled:
                    continue
                
                leaves = [l.name for l in node.leaves()]
                vals = [
                    leaf_map.get(re.sub('QUERY___','',l)) 
                    for l in leaves
                ]
                
                # If all leaves share the same value
                if len(set(vals)) == 1:
                    val = vals[0]
                    color = color_map.get(val, default_color)
                    
                    f.write(f"{node.name},range,{color},{val if val else 'Default'}\n")
                    
                    # Prevent redundant labeling of children
                    for desc in node.descendants():
                        nodes_labeled.add(desc.name)

def main():
    parser = argparse.ArgumentParser(description="Generate iTOL files with default color options")
    parser.add_argument("-t", "--tree", required=True, help="Newick tree file")
    parser.add_argument("-a", "--ann", required=True, help="CSV annotation file")
    parser.add_argument("-y", "--yaml", required=True, help="YAML color configuration")
    parser.add_argument("-d", "--out", required=True, help="Output directory")
    parser.add_argument("--default", default="#ffffff", help="Default HEX color (default: #ffffff)")
    args = parser.parse_args()

    # Load resources
    try:
        tree = Tree(args.tree)
        df = pd.read_csv(args.ann)
        with open(args.yaml, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading files: {e}")
        return

    master_colors = config.get('dict_cstrip_val_color', {})
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    bn = os.path.basename(args.tree)
    
    for col in ["Taxon", "Gene", "Domain"]:
        if col in df.columns:
            # Mapping Sequence_ID -> Column Value
            leaf_map = dict(zip(df["Sequence_ID"].astype(str), df[col].astype(str)))
            col_colors = master_colors.get(col, {})
            
            # Generate files with the specified default color
            create_itol_strip(tree, leaf_map, col_colors, col, args.out, bn, args.default)
        else:
            print(f"Column '{col}' missing from CSV. Skipping.")

    for col in ["Substrate", "Source"]:
        if col in df.columns:
            # Mapping Sequence_ID -> Column Value
            leaf_map = dict(zip(df["Sequence_ID"].astype(str), df[col].astype(str)))
            col_colors = master_colors.get(col, {})
            
            # Generate files with the specified default color
            create_itol_treecolors(tree, leaf_map, col_colors, col, args.out, bn, args.default)
        else:
            print(f"Column '{col}' missing from CSV. Skipping.")


if __name__ == "__main__":
    main()