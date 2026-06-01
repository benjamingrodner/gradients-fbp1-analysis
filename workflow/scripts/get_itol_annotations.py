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
import json


def label_unlabeled_nodes(input_jplace, output_jplace, prefix="node_"):
    # 1. Load the jplace file as JSON
    with open(input_jplace, "r") as f:
        data = json.load(f)

    raw_tree_string = data["tree"]

    # 2. Mask the jplace edge tags so ETE3 doesn't choke on them.
    # jplace uses '{0}', '{1}', etc., which violates standard Newick.
    # We will temporarily replace '{number}' with '__EDGE_number__'
    masked_tree_string = re.sub(r"(:\d+\.\d+)\{(\d+)\}", r"__EDGE_\2__\1", raw_tree_string)

    # 3. Parse the tree with ETE3 (format 1 reads internal node names)
    t = Tree(masked_tree_string, parser=1)

    # 4. Traverse and label unnamed internal nodes
    node_counter = 0
    for node in t.traverse("postorder"):
        if not node.is_leaf:
            # Check if it's unnamed (ignoring our edge mask)
            # If node.name only contains the edge mask or is empty
            if node.name is None:
                clean_name = ''
                name = ''
            else:
                clean_name = re.sub(r"__EDGE_\d+__", "", node.name).strip()
                name = node.name

            if not clean_name:
                # Generate a unique name
                while True:
                    potential_name = f"{prefix}{node_counter}"
                    # Ensure we don't accidentally overwrite an existing name
                    if not list(t.search_nodes(name=potential_name)):
                        break
                    node_counter += 1

                # Append the new name to any existing edge mask on that node
                node.name = potential_name + name
                node_counter += 1

    # 5. Write the tree back out to Newick format
    # format=1 preserves internal node names
    updated_tree_string = t.write(parser=1)

    # 6. Unmask the jplace edge tags back to '{number}'
    unmasked_tree_string = re.sub(
        r"__EDGE_(\d+)__(:\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", 
        r"\2{\1}", 
        updated_tree_string
    )

    # 7. Update the JSON object and save
    data["tree"] = unmasked_tree_string

    with open(output_jplace, "w") as f:
        json.dump(data, f, indent=4)


# Example usage:
# label_unlabeled_nodes("input.jplace", "output_labeled.jplace")

def write_fbp1_symbol(master_colors, output_dir):
    """Generates a DATASET_SYMBOL file."""
    output_fn = os.path.join(output_dir, f"fbp1_symbol.txt")
    
    header = f"DATASET_SYMBOL\nSEPARATOR COMMA\nDATASET_LABEL,fbp1\nCOLOR,#000000\nMAXIMUM_SIZE,20\nDATA\n"
    
    with open(output_fn, 'w') as f:
        f.write(header)
        name, color, symbol, size, position = master_colors['fbp1_symbol']
        f.write(f"{name},{symbol},{size},{color},1,{position}\n")


def write_itol_branch(tree, df, master_colors, output_dir):
    """Generates a TREE_COLORS file for branch colors."""
    output_fn = os.path.join(output_dir, f"Branch_treecolors.txt")
    
    header = f"TREE_COLORS\nSEPARATOR COMMA\nDATASET_LABEL,Branch\nDATA\n"
    smap = dict(zip(df["Sequence_ID"].astype(str), df['Source'].astype(str)))
    dmap = dict(zip(df["Sequence_ID"].astype(str), df['Domain'].astype(str)))

    scmap = master_colors['Color_branch']['Source']
    dcmap = master_colors['Color_branch']['Domain']
    
    with open(output_fn, 'w') as f:
        f.write(header)
        for node in tree.leaves():
            n = re.sub('QUERY___','',node.name)
            sval = smap.get(n)
            dval = dmap.get(n)
            if sval is not None:
                scolor = scmap.get(sval)
                dcolor = dcmap.get(dval)
                if scolor is not None:
                    f.write(f"{node.name},branch,{scolor},normal,1\n")
                elif dcolor is not None:
                    f.write(f"{node.name},branch,{dcolor},normal,1\n")
                else:
                    raise ValueError(f'Color is not defined for source {sval} or domain {dval}')
            else:
                raise ValueError(f'No source defined for leaf {n}')


def get_itol_bootstrap_symbol(fn_tree_jplace, tree_support, master_colors, output_dir):
    """Generates a DATASET_SYMBOL file and a new jplace file."""    

    # Generate new file with node labels
    fn_jplace_out = re.sub('.jplace','.nodelabel.jplace',fn_tree_jplace)
    label_unlabeled_nodes(fn_tree_jplace, fn_jplace_out)

    # Load tree to ete
    with open(fn_jplace_out, "r") as f:
        jplace_data = json.load(f)
    raw_tree_string = jplace_data["tree"]
    masked_tree_string = re.sub(r"\{(\d+)\}", r"", raw_tree_string)
    tj = Tree(masked_tree_string, parser=1)
    
    output_fn = os.path.join(output_dir, f"Bootstrap_symbol.txt")

    header = f"DATASET_SYMBOL\nSEPARATOR COMMA\nDATASET_LABEL,Bootstrap\nCOLOR,#000000\nMAXIMUM_SIZE,20\nDATA\n"

    symbol, size, color, position = [
        master_colors['Bootstrap'][i]
        for i in ['symbol', 'size', 'color', 'position']
    ]

    with open(output_fn, 'w') as f:
        f.write(header)
        for nj, ns in zip(tj.traverse(),tree_support.traverse()):
            name = nj.name
            sname = ns.name
            if sname is not None:
                if name != sname:
                    raise ValueError(f"Support tree node label {sname} does not match jplace tree node label {name}")
            supp = ns.support
            if supp is not None:
                if float(supp) >= master_colors['Bootstrap']['min_thresh']:
                    f.write(f'{name},{symbol},{size},{color},1,{position}\n')
    

def create_itol_symbol(tree, leaf_map, color_map, column_name, output_dir):
    """Generates a DATASET_SYMBOL file."""
    output_fn = os.path.join(output_dir, f"{column_name}_symbol.txt")
    
    header = f"DATASET_SYMBOL\nSEPARATOR COMMA\nDATASET_LABEL,{column_name}\nCOLOR,#000000\nMAXIMUM_SIZE,20\nDATA\n"
    
    with open(output_fn, 'w') as f:
        f.write(header)
        for node in tree.leaves():
            n = re.sub('QUERY___','',node.name)
            val = leaf_map.get(n)
            color_sym = color_map.get(val)
            if color_sym is not None:
                color, symbol, size, position = color_sym
                f.write(f"{node.name},{symbol},{size},{color},1,{position}\n")


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

def create_itol_treecolors_range(tree, leaf_map, color_map, column_name, output_dir, base_name, default_color="#ffffff"):
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
    parser.add_argument("-tj", "--tree_jplace", required=True, help="Newick tree file")
    parser.add_argument("-s", "--tree_support", required=True, help="Newick tree file with bootstrap support")
    parser.add_argument("-a", "--ann", required=True, help="CSV annotation file")
    parser.add_argument("-y", "--yaml", required=True, help="YAML color configuration")
    parser.add_argument("-d", "--out", required=True, help="Output directory")
    parser.add_argument("--default", default="#ffffff", help="Default HEX color (default: #ffffff)")
    args = parser.parse_args()

    # Load resources
    try:
        tree = Tree(args.tree)
        tree_support = Tree(args.tree_support)
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
    
    for col in ["Taxon", "Gene", "Domain","Substrate"]:
        if col in df.columns:
            # Mapping Sequence_ID -> Column Value
            leaf_map = dict(zip(df["Sequence_ID"].astype(str), df[col].astype(str)))
            col_colors = master_colors.get(col, {})
            
            # Generate files with the specified default color
            create_itol_strip(tree, leaf_map, col_colors, col, args.out, bn, args.default)
        else:
            print(f"Column '{col}' missing from CSV. Skipping.")

    for col in ["Substrate", "Source","Taxon"]:
        if col in df.columns:
            # Mapping Sequence_ID -> Column Value
            leaf_map = dict(zip(df["Sequence_ID"].astype(str), df[col].astype(str)))
            col_colors = master_colors.get(col, {})
            
            # Generate files with the specified default color
            create_itol_treecolors_range(tree, leaf_map, col_colors, col, args.out, bn, args.default)
        else:
            print(f"Column '{col}' missing from CSV. Skipping.")



    for col in ["Source"]:
        if col in df.columns:
            # Mapping Sequence_ID -> Column Value
            leaf_map = dict(zip(df["Sequence_ID"].astype(str), df[col].astype(str)))
            col_ = col + "_symbol"
            col_colors = master_colors.get(col_)
            if col_colors is None:
                raise ValueError(f"{col_} does not have a defined colormap for symbol annotation.")
            
            # Generate files with the specified default color
            create_itol_symbol(tree, leaf_map, col_colors, col, args.out)
        else:
            print(f"Column '{col}' missing from CSV. Skipping.")
    
    get_itol_bootstrap_symbol(args.tree_jplace, tree_support, master_colors, args.out)

    write_itol_branch(tree, df, master_colors, args.out)

    write_fbp1_symbol(master_colors, args.out)


if __name__ == "__main__":
    main()