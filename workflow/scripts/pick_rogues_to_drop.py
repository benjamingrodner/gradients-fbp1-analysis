# Generated with Gemini 3 Flash, operating in the Free tier
# 29 April 2026
# using the following prompt:

###
# "Write a Python script that identifies the 'elbow' point of taxon improvement from a RogueNaRok output file. The script should meet the following requirements:

# Input Handling:
# Accept a RogueNaRok_droppedRogues file and output filename as command-line arguments using argparse.
# Handle the file as a space-delimited or whitespace-delimited table.
# Standardize column headers to be case-insensitive to ensure 'taxon' and 'rawimprovement' are always found.
# Mathematical Logic (Elbow Detection):
# Calculate the 'elbow' or 'knee' of the rawimprovement curve using the geometric distance-to-chord method.
# This involves defining a line between the first and last data points and finding the point on the curve with the maximum orthogonal distance from that line.
# Handle edge cases where the file might have fewer than 3 rows.
# Output:
# The script should define a main function that write a list of the taxa names (one per line) up to and including the identified elbow point to the defined output file.
# Include an optional --verbose flag to display the count of taxa found.
# Dependencies:
# Use pandas for data manipulation and numpy for vector calculations.
# Please include clean error handling for missing files or malformed columns."
###

# Edited to debug

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import sys
import os

def find_elbow_point(data):
    """
    Finds the elbow point using the geometric distance-to-chord method.
    The elbow is the point with the maximum orthogonal distance from 
    the line connecting the first and last points.
    """
    n_points = len(data)
    if n_points < 3:
        # Not enough points to form a curve/elbow
        return n_points - 1

    # Coordinates of all points
    # x is the index (0 to N), y is the raw improvement value
    coords = np.column_stack((np.arange(n_points), data['score_sorted'].values))

    first_point = coords[0]
    last_point = coords[-1]

    # Vector representing the line (chord) from first to last point
    line_vec = last_point - first_point
    line_vec_norm = line_vec / np.sqrt(np.sum(line_vec**2))

    # Vector from first point to all points
    vec_from_first = coords - first_point

    # Scalar projection of vec_from_first onto line_vec
    scalar_proj = np.dot(vec_from_first, line_vec_norm)
    
    # Orthogonal vector: vec_from_first - parallel_projection
    vec_to_line = vec_from_first - np.outer(scalar_proj, line_vec_norm)

    # Distances are the norms of the orthogonal vectors
    distances = np.sqrt(np.sum(vec_to_line**2, axis=1))

    return np.argmax(distances)

def main():
    parser = argparse.ArgumentParser(
        description="Get the subset of RogueNaRok rogues worth dropping."
    )
    parser.add_argument("input", help="Path to the RogueNaRok_droppedRogues file")
    parser.add_argument("output", help="Path for the output list of taxa")
    parser.add_argument("-f", "--frac_of_total_improvement", default='elbow', help="If removing all of the rogues gives improvement of 1, what fraction of that improvement do you want. Default mode finds the elbow of the curve (where diminishing returns start).")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print the count of taxa found")

    args = parser.parse_args()

    # 1. Load Data
    try:
        # RogueNaRok files are usually whitespace-delimited
        df = pd.read_csv(args.input, sep=r'\s+', engine='python')
    except FileNotFoundError:
        print(f"Error: File '{args.input}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    # 2. Standardize Columns (Case-insensitive)
    df.columns = [c.lower() for c in df.columns]
    
    if 'taxon' not in df.columns or 'rawimprovement' not in df.columns:
        print("Error: Required columns 'taxon' and 'rawImprovement' not found.")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)

    if df.empty:
        print("Error: Input file is empty.")
        sys.exit(1)

    # remove nan rows for rawimprovement and taxon
    df = df[~df['taxon'].isna() & ~df['rawimprovement'].isna()]

    # Ensure data is sorted by improvement (RogueNaRok usually sorts descending)
    df = df.sort_values(by='rawimprovement', ascending=False).reset_index(drop=True)

    # Get a proxy for the score now that improvements are sorted
    score = 0
    score_sorted = []
    for imp in df['rawimprovement'].values:
        score += imp
        score_sorted.append(score)
    score_sorted = np.array(score_sorted) / np.max(score_sorted)
    df['score_sorted'] = score_sorted
    
    # 3. Pick index to filter by
    if args.frac_of_total_improvement == 'elbow':
        elbow_idx = find_elbow_point(df)
    else:
        frac = float(args.frac_of_total_improvement)
        elbow_idx = np.argmin(np.abs(score_sorted - frac))

    # Plot choice
    fig, ax = plt.subplots(figsize=[5,5])
    ax.plot(np.arange(df.shape[0]), score_sorted)
    ax.plot([0,df.shape[0]],[score_sorted[elbow_idx]]*2)
    ax.plot([elbow_idx]*2,[0,1])
    ax.set_ylabel('Frac of improvement')
    ax.set_xlabel('Number of rogues removed')
    plt.savefig(args.output + '.png')

    # 4. Filter and Output
    # "Up to and including the elbow point"
    rogues_to_drop = df.loc[:elbow_idx, 'taxon']

    try:
        with open(args.output, 'w') as f:
            for taxon in rogues_to_drop:
                f.write(f"{taxon}\n")
    except IOError as e:
        print(f"Error writing to output file: {e}")
        sys.exit(1)

    if args.verbose:
        print(f"Total taxa processed: {len(df)}")
        print(f"Elbow index identified: {elbow_idx}")
        print(f"Taxa included in drop-set: {len(rogues_to_drop)}")

if __name__ == "__main__":
    main()