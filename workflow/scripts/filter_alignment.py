#!/usr/bin/env python3

import argparse
import numpy as np


def fasta_to_dict(fasta_file):
    """
    Reads a FASTA file and returns a dictionary of {header: sequence}.

    Parameters:
        fasta_file (str): Path to the FASTA file.

    Returns:
        dict: Keys are headers (without '>'), values are sequences as strings.
    """
    fasta_dict = {}
    header = None
    seq_chunks = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # skip empty lines
            if line.startswith(">"):
                # Save the previous entry
                if header:
                    fasta_dict[header] = ''.join(seq_chunks)
                header = line[1:].strip()  # remove '>'
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # Don't forget the last entry
        if header:
            fasta_dict[header] = ''.join(seq_chunks)

    return fasta_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fn_aln", type=str, required=True, help="alignment file")
    parser.add_argument("-o", "--fn_filt", type=str, required=True, help="Output file")
    parser.add_argument("-f", "--frac_range", type=float, default=0.25, help="Within what fraction of the mean in either direction can be kept? e.g. 0.25 for mean length 100 keeps sequences between 125 and 75 bp")
    args = parser.parse_args()


    fn_outside = f'{args.fn_filt}.removed'
    dict_header_seq = fasta_to_dict(args.fn_aln)
    lens = []
    for s in dict_header_seq.values():
        ngaps = s.count('-')
        slen = len(s) - ngaps
        lens.append(slen)
    len_mean = np.mean(lens)
    t_low = len_mean - (len_mean*args.frac_range)
    t_high = len_mean + (len_mean*args.frac_range)
    print(t_low, len_mean, t_high)

    dict_out = {}
    outsides = []
    for h, s in dict_header_seq.items():
        ngaps = s.count('-')
        slen = len(s) - ngaps
        bool_l = slen < t_low
        bool_h = slen > t_high
        # print(slen, bool_l, bool_h)
        if bool_l | bool_h:
            outsides.append(h)
        else:
            dict_out[h] = s
    with open(args.fn_filt, 'w') as f:
        for h, s in dict_out.items():
            f.write(f'>{h}\n{s}\n')
    with open(fn_outside, 'w') as f:
        for h in outsides:
            f.write(f'{h}\n')
    print(f'Filtered {len(outsides)}')
    print(f'Kept {len(dict_header_seq) - len(outsides)}')
    return

if __name__ == "__main__":
    main()