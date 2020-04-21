#!/home/kmm5/anaconda2/envs/py36/bin/python
"""
Name: process_scores.py
Author: Kareem Mehrabiani
Date: 12th February 2020

Creates a contact map from DCA raw scoring matrix. Optional: Maps DCA residues
to PDB.

Usage: process_scores.py -s [dca scoring matrix] 
       If mapping pfam to pdb indices add:
       --hmmscan [hmmscan file]

optional arguments:
  -h, --help         show this help message and exit
  -s S               Raw DCA scores output.
  --hmmscan HMMSCAN  hmmscan file (if mapping)

Outputs:
    - three-column contact map (in same dir as input)
    - mapping reference (optional)
    - mapped three-column contact map (in current dir)

"""
import argparse
import os
import pandas as pd
import numpy as np
import sys
import re


def process_scan(scan_in):
    """Separates hmmscan results into Domain A and B."""
    f = open(scan_in, 'r')
    data = f.read()
    raw_data = (re.findall(r'== domain (.*?)Internal', data, re.DOTALL))
    len(raw_data)

    scan_A = raw_data[0].split('\n')
    if len(raw_data) > 1:
        scan_B = raw_data[1].split('\n')

    domain_A = scan_A[1].split()
    pdb_A = scan_A[3].split()
    domain_B = scan_B[1].split()
    pdb_B = scan_B[3].split()
    domain_seq = domain_A[2] + domain_B[2]
    pdb_seq = pdb_A[2] + pdb_B[2]
    return [domain_A, domain_B], [pdb_A, pdb_B], domain_seq, pdb_seq


###############################################################################################
def map_loop(scan_in):
    """Creates a dictionary object with DCA-to-PDB map."""
    print("\tMapping using %s file\n" % scan_in)
    # Initialize
    file_header = 'Domain\td_seq\tp_seq\tPDB'
    domain, pdb, domain_seq, pdb_seq = process_scan(scan_in)
    domain_A, domain_B = domain[0], domain[1]
    pdb_A, pdb_B = pdb[0], pdb[1]
    domain_begin_A = int(domain_A[1])
    pdb_begin_A = int(pdb_A[1])
    length_domain = len(domain_seq)
    length_pdb = len(pdb_seq)
    domain_num = np.zeros(length_domain, dtype=int)
    pdb_num = np.zeros(length_pdb, dtype=int)
    resid_domain = []
    resid_pdb = []
    insert_count = 0
    gap_count = 0

    # Loop through domain and pdb sequence, counting every non-insert
    for i, res in enumerate(domain_seq, domain_begin_A):
        resid_domain.append(res)
        if res == '.':
            insert_count = insert_count + 1
        if res != '.':
            domain_num[i - int(domain_A[1])] = i - insert_count
    for i, res in enumerate(pdb_seq, pdb_begin_A):
        resid_pdb.append(res)
        if i > int(pdb_A[3]):
            if res == '-':
                gap_count = gap_count + 1
            if res != '-':
                pdb_num[i - int(pdb_A[1])] = (i + 4) - gap_count
        else:
            if res == '-':
                gap_count = gap_count + 1
            if res != '-':
                pdb_num[i - int(pdb_A[1])] = i - gap_count

    name = os.path.split(scan_in)[-1]
    outfile_ref = "map_ref_" + name
    map_array = np.transpose([domain_num, resid_domain, resid_pdb, pdb_num])
    np.savetxt(outfile_ref, map_array, fmt='%s', delimiter='\t', header=file_header)

    map_dict = dict(zip(domain_num[domain_num != 0] - (domain_begin_A - 2),
                        pdb_num[pdb_num != 0] - (domain_begin_A - 2)))
    print("Output (current dir):\n\t%s" % outfile_ref)
    return map_dict


###############################################################################################
def map_dca(scan_in, dcafile):
    """Converts DCA pairs into PDB coordinates and outputs a file."""
    names = ['i', 'j', 'score']
    dca_data = pd.read_csv(dcafile, names=names, delim_whitespace=True, usecols=(0, 1, 2))
    name = os.path.split(scan_in)[-1]
    outfile_map = "map_" + name
    map_dict = map_loop(scan_in)
    N = len(dca_data)
    mapped_dca = []
    for i in range(N):
        dca_resi = dca_data['i'][i]
        dca_resj = dca_data['j'][i]
        fn = dca_data['score'][i]
        if dca_resi in map_dict and dca_resj in map_dict:
            mapped_dca.append([map_dict[dca_resi], map_dict[dca_resj], fn, dca_resi, dca_resj])

    header = 'i\tj\tfn\ti_pfam\tj_pfam'
    np.savetxt(outfile_map, mapped_dca, fmt='%d\t%d\t%f\t%d\t%d', header=header)
    print("\t%s" % outfile_map)


###############################################################################################
def cm_make(scores):
    # x is the output FN matrix
    filename = scores.strip(".txt") + "_CM.txt"
    x = np.loadtxt(scores)
    x_output = []
    N = x.shape[0]
    for i in range(N - 1):
        for j in range(i + 1, N):
            x_output.append([i + 1, j + 1, x[i, j]])

    header = ['i', 'j', 'score']
    df_x = pd.DataFrame(x_output, columns=header)
    df_x = df_x.sort_values(ascending=False, by=['score'])
    np.savetxt(filename, df_x, fmt='%d\t%d\t%f')
    print("\tWrote contact map to: %s" % filename)
    return filename


###############################################################################################
def main():
    print("Starting")
    parser = argparse.ArgumentParser(description='Creates a contact map from \
            DCA raw scoring matrix. Optional: Maps DCA residues to PDB.')
    parser.add_argument('-s', required=True, help="Raw DCA scores output.")
    parser.add_argument('--hmmscan', required=False, help='hmmscan file \
            (if mapping)')
    args = parser.parse_args()
    scores = args.s

    dcafile = cm_make(scores)

    if args.hmmscan:
        scan_in = args.hmmscan
        map_dca(scan_in, dcafile)


if __name__ == '__main__':
    main()
