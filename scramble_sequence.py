import os
import numpy as np


def split_header_seq(msa_name, chain_length):
    """
    Function to split fasta headers and sequences

    :param msa_name:
    :param chain_length:
    :return:
    """
    header_a = []
    header_b = []
    seq_a = []
    seq_b = []
    msadir = "PDB_benchmark_alignments\\"
    # msadir = "PDB_benchmark_alignments\\filtered\\"
    lines = open("{}{}.fas".format(msadir, msa_name), "r")
    next(lines)  # skips null template
    next(lines)
    for idx, line in enumerate(lines):
        line = line.rstrip()  # removes whitespace from the right
        if line[0] == '>':
            if len(header_a) == -1:
                break
            header_entry = line[1:].split('_')
            # header_entry = line[1:].split('|')
            if len(header_entry) < 3:
                header_a.append(header_entry[0])
                header_b.append(header_entry[1])
            else:
                header_a.append(header_entry[1])
                header_b.append(header_entry[2])
        else:
            seq_a.append(line[:chain_length])  # sequence A
            seq_b.append(line[chain_length:])  # sequence B

    lines.close()
    return np.array(header_a), np.array(header_b), np.array(seq_a), np.array(seq_b)


def permute_index(n_seqs, n_replicates):
    # creates 2 lists of random indices for seq A and B
    for seed in range(n_replicates):
        R1 = np.random.RandomState(seed)
        R2 = np.random.RandomState(seed + 2)
        yield R1.permutation(n_seqs), R2.permutation(n_seqs)


def scramble_sequence(msa_name, n_replicates):
    from get_region import get_dca_indices
    from read_db import get_lengths
    """

    :param msa_name:
    :param n_replicates:
    :return:
    """
    if msa_name[0] == '1':
        pdbid_start_number = 1
    elif msa_name[0] == '2':
        pdbid_start_number = 2
    elif msa_name[0] == '3':
        pdbid_start_number = 3
    elif msa_name[0] == '4':
        pdbid_start_number = 4
    elif msa_name[0] == '5':
        pdbid_start_number = 5

    results_dir = "scrambled_sequences_nots\\pdbid_{}\\{}\\".format(pdbid_start_number, msa_name)
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    print("\tScramble {}".format(msa_name))
    uniprot_lengths = get_lengths(msa_name)
    _, chain_length, _ = get_dca_indices(msa_name, uniprot_lengths[0])
    print(chain_length[0], uniprot_lengths[0])
    header_a, header_b, seq_a, seq_b = split_header_seq(msa_name, chain_length[0])
    nSeqs = len(seq_b)
    # creates 2 lists of random indices for seq A and B
    randomIndex = list(permute_index(nSeqs, n_replicates))
    outfile = []
    for rep in range(n_replicates):
        scramble_seq = []
        scramble_header = []
        for i in range(nSeqs):
            rand_index_1 = randomIndex[rep][0][i]
            rand_index_2 = randomIndex[rep][1][i]
            scramble_header.append(header_a[rand_index_1] + '_' + header_b[rand_index_2])
            scramble_seq.append(seq_a[rand_index_1] + seq_b[rand_index_2])
        scramble_msa_dict = dict(zip(scramble_header, scramble_seq))
    #     Write MSA replicates to file
        outfile.append('{}{}_rep{}_scrambled.fas'.format(results_dir, msa_name, rep))
        with open(outfile[rep], 'w', encoding='utf-8') as f:
            for key in scramble_msa_dict.keys():
                f.write(">{}\n{}\n".format(key, scramble_msa_dict[key]))
    return outfile


# nr = 5
# dimers = ['5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A']
# dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
#       "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
#       '5MU7_B_5MU7_A']
# pdbid = '1M56'
# c1 = 'C'
# c2 = 'D'
# msa = '{}_{}_{}_{}_filtered_25'.format(pdbid, c1, pdbid, c2)
# out = scramble_sequence(dimers[5], nr)
