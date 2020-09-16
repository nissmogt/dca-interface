def split_header_seq(msa_file, len_a, limit=-1):
    import numpy as np
    """Function to split fasta headers and sequences"""
    header_a = []
    header_b = []
    seq_a = []
    seq_b = []
    lines = open(msa_file, "r")
    # next(lines)  # skips null template
    # next(lines)
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
            seq_a.append(line[:len_a])  # sequence A
            seq_b.append(line[len_a:])  # sequence B

    lines.close()
    return np.array(header_a), np.array(header_b), np.array(seq_a), np.array(seq_b)


def permute_index(n_seqs, n_replicates):
    import numpy as np
    # creates 2 lists of random indices for seq A and B
    for seed in range(n_replicates):
        r1 = np.random.RandomState(seed)
        r2 = np.random.RandomState(seed + 2)
        yield r1.permutation(n_seqs), r2.permutation(n_seqs)


def scramble_sequence(msa_file, len_a, n_replicates):
    import os
    """Randomly pair sequences"""
    msa_name = os.path.basename(msa_file)
    header_a, header_b, seq_a, seq_b = split_header_seq(msa_file, len_a)
    n_seqs = len(seq_b)
    # creates 2 lists of random indices for seq A and B
    index = list(permute_index(n_seqs, n_replicates))
    outfile = []
    for rep in range(n_replicates):
        scramble_seq = []
        scramble_header = []
        for i in range(n_seqs):
            scramble_header.append(header_a[index[rep][0][i]] + '_' + header_b[index[rep][1][i]])
            scramble_seq.append(seq_a[index[rep][0][i]] + seq_b[index[rep][1][i]])
        scramble_msa = dict(zip(scramble_header, scramble_seq))

        # Write MSA replicates to file
        outfile.append(('MSA_rep%d_scrambled_' + msa_name) % rep)
        with open(outfile[rep], 'w') as f:
            for key in scramble_msa.keys():
                f.write(">%s\n%s\n" % (key, scramble_msa[key]))
    return outfile


import shannon
import numpy as np
nr = 1
len_a = 112
id = '1EM8'
c1 = 'D'
c2 = 'C'
# msa_directory = 'PDB_benchmark_alignments\\'
msa_directory = 'Sequences\\'
msa = msa_directory + '{}_{}_{}_{}.fas'.format(id, c1, id, c2)
out = scramble_sequence(msa, len_a, nr)

# Calc Shannon Entropy
# aln, l_seq, ind = shannon.parseMSA(msa, alnformat='fasta', verbose=1)
# la, gl = shannon.shannon_entropy_list_msa(aln)

# Compute avg of entropy of scrambled seqs
# array = []
# for i in range(nr):
#     scramble_msa = "MSA_rep{}_scrambled_{}_{}_{}_{}.fas".format(i, id, c1, id, c2)
#     saln, sl_seq, sind = shannon.parseMSA(scramble_msa, alnformat='fasta', verbose=1)
#     sla, sgl = shannon.shannon_entropy_list_msa(saln)
#     array.append(sla)
# avg_ent = np.average(array, axis=0)
#
# import ttest
# ttest.shapiro_test(la, avg_ent)
# p_t = ttest.f_test(la, avg_ent)
# ttest.t_test(la, avg_ent, p_t, 0.05)
# shannon.plot(ind, la, sla, verbose=1, msa_name='scramble_unscrambled_{}_{}_{}_{}'.format(id, c1, id, c2))
