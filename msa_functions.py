import time
import numpy as np


def read_first_sequence_in_msa(msa_name, split=False, len_a=None):
    """
    Reads first sequence in MSA which should be a template
    :return: List - MSA sequence
    """
    # Create MSA-seq-template file object
    msa_dir = "PDB_benchmark_alignments\\"
    msa_file_obj = open("{}{}.fas".format(msa_dir, msa_name), 'r')
    msa_header = msa_file_obj.readline().rstrip()
    msa_seq = msa_file_obj.readline().rstrip()
    msa_file_obj.close()
    if split:
        return msa_seq[:len_a], msa_seq[len_a:]
    else:
        return msa_seq


def parse_fasta(msa_name, skip_first=True, limit=-1):
    """
    function to parse fasta
    :return header and sequence of fasta
    """
    header = []
    sequence = []
    msa_path = "PDB_benchmark_alignments\\"

    lines = open("{}{}".format(msa_path, msa_name), 'r')
    if skip_first:
        # used to skip first fasta sequence and header
        skip_null = [next(lines) for x in range(2)]
    for line in lines:
        line = line.rstrip()
        if line[0] == ">":
            if len(header) == limit:
                break
            header.append(line[1:])
            sequence.append([])
        else:
            sequence[-1].append(line)
    lines.close()
    sequence = [''.join(seq) for seq in sequence]
    return np.array(header), np.array(sequence)


def effective_seq_calc(sequence, eff_cutoff=0.8):
    """compute effective weight for each sequence"""
    from scipy.spatial.distance import pdist, squareform

    # pairwise identity
    msa_sm = 1.0 - squareform(pdist(sequence, "hamming"))

    # weight for each sequence
    msa_w = (msa_sm >= eff_cutoff).astype(np.float)
    msa_w = 1 / np.sum(msa_w, -1)

    return msa_w


def msa_object(msaName, skip_first=True):
    """converts list of sequences to msa"""
    from aa2num import aa2num

    headers, sequences = parse_fasta(msaName, skip_first=skip_first)

    msa_in_number_form = []
    for seq in sequences:
        msa_in_number_form.append([aa2num(aa) for aa in seq])
    msa_in_number_form = np.array(msa_in_number_form)

    # remove positions with more than > 50% gaps
    # msa, v_idx = self.filt_gaps(msa_ori, 0.5)
    msa = msa_in_number_form

    # compute effective weight for each sequence
    msa_weights = effective_seq_calc(msa, 0.8)
    neff = np.sum(msa_weights)
    ncol = msa.shape[1]  # length of sequence
    # w_idx = v_idx[np.stack(np.triu_indices(ncol, 1), -1)]

    # "msa_ori": msa_ori,
    #     "msa": msa,
    #     "weights": msa_weights,
    # "ncol_ori": msa_ori.shape[1]}
    # "v_idx": v_idx,
    # "w_idx": w_idx
    # }
    return {"neff": neff, "nrow": msa.shape[0], "ncol": ncol}


def count_sequence(msaName):
    """
    Reads msa fasta file and counts number of '>' characters which equals the number of sequences
    :param msaName:
    :return: int Number of sequences in fasta file
    """
    msaDir = "PDB_benchmark_alignments\\"
    with open("{}{}".format(msaDir, msaName), "r", encoding='utf-8') as msa:
        count = msa.read().count(">")
    return count


def stats_for_msa_list(msa_list):
    """
    For every msa in msa list count number of sequences, calculate number of effective sequences,
    and add these values to a list with the msa name.
    :param msa_list:
    :return: list containing MSA name, number of sequences, number of effective sequences
    """
    stats_list = []
    with open(msa_list, "r", encoding='utf-8') as msaFile:
        start = time.time()
        for msa in msaFile:
            msaName = msa.rstrip()
            print(msaName)
            msaObject = msa_object(msaName, skip_first=True)
            nSeqs = msaObject["nrow"]
            nEffective = msaObject["neff"]
            seqLength = msaObject["ncol"]
            stats_list.append([msaName.strip('.fas'), nSeqs, nEffective, seqLength])
        print("Duration of loop: {}".format(time.time() - start))
    return stats_list


def write_msa_stats(msa_list):
    import csv
    stats_list = stats_for_msa_list(msa_list)
    outfile = "stats_{}.csv".format(msa_list.strip(".txt"))
    with open(outfile, "w", encoding='utf-8', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(["MSA", "Nseqs", "Neffective"])
        writer.writerows(stats_list)
    return stats_list


# msafiles = "benchmark_dimer_systems.txt"
# s = write_msa_stats(msafiles)

def make_monomer_msa_from_dimer(msa):
    from scramble_sequence import split_header_seq
    from read_db import get_lengths
    from get_region import get_dca_indices
    outDir = "monomer_alignments\\"
    # msa = "2OXG_Z_2OXG_Y"
    cid = [msa.split("_")[1], msa.split("_")[3]]
    ch = get_lengths(msa)
    _, dca_chains, _ = get_dca_indices(msa, length_a=ch[0])
    x = split_header_seq(msa, dca_chains[0])
    for i in range(2):
        a = np.array([x[i], x[i+2]])
        a = a.transpose()
        np.savetxt("{}{}_{}.fas".format(outDir, msa[:4], cid[i]), a, fmt=">%s\n%s")

dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
      "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
      '5MU7_B_5MU7_A']
msaname = dimers[0]
# make_monomer_msa_from_dimer(msaname)