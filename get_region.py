def get_dca_indices(msa_name, length_a=None):
    alignment = "PDB_benchmark_alignments\\a2m\\{}.a2m".format(msa_name)
    f = open(alignment, 'r')
    header = f.readline()
    seqs = f.readlines()
    target_sequence = []
    for i, seq in enumerate(seqs):
        if seq[0] != '>':
            target_sequence.append(seq.strip('\n'))
        else:
            break
    target_sequence = ''.join(target_sequence)
    if length_a:
        seq1 = target_sequence[:length_a]
        seq2 = target_sequence[length_a:]
        _, low_1 = get_lowercase_index(seq1)
        _, low_2 = get_lowercase_index(seq2)
        return [seq1, seq2], [len(seq1) - len(low_1), len(seq2) - len(low_2)], [low_1, low_2]
    else:
        idx, _ = get_lowercase_index(target_sequence)
        return idx


def get_lowercase_index(template_sequence):
    """
    Finds index for lowercase characters in any sequence.
    :param template_sequence:
    :return:
    """
    idx = []
    lower_idx = []
    count = 0
    for i, ch in enumerate(template_sequence):
        if not ch.islower():
            count += 1
            idx.append(count)
        else:
            lower_idx.append(i)

    return idx, lower_idx
