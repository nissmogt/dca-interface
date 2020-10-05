def read_msa(msa_name, split=False, len_a=None):
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


def parse_fasta(msa_name, null=True, limit=-1):
    """
    function to parse fasta
    :return header and sequence of fasta
    """
    import numpy as np
    header = []
    sequence = []
    msa_path = "PDB_benchmark_alignments\\"

    lines = open("{}{}".format(msa_path, msa_name), 'r')
    if null:
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
