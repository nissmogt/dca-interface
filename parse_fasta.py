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
