def cat_dimer(msa_file, limit=-1):
    msa_dir = "Sequences\\"
    msa_name = "PF01746_1UAL.txt_filtered_25p"
    import numpy as np
    """Function to split fasta headers and sequences"""
    header = []
    seq = []
    lines = open(msa_dir + msa_name, "r")
    for idx, line in enumerate(lines):
        line = line.rstrip()  # removes whitespace from the right
        if line[0] == '>':
            header_entry = line[1:].split('_')
            header.append(header_entry[1].split('/')[0] + "_" + header_entry[1].split('/')[0])
        else:
            len_a = len(line)
            seq.append(line + line)
    lines.close()
    new_dimer = dict(zip(header, seq))
    outfile = "{}dimer_{}".format(msa_dir, msa_name)
    with open(outfile, 'w') as f:
        for key in new_dimer.keys():
            f.write(">%s\n%s\n" % (key, new_dimer[key]))
    return outfile
