def msa_first_seq_reader(msa_name, split=False, len_a=None):
    """
    Reads first sequence in MSA which should be a template
    :return: List - MSA sequence
    """
    # Create MSA-seq-template file object
    msa_path = "PDB_benchmark_alignments\\"
    msa_file_obj = open("{}{}".format(msa_path, msa_name), 'r')
    msa_header = msa_file_obj.readline().rstrip()
    msa_seq = msa_file_obj.readline().rstrip()
    msa_file_obj.close()
    if split:
        return msa_seq[:len_a], msa_seq[len_a:]
    else:
        return msa_seq
