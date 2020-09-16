def combine_chain_sequence(msa_name, split=False):
    """
    Concatenates PDB sequences
    # writes file of pdbid and length of first chain and pdb fasta file
    :return: str - combined chain 1 and chain 2 protein sequence from pdb
    """
    from read_pdb import read_pdb
    fname = "(combine_pdb_seq)"
    print("{}\tcat-ting sequences...".format(fname))
    print(fname + "\t\tFasta file\tNumber of chains")

    chain_ids = [(msa_name.strip(".fas")).split("_")[1], (msa_name.strip(".fas")).split('_')[3]]
    pdb_chain, pdb_fasta_seq = read_pdb(msa_name)
    header_seq_dict = dict(zip(pdb_chain, pdb_fasta_seq))

    # For each chain in msa find relevant pdb seq and cat the two seqs
    first_seq = header_seq_dict[chain_ids[0]]
    second_seq = header_seq_dict[chain_ids[1]]
    full_seq = first_seq + second_seq
    print("{}\t\tPDB seq length: {}".format(fname, len(full_seq)))
    print("\tFinished cat-ting sequences.".format(fname))
    if split:
        return first_seq, second_seq
    else:
        return full_seq
