def read_pdb(pdb_id):
    """
    Reads PDB file and extracts relevant chain sequence.
    :return: Two lists - Chain ids and sequence
    """
    fname = "(read_pdb)"
    pdb_path = "PDB_benchmark_structures\\"
    pdb_fasta_file = "{}{}.fasta".format(pdb_path, pdb_id)
    pdb_fasta_obj = open(pdb_fasta_file, 'r')
    num_of_chains = int(pdb_fasta_obj.read().count('>'))
    pdb_fasta_obj.seek(0)  # return to top of file

    pdb_chains = []
    pdb_seq = []
    # Loop thru each chain and append header and seq to a list, then make into lookup table
    for line in range(num_of_chains):
        pdb_chains.append(pdb_fasta_obj.readline().rstrip()[-1])
        pdb_seq.append(pdb_fasta_obj.readline().rstrip())
    print("{}\t\t{}\tnum of chains: {}".format(fname, pdb_fasta_file, num_of_chains))
    pdb_fasta_obj.close()

    # Error check
    if len(pdb_seq) != len(pdb_chains):
        print("Error! Number of seqs not equal to number of chains\n".format(fname))
    return pdb_chains, pdb_seq
