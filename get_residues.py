def get_residues(msa, seq=False):
    """
    Build a simple list of residues from a single chain of a PDB file.
    :param msa: String - MSA filename
    :param seq: Boolean (Default: False) - Outputs sequence if True.
    :return: A list of Bio.PDB.Residue objects.
    """
    import Bio.PDB
    import sys
    import os
    from three2one import three2one

    pdb_path = "PDB_benchmark_structures\\"
    msa_name = msa.strip(".fas")
    pdbfile = "{}.pdb".format(msa_name[:4])
    chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]
    chain_lengths = []

    fname = "(get_residues)"
    print("\t{}\tprocessing {} file...".format(fname, pdbfile))
    # if pdb_id == 'pdb':
    parser = Bio.PDB.PDBParser()
    # else:
    #     parser = Bio.PDB.MMCIFParser()

    struct = parser.get_structure(pdbfile.strip(".pdb"), pdb_path + pdbfile)
    model = struct[0]
    # if len(self.chain_ids) == 0:
    # get residues from every chain.
    #    chains = model.get_list()
    # else:
    chains = [model[ch_id] for ch_id in chain_ids]
    print("\t{} CHAIN IDs:\t{}".format(fname, chain_ids))

    residues = []
    sequence = []
    for ch in chains:
        # make sure res are standard AA
        num_residues = 0
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            # if Bio.PDB.is_aa(res, standard=True):
            is_regular_res = res.has_id('CA') and res.has_id('O')
            res_id = res.get_id()[0]
            if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
                residues.append(res)
                sequence.append(res.get_resname())
                num_residues += 1
            else:
                sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                 (res.get_id(), os.linesep))
        chain_lengths.append(num_residues)  # P1: BUG if runs twice, chain_length list will be twice as long.

    if seq:
        sequence = three2one(sequence)
        seq_a = sequence[:chain_lengths[0]]
        seq_b = sequence[chain_lengths[0]:]
        return seq_a, seq_b
    else:
        return residues, chain_lengths
