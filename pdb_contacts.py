"""
Functions for calculating a PDB distance matrix.
    get_residues(pdbfilename, ['Chain A', 'Chain B'])
        Returns: list of all residues in PDB chains
    calc_distances(residue a, residue b)
        Returns: distance between residue a and residue b
    pdb_map(pdbfilename, ['Chain A', 'Chain B'], distance cutoff)
        Returns:
"""
# test commit text
def get_residues(pdb_fn, chain_ids=None):
    import Bio.PDB
    import os
    import sys
    """Build a simple list of residues from a single chain of a PDB file.
    Args:
        pdb_fn: The path to a PDB file.
        chain_ids: A list of single-character chain identifiers.
    Returns:
        A list of Bio.PDB.Residue objects.
    """

    pdb_id = os.path.splitext(os.path.basename(pdb_fn))[0]

    parser = Bio.PDB.PDBParser(pdb_id, pdb_fn)
    struct = parser.get_structure(pdb_id, pdb_fn)
    model = struct[0]

    if chain_ids is None:
        # get residues from every chain.
        chains = model.get_list()
    else:
        chains = [model[ch_id] for ch_id in chain_ids]

    residues = []

    for ch in chains:
        # make sure res are standard AA
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            if not Bio.PDB.is_aa(res, standard=True):
                sys.stderr.write("WARNING: non-standard AA at %r%s"
                                 % (res.get_id(), os.linesep))
            residues.append(res)

    return residues


def calc_distance(res_a, res_b):
    import numpy as np
    """Calculate the Euclidean distance between a pair of residues
    according to a given distance metric.
    Args:
        res_a: A ``Bio.PDB.Residue`` object.
        res_b: A ``Bio.PDB.Residue`` object.
    Returns:
        The distance between ``res_a`` and ``res_b``.
    """
    # P1: add CB-CB distanct calculation (except for GLY)
    A = res_a["CA"].get_coord()
    B = res_b["CA"].get_coord()
    dist = np.linalg.norm(A - B)
    return dist


def pdb_map(pdbfile, chain_ids, cutoff, ):
    import os
    import numpy as np
    from pandas import read_csv
    from itertools import combinations_with_replacement
    pdbname = os.path.splitext(pdbfile)[0]
    #     cutoff = 10
    # chain_ids = ['C','D']
    # output filename
    filename = pdbname + "_" + str(cutoff) + "A.txt"
    fileout = open(filename, 'w')
    # create list of residues from pdb
    residues = get_residues(pdbfile, chain_ids=chain_ids)
    # make each possible pairs of residues
    pair_indices = combinations_with_replacement(range(len(residues)), 2)
    # write header for output file
    fileout.write("i\tj\td\tchain_1\tchain_2\n")
    for i, j in pair_indices:
        res_a = residues[i]
        res_b = residues[j]
        # get chain id
        chain_a = res_a.get_parent().id
        chain_b = res_b.get_parent().id
        dist = calc_distance(res_a, res_b)
        if cutoff >= dist > 0:
            fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i + 1, j + 1, dist, chain_a, chain_b))
    # makes a pandas dataframe
    df_pdb = read_csv(filename, delim_whitespace=True)
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
    return df_pdb, df_mon, df_inter
