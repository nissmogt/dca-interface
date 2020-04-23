"""
Functions for calculating a PDB distance matrix.
    get_residues(pdbfilename, ['Chain A', 'Chain B'])
        Returns: list of all residues in PDB chains
    calc_distances(residue a, residue b)
        Returns: distance between residue a and residue b
    pdb_map(pdbfilename, ['Chain A', 'Chain B'], distance cutoff)
        Returns:
"""


# dev branch test
def get_residues(pdbfile, chain_ids=None):
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

    pdb_id = os.path.splitext(os.path.basename(pdbfile))[0]
    if os.path.basename(pdbfile).split('.')[-1] == 'pdb':
        print('.pdb file...')
        parser = Bio.PDB.PDBParser()
        struct = parser.get_structure(pdb_id, pdbfile)
        model = struct[0]
    else:
        print('.cif file...')
        parser = Bio.PDB.MMCIFParser()
        struct = parser.get_structure(pdb_id, pdbfile)
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


def calc_ha_distance(res_a, res_b):
    """Calculate the Euclidean distance between a pair of heavy atoms
    Args:
        res_a: A ``Bio.PDB.Residue`` object.
        res_b: A ``Bio.PDB.Residue`` object.
    Returns:
        The distance between ``res_a`` and ``res_b``.
    """
    import numpy as np
    min_dist = 8.0
    for a in res_a.get_iterator():
        for b in res_b.get_iterator():
            coord_a = a.get_coord()
            coord_b = b.get_coord()
            dist = np.linalg.norm(coord_a - coord_b)
            if dist < min_dist:
                min_dist = dist

    return min_dist


def calc_ca_distance(res_a, res_b):
    import numpy as np
    """Calculate the Euclidean distance between a pair of residues
    according to a given distance metric.
    Args:
        res_a: A ``Bio.PDB.Residue`` object.
        res_b: A ``Bio.PDB.Residue`` object.
    Returns:
        The distance between ``res_a`` and ``res_b``.
    """
    a = res_a["CA"].get_coord()
    b = res_b["CA"].get_coord()
    dist = np.linalg.norm(a - b)
    return dist


def pdb_map(pdbfile, chain_ids, cutoff, ):
    import os
    from pandas import read_csv
    from itertools import combinations_with_replacement
    import time

    pdbname = os.path.splitext(pdbfile)[0]
    # output filename
    filename = pdbname + "_" + str(cutoff) + "A.txt"
    fileout = open(filename, 'w')
    # write header for output file
    fileout.write("i\tj\td\tchain_1\tchain_2\n")

    # create list of residues from pdb
    residues = get_residues(pdbfile, chain_ids=chain_ids)

    # make each possible pairs of residues
    pair_indices = combinations_with_replacement(range(len(residues)), 2)
    start_time = time()
    for i, j in pair_indices:
        res_a = residues[i]
        res_b = residues[j]
        # get chain id
        chain_a = res_a.get_parent().id
        chain_b = res_b.get_parent().id
        dist = calc_ha_distance(res_a, res_b)
    #    if cutoff >= dist > 0:
        fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i + 1, j + 1, dist, chain_a, chain_b))

    print("Loop time: %s" % (time() - start_time))
    # makes a pandas dataframe
    df_pdb = read_csv(filename, delim_whitespace=True)
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]

    return df_pdb, df_mon, df_inter


