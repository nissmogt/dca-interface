def calc_min_dist(res_a, res_b):
    import numpy as np
    dist = []
    for atom_i in res_a:
        for atom_j in res_b:
            dist.append(np.linalg.norm(atom_i.get_coord() - atom_j.get_coord()))
    dist_array = np.array(dist)
    mindist = min(dist_array[np.nonzero(dist_array)])
    return mindist


def calc_ca_distance(res_a, res_b):
    """
    Calculates the distance between a pair of CA atoms
    :param res_a: Biopython residue object - residue a
    :param res_b: Biopython residue object - residue b
    :return: Distance between CA atoms
    """
    import numpy as np
    a = res_a["CA"].get_coord()
    b = res_b["CA"].get_coord()
    dist = np.linalg.norm(a - b)
    return dist


def distance_dca(msa, n_pairs, cutoff, all_atom=False):
    """
    Calculates distance matrix.
    :param n_pairs:
    :param all_atom:
    :return: Three Dataframe objects.
    """
    from Bio.PDB import is_aa
    import numpy as np
    from pandas import read_csv
    from itertools import combinations_with_replacement
    import time
    from get_residues import get_residues

    # PARAMS
    fname = "(DCA DISTANCE) "
    pdb_path = "PDB_benchmark_structures\\"
    msa_name = msa.strip(".fas")
    chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]
    all_atom = False

    if all_atom:
        filename = "{}dca_heavy_atom_distance_matrix_{}_{}A.txt".format(pdb_path, msa_name, cutoff)
    else:
        filename = "{}dca_distance_matrix_{}_ca{}A.txt".format(pdb_path, msa_name, cutoff)

    fileout = open(filename, 'w')
    fileout.write("i\tj\tscore\td\tchain_1\tchain_2\n")

    # LOAD DCA CONTACTS AND GET LIST OF ALL RESIDUES IN STRUCTURE
    dca_residues = np.loadtxt("results\\mapped_cm_FNi_apc_{}.txt".format(msa_name))
    residues = get_residues(msa)

    # THIS CAN BE REMOVED IN the FUTURE
    chain1_residues = []
    chain2_residues = []
    for res in residues:
        if res.get_parent().id == chain_ids[0]:
            if is_aa(res, standard=True):
                chain1_residues.append(res)
        elif res.get_parent().id == chain_ids[1]:
            if is_aa(res, standard=True):
                chain2_residues.append(res)
    # make each possible pairs of residues
    # pair_indices = combinations_with_replacement(range(len(residues)), 2)
    pair_indices = dca_residues[:n_pairs]
    start_time = time.time()
    for i, j, score in pair_indices:
        print(i, j)
        res_a = residues[int(i)]
        res_b = residues[int(j)]
        # res_a = chain1_residues[int(i)]
        # res_b = chain2_residues[int(j) - len(chain1_residues)]
        # get chain id
        if all_atom:
            chain_a = res_a.get_parent().id
            chain_b = res_b.get_parent().id
            mindist = calc_min_dist(res_a, res_b)
            if mindist <= cutoff:
                fileout.write("%d\t%d\t%f\t%s\t%s\n" % (int(i), int(j), mindist, chain_a, chain_b))
        else:
            if res_a.has_id("CA") and res_b.has_id("CA"):
                chain_a = res_a.get_parent().id
                chain_b = res_b.get_parent().id
                dist = calc_ca_distance(res_a, res_b)
                if cutoff >= dist > 0.0:
                    fileout.write("%d\t%d\t%f\t%f\t%s\t%s\n" % (int(i), int(j), score, dist, chain_a, chain_b))
                else:
                    print("{} NOTE! Res {} \n\tor {} not calculated! (missing CA)\n".format(fname,
                                                                                            res_a.get_full_id(),
                                                                                            res_b.get_full_id()))
    fileout.close()
    print("{}\t -- MAIN LOOP TIME -- {}".format(fname, time.time() - start_time))
    # Load data into a pandas dataframe
    df_pdb = read_csv(filename, delim_whitespace=True)
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
    return df_inter
