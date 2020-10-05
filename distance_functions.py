def calc_min_dist(res_a, res_b):
    import numpy as np
    na_atoms = len(res_a.get_list())
    nb_atoms = len(res_b.get_list())
    atomic_matrix = np.zeros((na_atoms, nb_atoms))
    for i, atom_i in enumerate(res_a):
        for j, atom_j in enumerate(res_b):
            atomic_matrix[i][j] = np.linalg.norm(atom_i.get_coord() - atom_j.get_coord())
    mindist = np.amin(atomic_matrix)
    atom_i_index, atom_j_index = np.where(atomic_matrix == mindist)    # gives index of mindist atoms i,j
    min_atom_i = res_a.get_list()[atom_i_index[0]]
    min_atom_j = res_b.get_list()[atom_j_index[0]]
    return mindist, (min_atom_i.id, min_atom_j.id)


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
