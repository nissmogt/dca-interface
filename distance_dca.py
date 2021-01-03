
def distance_dca(df_dca, msa_name, atom, other_name=False):
    """
    Calculates distance matrix.
    :param other_name:
    :param atom:
    :param msa_name:
    :param df_dca:
    :return: Three Dataframe objects.
    """
    import numpy as np
    import time
    from get_residues import get_residues
    from distance_functions import calc_ca_distance, calc_min_dist

    # PARAMS
    fname = "(DCA DISTANCE) "
    print(fname)
    pdb_path = "PDB_benchmark_structures\\"
    chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]

    # output list of residues from pdb
    residues, chain_lengths = get_residues(msa_name)

    # make each possible pairs of residues
    n_pairs = len(df_dca)
    # pair_indices = combinations_with_replacement(range(len(residues)), 2)
    pair_indices = np.array([df_dca["i"], df_dca["j"], df_dca["score"]]).transpose()
    start_time = time.time()
    distance_list = []
    chain_1_list = []
    chain_2_list = []
    atom_id_list = []
    residue_list = []

    for i, j, score in pair_indices:
        # subtract index by one because index for residues begins at 0 and i,j start at 1
        res_a = residues[int(i)-1]
        res_b = residues[int(j)-1]
        # get chain id and append to list
        chain_1_list.append(res_a.get_parent().id)
        chain_2_list.append(res_b.get_parent().id)
        residue_list.append((res_a.resname, res_b.resname))
        if atom == 'aa':
            mindist, atom_ids = calc_min_dist(res_a, res_b)
            distance_list.append(mindist)
            atom_id_list.append(atom_ids)
        else:
            if res_a.has_id("CA") and res_b.has_id("CA"):
                distance_list.append(calc_ca_distance(res_a, res_b))
            else:
                print("{} NOTE! Res {} \n\tor {} not calculated! (missing CA)\n".format(fname, res_a.get_full_id(),
                                                                                        res_b.get_full_id()))
    # fileout.close()
    print("{}\t -- MAIN LOOP TIME -- {}".format(fname, time.time() - start_time))

    # Saves results to file depending on cutoff type (all atom or c-alpha)
    if not other_name:
        if atom == 'aa':
            df_new = df_dca.assign(dist=distance_list, chain_1=chain_1_list, chain_2=chain_2_list,
                                   resnames=residue_list, atom_id=atom_id_list)
            # # outfile = "results\\FN_{}_inter_mapped_aa_dist_top{}.txt".format(msa_name, n_pairs)
            # outfile = "results\\FN_{}_mapped_aa_dist_top{}.txt".format(msa_name, n_pairs)
            # np.savetxt(outfile, df_new, header="i\tj\tscore\tdist_aa\tchain_1\tchain_2\tresnames\tatom_id",
            #            fmt='%d\t%d\t%f\t%f\t%s\t%s\t%s\t%s', comments='')
        else:
            df_new = df_dca.assign(dist=distance_list, chain_1=chain_1_list,
                                   chain_2=chain_2_list, resnames=residue_list)
            # outfile = "results\\FN_{}_inter_mapped_ca_dist_top{}.txt".format(msa_name, n_pairs)
            outfile = "results\\FN_{}_mapped_ca_dist_top{}.txt".format(msa_name, n_pairs)
            np.savetxt(outfile, df_new, header="i\tj\tscore\tdist_ca\tchain_1\tchain_2\tresnames",
                       fmt='%d\t%d\t%f\t%f\t%s\t%s\t%s', comments='')
    else:
        # outDir = "scrambled_results\\from_averaging_jmatrices\\"
        df_new = df_dca.assign(dist=distance_list, chain_1=chain_1_list, chain_2=chain_2_list,
                               resnames=residue_list, atom_id=atom_id_list)
        # outfile = "{}FNi_apc_{}_inter_mapped_aa_dist_top{}.txt".format(outDir, msa_name, n_pairs)
        # np.savetxt(outfile, df_new, header="i\tj\tscore\tdist_aa\tchain_1\tchain_2\tresnames\tatom_id",
        #            fmt='%d\t%d\t%f\t%f\t%s\t%s\t%s\t%s', comments='')

    return df_new


def read_dca_distance_matrix(msa_name, n_pairs, atom, other_name=False):
    from pandas import read_csv
    if not other_name:
        if atom:
            filename = "results\\FN_{}_inter_mapped_aa_dist_top{}.txt".format(msa_name, n_pairs)
            df_dca = read_csv(filename, delimiter="\t")
        else:
            filename = "results\\FN_{}_inter_mapped_ca_dist_top{}.txt".format(msa_name, n_pairs)
            df_dca = read_csv(filename, delimiter="\t")
    else:
        # filename = "scrambled_results\\FNi_apc_{}_inter_mapped_aa_dist_top{}.txt".format(msa_name, n_pairs)
        filename = "scrambled_results\\fni_matrices\\FNi_apc_after_{}_inter_mapped_aa_dist_top{}.txt".format(msa_name, n_pairs)
        df_dca = read_csv(filename, delimiter="\t")

    return df_dca


def vectorize_dca_contacts(df_dca, pdb_total_length):
    import numpy as np
    """
    Converts DCA pairs into binary matrix and then flattens into 1-D array
    :param df_dca: Dataframe object - DCA pairs NOTE: Should already be sliced to desired length
    :param dimer_length: int - Length of dimer
    :return: 1-D flat array
    """
    import pandas
    # Initialize contact matrices
    dca_matrix = np.zeros((pdb_total_length, pdb_total_length))

    dca_array = df_dca.iloc[:, :2].to_numpy()  # convert column i and j to numpy array
    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in dca_array:
        dca_matrix[int(i)-1, int(j)-1] = 1
    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    dca_flat_matrix = dca_matrix.ravel()
    return dca_flat_matrix


