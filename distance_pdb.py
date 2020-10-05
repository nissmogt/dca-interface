def distance_matrix(msa_name, all_atom=False):
    """
    Calculates distance matrix.
    :param df_dca:
    :param msa_name:
    :param all_atom:
    :return: Three Dataframe objects.
    """
    # msa = "3RPF_A_3RPF_D.fas"
    # cutoff = 12
    import time
    import pandas as pd
    import numpy as np
    from itertools import combinations_with_replacement
    from get_residues import get_residues
    from distance_functions import calc_min_dist, calc_ca_distance

    fname = "(distance matrix)"
    out_path = "distance_matrix\\"
    chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]

    # output list of residues from pdb
    residues, chain_lengths = get_residues(msa_name)
    # make each possible pairs of residues NOTE: INDEX BEGINS AT 0 BUT 1 IS ADDED BELOW
    pair_list = combinations_with_replacement(range(len(residues)), 2)
    start_time = time.time()
    resi_list = []
    resj_list = []
    actual_i_list = []
    actual_j_list = []
    distance_list = []
    chain_1_list = []
    chain_2_list = []
    atom_id_list = []
    residue_list = []
    print("\t{} begin loop".format(fname))

    for i, j in pair_list:
        if i != j:    # ensure residue i not equal to j
            res_a = residues[int(i)]
            res_b = residues[int(j)]
            actual_i_list.append(res_a.id[1])
            actual_j_list.append(res_b.id[1])
            # get chain id
            chain_1_list.append(res_a.get_parent().id)
            chain_2_list.append(res_b.get_parent().id)
            # resets res index to 1 SEE NOTE ABOVE.
            resi_list.append(i+1)
            resj_list.append(j+1)
            residue_list.append((res_a.resname, res_b.resname))
            if all_atom:
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
    print("{}\t -- LOOP TIME -- {}".format(fname, time.time() - start_time))
    # makes a pandas dataframe
    if all_atom:
        df_pdb = pd.DataFrame({'i': resi_list, 'j': resj_list, 'd': distance_list,
                               'si': actual_i_list, 'sj': actual_j_list, 'chain_1': chain_1_list,
                               'chain_2': chain_2_list, 'resnames': residue_list, 'atom_id': atom_id_list})
        filename = "{}heavy_atom_distance_matrix_{}.txt".format(out_path, msa_name)
        np.savetxt(filename, df_pdb, header="i\tj\td\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id",
                   fmt='%d\t%d\t%f\t%d\t%d\t%s\t%s\t%s\t%s', comments='')
    else:
        df_pdb = pd.DataFrame({'i': resi_list, 'j': resj_list, 'd': distance_list,
                               'si': actual_i_list, 'sj': actual_j_list,
                               'chain_1': chain_1_list, 'chain_2': chain_2_list, 'resnames': residue_list})
        filename = "{}ca_distance_matrix_{}.txt".format(out_path, msa_name)
        np.savetxt(filename, df_pdb, header="i\tj\td\tsi\tsj\tchain_1\tchain_2\tresnames",
                   fmt='%d\t%d\t%f\t%d\t%d\t%s\t%s\t%s', comments='')

    print("{} wrote {}".format(fname, filename))
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
    return df_pdb, df_mon, df_inter, chain_lengths


def read_distance_matrix_file(msa_name, all_atom=False):
    import os
    import pandas as pd
    fname = "(read dist mat)"
    out_path = "distance_matrix\\"
    c1 = msa_name.split("_")[1]
    c2 = msa_name.split("_")[-1]
    if all_atom:
        filename = "{}heavy_atom_distance_matrix_{}.txt".format(out_path, msa_name)
    else:
        filename = "{}ca_distance_matrix_{}.txt".format(out_path, msa_name)

    if os.stat(filename).st_size != 0:                          # checks if file is not empty
        print("\t{} reading from {}".format(fname, filename))
        df_pdb = pd.read_csv(filename, delimiter="\t")
        df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
        df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
        total_length = max(df_mon[df_mon["chain_2"] == c2].j)
        chain_1_length = max(df_mon[df_mon["chain_1"] == c1].j)
        chain_2_length = total_length - chain_1_length
    else:
        df_pdb, df_mon, df_inter, chain_lengths = distance_matrix(msa_name, all_atom)
    return df_mon, df_inter, [chain_1_length, chain_2_length]


def plot_cm(pdb_df_list, cutoff, length_a, length, atom, df_dca, msa_name=None, other_dca=None):
    """
    :param df_dca:
    :param pdb_df_list:
    :param msa_name:
    :param cutoff:
    :param length_a:
    :param length:
    :param atom:
    :param other_dca:
    :return:
    """
    import matplotlib.pylab as plt
    from numpy.random import randint

    print("\t\tPlotting pdb monomer and interface...")
    count = ord(msa_name[1]) * randint(100) + ord(msa_name[2]) * randint(101, 1000)
    df_mon = pdb_df_list[0]
    df_inter = pdb_df_list[1]

    # Plotting
    fig = plt.figure(num=count, figsize=(10, 10), dpi=100)
    ax = fig.add_subplot(1, 1, 1)

    # monomer
    mon_cutoff = 8
    ax.scatter('i', 'j', data=df_mon[df_mon["d"] <= mon_cutoff],
               label='{} monomer, {}$\AA$'.format(msa_name[:4], mon_cutoff),
               c='xkcd:grey', marker='o', s=25, alpha=0.7)
    # interface
    ax.scatter('i', 'j', data=df_inter[df_inter["d"] <= cutoff],
               label='{} dimer, {}$\AA$'.format(msa_name[:4], cutoff),
               c='xkcd:azure', marker='o', s=25, alpha=0.7)

    # Vanilla DCA
    if not df_dca.empty:
        ax.scatter('i', 'j', data=df_dca, label='Vanilla DCA interface top {}'.format(len(df_dca)), c='black', s=30)

    if not other_dca.empty:
        ax.scatter('i', 'j', data=other_dca, label='DCAi scrambled top {}'.format(len(other_dca)),
                   c='red', s=30, alpha=0.6)

    # Plot dimer separator line
    extend = 2
    ax.hlines(length_a, 0, length, linestyles=(0, (5, 10)), alpha=0.5, colors='black')
    ax.vlines(length_a, 0, length+extend, linestyles=(0, (5, 10)), alpha=0.5, colors='black')

    # plot design
    ax.legend(loc='lower right')
    ax.set_facecolor('ivory')
    if atom == 'aa':
        plt.title("PDB: {} (aa)".format(msa_name))
    else:
        plt.title("PDB: {} (c$\alpha$)".format(msa_name))
    plt.minorticks_on()
    plt.grid(alpha=0.3)
    plt.xlabel("residue i"), plt.ylabel("residue j")
    ax.grid(which='major', alpha=0.4, c='black')
    ax.grid(which='minor', linestyle=':', alpha=0.5, c='gray')
    if not df_dca.empty:
        img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\plmDCA_scrambled\\"
        # imgname = "{}{}_top{}_{}{}.png".format(img_dir, msa_name, len(df_dca), atom, cutoff)
        # imgname = "{}{}_top{}_inter_{}{}.png".format(img_dir, msa_name, len(df_dca), atom, cutoff)     # interface-only
        imgname = "{}{}_top{}_FNi_{}{}.png".format(img_dir, msa_name, len(df_dca), atom, cutoff)     # FNi
    else:
        img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\distance_matrix\\"
        imgname = "{}{}_len_{}_{}{}.png".format(img_dir, msa_name, length, atom, cutoff)
    plt.savefig(imgname, dpi=900, bbox_inches='tight')
    plt.show()
    # plt.close()


def vectorize_pdb_contacts(pdb_df, dimer_length):
    """
    Converts PDB pairs into binary matrix and then flattens into 1-D array
    :param pdb_df: Dataframe object - PDB pairs
    :param dimer_length: int - Length of dimer
    :return: 1-D flat array
    """
    import pandas
    import numpy as np
    pdb_array = pdb_df.iloc[:, :2].to_numpy()  # convert column i and j to numpy array
    # Initialize contact matrices
    pdb_matrix = np.zeros((dimer_length, dimer_length))

    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in pdb_array:
        pdb_matrix[int(i)-1, int(j)-1] = 1

    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    pdb_flat_matrix = pdb_matrix.ravel()
    return pdb_flat_matrix
