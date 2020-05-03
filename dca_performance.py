from sklearn.metrics import precision_score
import numpy as np
import logging


def run_analysis(msa_name, n_pairs, cutoff, results_dir, msa_dir, pdb_dir, img_dir, sasa_calc=False):
    from plot_cm import plot_cm, draw_dca_tcl
    from pdb import pdb_map, map_msa_to_pdb, cat_seq, cm_make, read_length_file, apply_map
    # P2: There are a lot of redundancies that can be called once, saved to file and then read file
    #  instead of running every time.
    import pandas as pd
    gremlin_dir = "gremlin\\"
    #gremlin_file = "{}GRM_{}.txt_inter".format(gremlin_dir, msa_name)
    #print("Gremlin {}".format(gremlin_file))
    #df_gremlin = pd.read_csv(gremlin_file, delim_whitespace=True, names=['i', 'j', 'score'], usecols=(1, 2, 3), skiprows=2)

    try:
        vanilla_dca = '{}fn_{}_plmdca_rt2.txt'.format(results_dir, msa_name)
        # df_vanilla = pd.read_csv(vanilla_dca, delimiter=',', names=['i', 'j', 'score'],
        #                         usecols=(0, 1, 2))
        # df_vanilla = df_vanilla.sort_values(ascending=False, by=['score'])
        dca_score_matrix = '{}FNi_{}.txt'.format(results_dir, msa_name)

        # Cat PDB chain seqs together and find msa to pdb alignment
        print("\tcatting seqs...")
        full_seq, msa_seq, chains = cat_seq(msa_name, msa_dir, pdb_dir, get_length=None)
        print("\tcreating map...")
        dca_indices, pdb_indices, map_to_dca, map_to_pdb = map_msa_to_pdb(full_seq, msa_seq)
        total_length = len(pdb_indices)

        # PDB contact map
        pdbid = msa_name[:4]
        pdbfile = "{}{}.cif".format(pdb_dir, pdbid)
        df_pdb, df_mon, df_inter = pdb_map(pdbfile, chains, cutoff)
        pdb_df_list = [df_mon, df_inter]

        # make contact map and apply backmapping
        # df_umap = cm_make(dca_score_matrix)
        print("\tapplying map to {}...".format(dca_score_matrix))
        df_map = cm_make(dca_score_matrix, map_to_pdb, dca_indices[0])
        # vanilla_array = apply_map(df_vanilla.to_numpy(), map_to_pdb, dca_indices[0])
        # df_map_v = pd.DataFrame(vanilla_array, columns=['i', 'j', 'score'])
        #gremlin_array = apply_map(df_gremlin.to_numpy(), map_to_pdb, dca_indices[0])
        #df_map_g = pd.DataFrame(gremlin_array, columns=['i', 'j', 'score'])
        #df_map_g = df_map_g.sort_values(ascending=False, by=['score'])

        # Get length of first chain
        length_a = read_length_file(msa_name)

        if sasa_calc:
            # run sasa program
            from sasa_calc import sasa_program
            df_sasa_dca = sasa_program(df_map, n_pairs, pdbid, chains, cutoff, pdb_dir)

        # Plot for a range of number of DCA pairs
        step = 10
        for i in range(10, n_pairs + step, step):
            print("Plotting top {}...".format(i))
            contact_map_file = "{}cm_FNi_{}.txt".format(results_dir, msa_name)
            draw_dca_tcl(contact_map_file, i, length_a, chains)
            plot_cm(pdb_df_list, df_map, i, cutoff, length_a, total_length, img_dir=img_dir, gremlin_pairs=None,
                    vanillafile=None, title="DCAi", msa_name=msa_name)

        # Plots TPR for top DCA predictions (compared to PDB interface pairs)
        print("Plotting TPR...")
        plot_performance(df_inter, df_map, n_pairs, total_length, msa_name, gremlin_pairs=None, sasa=None,
                         cutoff=cutoff, img_dir=img_dir)
        print("Finished analyzing {}".format(msa_name))

    except OSError:
        logging.debug("\t\t{}.fas has not been analysed.".format(msa_name))
        print("\t\t{}.fas NOT analysed yet.".format(msa_name))


def vectorize_dca_contacts(df_dca, dimer_length):
    """
    Converts DCA pairs into binary matrix and then flattens into 1-D array
    :param df_dca: Dataframe object - DCA pairs NOTE: Should already be sliced to desired length
    :param dimer_length: int - Length of dimer
    :return: 1-D flat array
    """
    import pandas
    # Initialize contact matrices
    dca_matrix = np.zeros((dimer_length, dimer_length))

    dca_array = df_dca.iloc[:, :2].to_numpy()  # convert column i and j to numpy array
    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in dca_array:
        dca_matrix[int(i), int(j)] = 1
    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    dca_flat_matrix = dca_matrix.ravel()
    return dca_flat_matrix


def vectorize_pdb_contacts(pdb_df, dimer_length):
    """
    Converts PDB pairs into binary matrix and then flattens into 1-D array
    :param pdb_df: Dataframe object - PDB pairs
    :param dimer_length: int - Length of dimer
    :return: 1-D flat array
    """
    from pdb import pdb_map
    import pandas
    pdb_array = pdb_df.iloc[:, :2].to_numpy()  # convert column i and j to numpy array
    # Initialize contact matrices
    pdb_matrix = np.zeros((dimer_length, dimer_length))

    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in pdb_array:
        pdb_matrix[int(i), int(j)] = 1

    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    pdb_flat_matrix = pdb_matrix.ravel()
    return pdb_flat_matrix


def tpr_top_pairs(pdb_flat_matrix, dca_df, n_contacts, dimer_length):
    """
    :param pdb_flat_matrix:
    :param dca_df:
    :param n_contacts:
    :param dimer_length:
    :return:
    """
    tpr_list = []
    for i in range(n_contacts - 1):
        dca_flat_matrix = vectorize_dca_contacts(dca_df[:i + 1], dimer_length)
        tpr = precision_score(pdb_flat_matrix, dca_flat_matrix, zero_division=1)
        tpr_list.append(tpr)
    return tpr_list


def plot_performance(pdb_df, dca_df, n_contacts, dimer_length, msa_name=None, gremlin_pairs=None, sasa=None,
                     cutoff=None, atom="ca", img_dir=None):
    """
    Plots True positive rate for each top DCA prediction
    :param sasa:
    :param gremlin_pairs:
    :param img_dir:
    :param pdb_df:
    :param dca_df:
    :param n_contacts:
    :param dimer_length:
    :param msa_name:
    :param cutoff:
    :param atom:
    :return:
    """
    import matplotlib.pylab as plt

    # make flat array of pdb and calculate tpr for top pairs
    pdb_flat_array = vectorize_pdb_contacts(pdb_df, dimer_length)
    tpr_list = tpr_top_pairs(pdb_flat_array, dca_df, n_contacts, dimer_length)

    plt.figure(0, figsize=(10, 10))
    # -- DCA --
    plt.plot(range(1, n_contacts), tpr_list, label=msa_name)

    # -- Gremlin --
    #if not gremlin_pairs.empty:
    #    tpr_list_g = tpr_top_pairs(pdb_flat_array, gremlin_pairs, n_contacts, dimer_length)
    #    plt.plot(range(1, n_contacts), tpr_list_g, label="gremlin_{}".format(msa_name), linestyle="dashed")

    plt.hlines(.5, 1, n_contacts, linestyles='dashed', alpha=0.6)
    plt.vlines(20, 0, max(tpr_list), linestyles='dashed', alpha=0.6)

    plt.title("{}-cutoff = {}".format(atom, cutoff))
    plt.xlabel('number of dca predictions')
    plt.ylabel('tpr')
    plt.grid(axis='both')
    plt.legend(loc='top-right')

    if img_dir:
        imgname = "tpr_{}_top{}_{}{}.png".format(msa_name, n_contacts, atom, cutoff)
        plt.savefig(img_dir + imgname, dpi=500)
