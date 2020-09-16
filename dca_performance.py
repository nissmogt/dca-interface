from sklearn.metrics import precision_score, confusion_matrix
import numpy as np
import logging


def run_analysis(msa_file, n_pairs, cutoff, results_dir, msa_dir, pdb_dir, sasa_calc=False, count=0):
    import os
    import PreProcess as pp
    import pandas as pd
    from distance_pdb import distance_matrix
    from get_residues import get_residues
    from draw_contacts import draw_dca_tcl
    from load_dca import load_dca
    from plot_cm import plot_cm
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"

    msa_name = os.path.basename(msa_file).strip(".fas")
    print("\n\t-- LOADING: {} --".format(msa_name))

    # -- VANILLA DCA --
    df_vanilla = load_dca(msa_name, results_dir)

    # -- PDB --
    atom = 'ca'
    # pdbid = msa_name[:4].lower()
    # pdbfile = "{}{}.cif".format(pdb_dir, pdbid)
    p = pp.Preprocess(msa_file, cutoff)
    # print("(pdb filename)\t{}".format(os.path.basename(p.pdbfile)))
    # df_mon, df_inter = p.read_distance_matrix_file()
    df_pdb, df_mon, df_inter, cl = distance_matrix(msa_name, cutoff)
    pdb_df_list = [df_mon, df_inter]

    # -- GET MAP FROM MSA TO PDB --
    pdbseq_1, pdbseq_2 = p.get_residues(seq=True)
    pdbseq = pdbseq_1 + pdbseq_2
    length_a = p.chain_lengths[0]
    # msaseq1, msaseq2 = p.msa_template(split=True, len_a=length_a+1)
    msaseq = p.msa_template()

    a, dca_indices, pdb_indices, map_to_dca, map_to_pdb = p.get_msa_to_pdb_dictionary(pdbseq, msaseq)
    # total_length = len(pdb_indices)
    total_length = p.chain_lengths[0] + p.chain_lengths[1]
    print("\nChain {} len: {}\tTotal Length: {}".format(p.chain_ids[0], length_a, total_length))
    print("(map dictionary) {}".format(map_to_pdb))
    # mapped_dca_file = "{}mapped_cm_FNi_{}.txt".format(results_dir, msa_name)
    # print("(dca_map_filename)\t{}".format(mapped_dca_file))
    # dca_score_matrix = "{}FNi_apc_{}.txt".format(results_dir, msa_name)
    # dca_score_matrix = "GRM_{}_APC_prior.txt".format(msa_name)
    # print("(dca score matrix) {}".format(dca_score_matrix))
    # df_umap = p.cm_make(dca_score_matrix)
    # df_map = p.cm_make(dca_score_matrix, backmap=True, map_dictionary=map_to_pdb)
    # df_map = df_map.loc[(df_map["i"] <= length_a) & (df_map["j"] >= length_a)]
    # VANILLA plmDCA
    vanilla_map = p.apply_map(df_vanilla.to_numpy(), map_to_pdb)
    df_map_v = pd.DataFrame(vanilla_map, columns=['i', 'j', 'score'])
    # df_map_v = df_map.loc[(df_map_v["i"] <= length_a) & (df_map_v["j"] >= length_a)]
    # np.savetxt("{}mapped_cm_plmDCA_apc_{}.txt".format(results_dir, msa_name), df_map_v, fmt="%d\t%d\t%f")

    # -- GREMLIN --
    # gremlin_dir = "gremlin\\"
    # gremlin_file = "{}GRM_{}.txt".format(gremlin_dir, msa_name)
    # gremlin_file = "GRM_{}_APC.txt".format(msa_name)
    # print("\n(gremlin filename)\tGremlin GRM_{}.txt".format(msa_name))
    # df_gremlin = pd.read_csv(gremlin_file, delim_whitespace=True, names=['i', 'j', 'score', 'i_score'],
    #                          usecols=(0, 1, 5, 7), skiprows=1)
    # df_gremlin = pd.read_csv(gremlin_file, delim_whitespace=True, names=['i', 'j', 'score'],
    #                          usecols=(0, 1, 5), skiprows=1)
    # df_map_g = p.cm_make(gremlin_file, backmap=True, map_dictionary=map_to_pdb)
    # df_map_g = df_map.loc[(df_map_g["i"] <= length_a) & (df_map_g["j"] >= length_a)]

    # df_gremlin = df_gremlin.dropna() - 1  # subtracting one starts i at 0
    # df_gremlin['score'] = df_gremlin['score'] + 1
    # df_gremlin['i_score'] = df_gremlin['score'] + 1
    # gremlin_array = p.apply_map_g(df_gremlin.to_numpy(), map_to_pdb)
    # df_map_g = pd.DataFrame(gremlin_array, columns=['i', 'j', 'score', 'i_score'])

    # df_map_g = pd.DataFrame(gremlin_array, columns=['i', 'j', 'score'])

    # df_map_g = df_map_g.sort_values(ascending=False, by=['i_score'])

    # df_map_g = df_map_g.sort_values(ascending=False, by=['score'])
    # df_map_g = df_map_g.loc[(df_map_g["i"] <= length_a) & (df_map_g["j"] >= length_a)]
    # np.savetxt("{}map_{}".format(gremlin_dir, os.path.basename(gremlin_file)), df_map_g, fmt='%d\t%d\t%f\t%f')

    pdb_flat_array = p.vectorize_pdb_contacts(df_inter, total_length)
    df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm

    # if sasa_calc:
    # run sasa program
    # from sasa_calc import sasa_program
    # df_sasa_dca = sasa_program(df_map, n_pairs, pdbid, chains, cutoff, pdb_dir)

    # -- PLOTS --
    step = 10
    # for i in range(10, n_pairs + step, step):
    # top_folder = "top{}\\".format(i)
    i = n_pairs
    draw_dca_tcl(df_map.to_numpy(), i, length_a+1, p.chain_ids, msa_name)
    print("Plotting top {}...".format(i))
    plot_cm(pdb_df_list, df_map_v, i, cutoff, length_a, total_length, df_dca_umap=df_empty, gremlin_pairs=df_empty,
            vanillafile=df_map_v, title="plmDCA", atom=atom, msa_name=msa_name)

    # Plots TP for top DCA predictions (compared to PDB interface pairs)
    print("Plotting True positives...")
    # tp_list = tp_top_pairs(pdb_flat_array, df_map, n_pairs, total_length)
    # np.savetxt("{}TP_{}.txt".format(results_dir, msa_name), tp_list, fmt='%d')
    # plot_tp(pdb_flat_array, df_map, n_pairs, total_length, msa_name, gremlin_pairs=df_map_g, sasa=None,
    #         cutoff=cutoff, atom=atom, img_dir=img_dir, count=count)
    # Plots TPR for top DCA predictions (compared to PDB interface pairs)
    print("(plot)\tPlotting TPR...")
    # plot_performance(pdb_flat_array, df_map, n_pairs, total_length, msa_name, gremlin_pairs=df_map_g, sasa=None,
    #                  cutoff=cutoff, atom=atom, img_dir=img_dir, count=count)

    print("(PROGRAM END)\tFinished analyzing {}".format(msa_name))

    return df_map


def map_cm(msa_file, cutoff, results_dir, pdb_dir):
    from pdb import cm_make
    import os
    import PreProcess as pp
    import pandas as pd
    msa_name = os.path.basename(msa_file).strip(".fas")

    try:
        dca_score_matrix = '{}FNi_{}.txt'.format(results_dir, msa_name)  # PDB contact map
        pdbid = msa_name[:4]
        pdbfile = "{}{}.cif".format(pdb_dir, pdbid)
        p = pp.Preprocess(pdbfile, msa_file, cutoff)
        chains = p.chain_ids  # Cat PDB chain seqs together and find msa to pdb alignment
        print("\tcatting seqs...")
        # full_seq, msa_seq, chains = cat_seq(msa_name, msa_dir, pdb_dir, get_length=None)
        pdbseq = p.combine_chain_sequence()
        msaseq = p.msa_template()
        print("\treading distance matrix...")
        dca_indices, pdb_indices, map_to_dca, map_to_pdb = p.get_msa_to_pdb_dictionary(pdbseq, msaseq)
        # dca_indices, pdb_indices, map_to_dca, map_to_pdb = map_msa_to_pdb(full_seq, msa_seq)
        total_length = len(pdb_indices)  # make contact map and apply backmapping
        df_umap = cm_make(dca_score_matrix)
        # print("\tapplying map to {}...".format(dca_score_matrix))
        # df_map = cm_make(dca_score_matrix, map_to_pdb, dca_indices[0])
    except OSError:
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


def tpr_top_pairs(pdb_flat_matrix, dca_df, n_contacts, dimer_length):
    """
    :param pdb_flat_matrix:
    :param dca_df:
    :param n_contacts:
    :param dimer_length:
    :return:
    """
    tpr_list = []
    for i in range(n_contacts):
        if i > 0:
            dca_flat_matrix = vectorize_dca_contacts(dca_df[:i], dimer_length)
            tpr = precision_score(pdb_flat_matrix, dca_flat_matrix, zero_division=1)
            tpr_list.append(tpr)
    return tpr_list


def tp_top_pairs(pdb_flat_matrix, dca_df, n_contacts, dimer_length):
    """
    :param pdb_flat_matrix:
    :param dca_df:
    :param n_contacts:
    :param dimer_length:
    :return:
    """
    tp_list = []
    for i in range(n_contacts):
        if i > 0:
            dca_flat_matrix = vectorize_dca_contacts(dca_df[:i], dimer_length)
            tn, fp, fn, tp = confusion_matrix(pdb_flat_matrix, dca_flat_matrix).ravel()
            tp_list.append(tp)
    return tp_list


def plot_performance(pdb_flat_array, dca_df, n_contacts, dimer_length, msa_name=None, gremlin_pairs=None, sasa=None,
                     cutoff=None, atom="ca", img_dir=None, count=0):
    """
    Plots True positive rate for each top DCA prediction
    :param count:
    :param sasa:
    :param gremlin_pairs:
    :param img_dir:
    :param pdb_flat_array:
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
    # pdb_flat_array = vectorize_pdb_contacts(pdb_flat_array, dimer_length)
    tpr_list = tpr_top_pairs(pdb_flat_array, dca_df, n_contacts, dimer_length)

    plt.figure(count + 100000, figsize=(10, 10))
    # -- DCA --
    plt.plot(range(1, n_contacts), tpr_list, label=msa_name)

    # -- Gremlin --
    if not gremlin_pairs.empty:
        tpr_list_g = tpr_top_pairs(pdb_flat_array, gremlin_pairs, n_contacts, dimer_length)
        plt.plot(range(1, n_contacts), tpr_list_g, label="gremlin_{}".format(msa_name), linestyle="dashed")

    plt.title("{}-cutoff = {}".format(atom, cutoff))
    plt.xlabel('number of dca predictions')
    plt.ylabel('tpr')
    plt.grid(axis='both')
    plt.legend(loc='best')

    if img_dir:
        imgname = "tpr_apc_vs_GRM_{}_top{}_{}{}.png".format(msa_name, n_contacts, atom, cutoff)
        plt.savefig(img_dir + imgname, dpi=500, bbox_inches='tight')
        plt.close()


def plot_tp(pdb_flat_array, dca_df, n_contacts, dimer_length, msa_name=None, gremlin_pairs=None, sasa=None, cutoff=None,
            atom="ca", img_dir=None, count=0):
    """
    Plots True positive rate for each top DCA prediction
    :param count:
    :param sasa:
    :param gremlin_pairs:
    :param img_dir:
    :param pdb_flat_array:
    :param dca_df:
    :param n_contacts:
    :param dimer_length:
    :param msa_name:
    :param cutoff:
    :param atom:
    :return:
    """
    import matplotlib.pylab as plt
    import random

    # make flat array of pdb and calculate tpr for top pairs
    # pdb_flat_array = vectorize_pdb_contacts(pdb_flat_array, dimer_length)
    tp_list = tp_top_pairs(pdb_flat_array, dca_df, n_contacts, dimer_length)

    plt.figure(count + 1000, figsize=(10, 10))
    # y = x
    plt.plot(range(n_contacts), range(n_contacts), color='black')
    # -- DCA --
    plt.plot(range(1, n_contacts), tp_list, label=msa_name + "DCAi")

    # -- Gremlin --
    if not gremlin_pairs.empty:
        tp_list_g = tp_top_pairs(pdb_flat_array, gremlin_pairs, n_contacts, dimer_length)
        plt.plot(range(1, n_contacts), tp_list_g, label="gremlin_{}".format(msa_name), linestyle="dashed")

    plt.title("{}-cutoff = {}".format(atom, cutoff))
    plt.xlabel('number of dca predictions')
    plt.ylabel('True positives')
    plt.grid(axis='both')
    plt.legend(loc='best')

    if img_dir:
        imgname = "TP_plmDCA_GRM_{}_top{}_{}{}.png".format(msa_name, n_contacts, atom, cutoff)
        plt.savefig(img_dir + imgname, dpi=500, bbox_inches='tight')
        plt.close()




