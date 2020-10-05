import matplotlib.pylab as plt
from sklearn.metrics import precision_score, confusion_matrix
from distance_dca import vectorize_dca_contacts


def plot_true_positives(dca_df, atom, msa_name=None, cutoff=None):
    """
    Plots True positive rate for each top DCA prediction
    :param dca_df:
    :param msa_name:
    :param cutoff:
    :param atom:
    :return:
    """
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\true_positives_interface\\"
    from numpy.random import randint

    count = ord(msa_name[1]) * randint(500, 900) + ord(msa_name[2]) * randint(1001, 2000)
    n_pairs = len(dca_df)
    plt.figure(count, figsize=(10, 10))
    # -- DCA --
    plt.plot(range(1, n_pairs+1), dca_df["tp"][:n_pairs], label="{} Vanilla DCA interface".format(msa_name))

    plt.title("{}-cutoff = {}".format(atom, cutoff))
    plt.xlabel('ranked DCA pairs')
    plt.ylabel('True positives')
    plt.grid(axis='both')
    plt.legend(loc='best')

    imgname = "TP_inter_plmDCA_{}_top{}_{}{}.png".format(msa_name, n_pairs, atom, cutoff)
    plt.savefig(img_dir + imgname, dpi=500, bbox_inches='tight')
    # plt.show()
    plt.close()


def plot_tpr_as_function_of_pairs(n_pairs, tpr_unbinned, threshold, nSystems):
    """
    Plot TPR as function of top ranked DCA pairs

    :param n_pairs:
    :param tpr_unbinned:
    :param threshold:
    :param nSystems:
    :return:
    """
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"
    from numpy.random import randint
    plotNumber = randint(nSystems*n_pairs)
    plt.figure(plotNumber)
    plt.plot(range(n_pairs), tpr_unbinned, color='black', label='{}$\AA$'.format(threshold))
    plt.xlabel("top ranked DCA intermonomer predictions")
    plt.ylabel("TPR")
    imgname_tpr_vs_top_pairs = "TPR_as_function_top{}_pairs_{}systems_inter_plmDCA_aa{}.png".format(n_pairs, nSystems,
                                                                                                    threshold)
    plt.savefig(img_dir + imgname_tpr_vs_top_pairs, dpi=600, bbox_inches='tight')
    plt.close()


def plot_tpr_as_function_of_scores(TPR_norm, binned_scores, width, n_pairs, threshold, nSystems):
    """
    Plot TPR vs Binned FN-apc scores

    :param TPR_norm:
    :param binned_scores:
    :param width:
    :param n_pairs:
    :param threshold:
    :param nSystems:
    :return:
    """
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"
    from numpy.random import randint
    plotNumber = randint(nSystems*n_pairs)
    plt.figure(plotNumber)
    plt.bar(binned_scores, TPR_norm, width=width, edgecolor="black", color="teal", alpha=0.7)
    plt.xlabel("binned FN-apc scores")
    plt.ylabel("TPR")
    imgname_tpr_vs_binned_fn = "TPR_as_function_FN_dx{}_{}systems_inter_plmDCA_top{}_aa{}.png".format(width, nSystems,
                                                                                                      n_pairs,
                                                                                                      threshold)
    plt.savefig(img_dir + imgname_tpr_vs_binned_fn, dpi=600, bbox_inches='tight')
    # plt.show()
    plt.close()


def plot_npairs_as_function_of_scores(norm, binned_scores, width, n_pairs, threshold, nSystems):
    # Plot Normalization (eg Number of pairs) vs Binned FN-apc scores
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"
    plt.figure(2)
    plt.bar(binned_scores, norm, width=width, edgecolor="black", color="navy", alpha=0.7)
    plt.yscale(value="log")
    plt.xlabel("binned FN-apc scores")
    plt.ylabel("number of pairs")
    imgname_pairs_vs_fn = "number_of_pairs_as_function_FN_dx{}_{}systems_top{}_aa{}.png".format(width, nSystems,
                                                                                                n_pairs, threshold)
    plt.savefig(img_dir + imgname_pairs_vs_fn, dpi=600, bbox_inches='tight')
    # plt.show()
    plt.close()


def tpr_top_pairs(pdb_flat_matrix, dca_df, n_contacts, pdb_total_length):
    """
    # This uses confusion matrix method. Don't use this. There is a better way by defining a cutoff threshold.
    :param pdb_total_length:
    :param pdb_flat_matrix:
    :param dca_df:
    :param n_contacts:
    :return:
    """
    tpr_list = []
    for i in range(n_contacts):
        if i > 0:
            dca_flat_matrix = vectorize_dca_contacts(dca_df[:i], pdb_total_length)
            tpr = precision_score(pdb_flat_matrix, dca_flat_matrix, zero_division=1)
            tpr_list.append(tpr)
    return tpr_list


def confusion_matrix_list(pdb_flat_matrix, dca_df, pdb_total_length):
    """
    # This uses confusion matrix method to calculate tp. Dont use this.
    :param pdb_flat_matrix:
    :param dca_df:
    :param pdb_total_length:
    :return:
    """
    tp_list = []
    fp_list = []
    fn_list = []
    tn_list = []
    n_contacts = len(dca_df)
    for i in range(n_contacts + 1):
        if i > 0:
            dca_flat_matrix = vectorize_dca_contacts(dca_df[:i], pdb_total_length)
            tn, fp, fn, tp = confusion_matrix(pdb_flat_matrix, dca_flat_matrix, normalize="pred").ravel()
            tp_list.append(tp)
            fp_list.append(fp)
            fn_list.append(fn)
            tn_list.append(tn)
    return tp_list, fp_list, fn_list, tn_list


