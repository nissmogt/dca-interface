import matplotlib.pylab as plt
from sklearn.metrics import precision_score, confusion_matrix
from distance_dca import vectorize_dca_contacts
img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\DEC_2020\\"


def plot_tpr_as_function_of_pairs(n_pairs, true_positives_array, threshold, nSystems, fni=False):
    """
    Plot TPR as function of top ranked DCA pairs

    :param fni:
    :param n_pairs:
    :param true_positives_array:
    :param threshold:
    :param nSystems:
    :return:
    """
    from numpy.random import randint
    plotNumber = randint(nSystems * n_pairs)
    plt.figure(plotNumber)
    plt.plot(range(n_pairs), true_positives_array, color='black', label='{}$\AA$'.format(threshold))
    plt.xlabel("top ranked DCA intermonomer predictions")
    plt.ylabel("TPR")
    if fni:
        imgname_tpr_vs_top_pairs = "TPR_as_function_top{}_pairs_{}systems_FNi_inter_plmDCA_aa{}.png".format(n_pairs,
                                                                                                            nSystems,
                                                                                                            threshold)
    else:
        imgname_tpr_vs_top_pairs = "TPR_as_function_top{}_interface_pairs_{}systems_FN_plmDCA_aa{}.png".format(n_pairs,
                                                                                                               nSystems,
                                                                                                               threshold)
    plt.savefig(img_dir + imgname_tpr_vs_top_pairs, dpi=600, bbox_inches='tight')
    plt.close()


def compare_two_tpr(n_pairs, tpr_unbinned_1, tpr_unbinned_2, threshold, nSystems):
    """
    Plot TPR as function of top ranked DCA pairs

    :param tpr_unbinned_2:
    :param n_pairs:
    :param tpr_unbinned_1:
    :param threshold:
    :param nSystems:
    :return:
    """
    from numpy.random import randint
    plotNumber = randint(nSystems * n_pairs)
    plt.figure(plotNumber)
    plt.plot(range(n_pairs), tpr_unbinned_1, color='black', label='FNi {}$\AA$'.format(threshold))
    plt.plot(range(n_pairs), tpr_unbinned_2, color='xkcd:red', label='FN {}$\AA$'.format(threshold))
    plt.xlabel("top ranked DCA intermonomer predictions")
    plt.ylabel("fraction of dca pairs with cutoff <= {}$\AA$".format(threshold))
    plt.legend(loc="best")
    imgname_tpr_vs_top_pairs = "compare_TP FNi_FN_top{}_{}systems_aa{}.png".format(n_pairs, nSystems, threshold)
    plt.savefig(img_dir + imgname_tpr_vs_top_pairs, dpi=600, bbox_inches='tight')
    plt.close()


def plot_tp_as_function_of_scores(TPR_norm, binned_scores, width, n_pairs, threshold, nSystems, fni=False):
    """
    Plot TPR vs Binned FN-apc scores

    :param TPR_norm:
    :param binned_scores:
    :param width:
    :param n_pairs:
    :param threshold:
    :param nSystems:
    :param fni:
    :return:
    """
    from numpy.random import randint
    plotNumber = randint(nSystems * n_pairs)
    plt.figure(plotNumber)
    plt.bar(binned_scores, TPR_norm, width=width, align='edge', edgecolor="black", color="teal", alpha=0.7)
    plt.ylabel("True Positives".format(threshold))

    if fni:
        plt.xlabel("FNi score")
        imgname_tp_vs_binned_fn = "FNi_TP_as_function_dx{}_{}systems_inter_plmDCA_top{}_aa{}.png".format(width,
                                                                                                           nSystems,
                                                                                                           n_pairs,
                                                                                                           threshold)
    else:
        plt.xlabel("FN scores")
        imgname_tp_vs_binned_fn = "FN_TP_as_function_all_dx{}_{}systems_plmDCA_top{}_aa{}.png".format(width,
                                                                                                              nSystems,
                                                                                                              n_pairs,
                                                                                                              threshold)
    plt.savefig(img_dir + imgname_tp_vs_binned_fn, dpi=600, bbox_inches='tight')
    # plt.xlim(0, 1.4)
    # plt.ylim(0, 1)
    plt.show()
    plt.close()


def plot_fp_as_function_of_scores(TPR_norm, binned_scores, width, n_pairs, threshold, nSystems, fni=False):
    """
    Plot TPR vs Binned FN-apc scores

    :param TPR_norm:
    :param binned_scores:
    :param width:
    :param n_pairs:
    :param threshold:
    :param nSystems:
    :param fni:
    :return:
    """
    from numpy.random import randint
    plotNumber = randint(nSystems * n_pairs)
    plt.figure(plotNumber)
    plt.bar(binned_scores, TPR_norm, width=width, align='edge', edgecolor="black", color="purple", alpha=0.7)
    plt.ylabel("False Positives".format(threshold))

    if fni:
        plt.xlabel("FNi scores")
        imgname_fp_vs_binned_fn = "FNi_FP_as_function_dx{}_{}systems_inter_plmDCA_top{}_aa{}.png".format(width,
                                                                                                           nSystems,
                                                                                                           n_pairs,
                                                                                                           threshold)
    else:
        plt.xlabel("FN scores")
        imgname_fp_vs_binned_fn = "FN_FP_as_function_interface_dx{}_{}systems_plmDCA_top{}_aa{}.png".format(width,
                                                                                                              nSystems,
                                                                                                              n_pairs,
                                                                                                              threshold)
    # plt.savefig(img_dir + imgname_fp_vs_binned_fn, dpi=600, bbox_inches='tight')
    plt.xlim(0, 1.4)
    # plt.ylim(0, 1)
    plt.show()
    # plt.close()


def plot_npairs_as_function_of_scores(norm, binned_scores, width, n_pairs, threshold, nSystems, fni=False):
    # Plot Normalization (eg Number of pairs) vs Binned FN-apc scores
    plt.figure(2)
    plt.bar(binned_scores, norm, width=width, align='edge', edgecolor="black", color="navy", alpha=0.7, label="total contacts")
    # plt.yscale(value="log")
    plt.ylabel("P(TP)")
    if fni:
        plt.xlabel("binned FNi scores")
        imgname_pairs_vs_fn = "FNi_pairs_as_function_dx{}_{}systems_top{}_aa{}.png".format(width, nSystems,
                                                                                           n_pairs, threshold)
    else:
        plt.xlabel("binned FN-apc scores")
        imgname_pairs_vs_fn = "FN_pairs_as_function_interface_dx{}_{}systems_top{}_aa{}.png".format(width, nSystems,
                                                                                                    n_pairs, threshold)

    plt.legend(loc='best')
    # plt.savefig(img_dir + imgname_pairs_vs_fn, dpi=600, bbox_inches='tight')
    plt.show()
    # plt.close()


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
