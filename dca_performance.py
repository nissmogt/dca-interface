import numpy as np
from sklearn.metrics import confusion_matrix, precision_score, recall_score, roc_curve
from sklearn.metrics import precision_recall_curve, average_precision_score


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


def plot_performance(pdb_df, dca_df, n_contacts, dimer_length, msa_name=None, cutoff=None, atom="ca"):
    """
    Plots True positive rate for each top DCA prediction
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

    pdb_flat_array = vectorize_pdb_contacts(pdb_df, dimer_length)
    tpr_list = tpr_top_pairs(pdb_flat_array, dca_df, n_contacts, dimer_length)
    plt.figure(0)
    plt.plot(range(len(tpr_list)), tpr_list)
    plt.xlabel('number of dca predictions')
    plt.ylabel('tpr')
    plt.grid(axis='both')
    if msa_name:
        img_dir = "C:\\Users\\kmehr\\Google Drive\\work\\images\\2020\\APR_2020\\Week_4\\"
        imgname = "tpr_{}_top{}_{}{}.png".format(msa_name, n_contacts, atom, cutoff)
        plt.savefig(img_dir + imgname, dpi=500)
