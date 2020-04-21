import numpy as np
from sklearn.metrics import confusion_matrix, precision_score, recall_score, roc_curve
from sklearn.metrics import precision_recall_curve, average_precision_score


def vectorize_dca_contacts(dca_array, dimer_length, n_contacts):
    """

    :param dca_array:
    :param dimer_length:
    :param n_contacts:
    :return:
    """
    # Initialize contact matrices
    dca_matrix = np.zeros((dimer_length, dimer_length))

    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in dca_array[:n_contacts]:
        print('dca_array:', i, j)
        dca_matrix[int(i), int(j)] = 1
    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    dca_flat_matrix = dca_matrix.ravel()
    return dca_flat_matrix


def vectorize_pdb_contacts(pdb_file, dimer_length, chain, cutoff):
    """

    :param pdb_file:
    :param dimer_length:
    :param chain:
    :param cutoff:
    :return:
    """
    from pdb_contacts import pdb_map
    df_pdb, df_mon, df_inter = pdb_map(pdb_file, chain, cutoff)
    pdb_array = df_inter.iloc[:, :2].to_numpy()
    # Initialize contact matrices
    pdb_matrix = np.zeros((dimer_length, dimer_length))

    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in pdb_array:
        pdb_matrix[int(i), int(j)] = 1

    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    pdb_flat_matrix = pdb_matrix.ravel()
    return pdb_flat_matrix


def loop_tpr(pdb_flat_matrix, dca_array, n_contacts, dimer_length, chain, cutoff):
    """
    :param pdb_flat_matrix:
    :param dca_array:
    :param n_contacts:
    :param dimer_length:
    :param chain:
    :param cutoff:
    :return:
    """
    tpr_list = []
    for i in range(n_contacts):
        dca_flat_matrix = vectorize_dca_contacts(dca_array, dimer_length, i)
        tpr = precision_score(pdb_flat_matrix, dca_flat_matrix, zero_division=1)
        tpr_list.append(tpr)
    return tpr_list


def plot_performance(pdb_file, dca_array, n_contacts, dimer_length, chain, cutoff):
    import matplotlib.pylab as plt
    # grem_array = np.loadtxt(gremlinfile, usecols=(0, 1), dtype=int)
    # plot_cm(pdbfile, dcafile, chain, 10, 15)

    tpr_list = loop_tpr(pdb_file, dca_array, n_contacts, dimer_length, chain, cutoff)
    plt.figure(0)
    plt.plot(range(len(tpr_list)), tpr_list)
    plt.xlabel('number of dca predictions')
    plt.ylabel('tpr')
    plt.grid(axis='both')
