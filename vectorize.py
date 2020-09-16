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
        pdb_matrix[int(i), int(j)] = 1

    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    pdb_flat_matrix = pdb_matrix.ravel()
    return pdb_flat_matrix


def vectorize_dca_contacts(df_dca, dimer_length):
    """
    Converts DCA pairs into binary matrix and then flattens into 1-D array
    :param df_dca: Dataframe object - DCA pairs NOTE: Should already be sliced to desired length
    :param dimer_length: int - Length of dimer
    :return: 1-D flat array
    """
    import pandas
    # Initialize contact matrices
    import numpy as np
    dca_matrix = np.zeros((dimer_length, dimer_length))

    dca_array = df_dca.iloc[:, :2].to_numpy()  # convert column i and j to numpy array
    # Vectorize contact pairs into binary array of shape (L,L)
    for i, j in dca_array:
        dca_matrix[int(i), int(j)] = 1
    # Flatten binary array to shape (L*L) for use in confusion matrix
    # Note: Index of pair(i,j) = L*i + j
    dca_flat_matrix = dca_matrix.ravel()
    return dca_flat_matrix

