import numpy as np


def cm_make(msa_name, score_matrix, apc=True):
    """
    Makes a contact map from a DCA score matrix file. Also maps dca indices to pdb.
    Rank by scores.
    :param apc:
    :param msa_name:
    :param score_matrix: Frobenius norm matrix
    :return: Three-column dataframe composed of pair i, j, and fn score
    """
    import os
    import pandas as pd
    fname = "(cm_make)"
    results = "scrambled_results\\"
    L = score_matrix.shape[0]
    dca_scores = []
    for i in range(L - 1):
        for j in range(i + 1, L):
            dca_scores.append([i+1, j+1, score_matrix[i, j]])
    df_dca = pd.DataFrame(np.array(dca_scores), columns=["i", "j", "score"])
    df_dca = df_dca.sort_values(ascending=False, by=["score"])
    df_dca = df_dca[abs(df_dca["i"] - df_dca["j"]) > 5]    # ensure sequence separation of 5
    if not apc:
        np.savetxt('{}FNi_{}.txt'.format(results, msa_name), df_dca, delimiter='\t')
    else:
        np.savetxt('{}FNi_apc_{}.txt'.format(results, msa_name), df_dca, delimiter='\t')
    return df_dca


def calculate_fni_score(msa_name, j_matrix, j_matrix_null):
    # Function that calculates Frobenius Norm score of coupling matrices.
    q = j_matrix.shape[0]
    L = j_matrix_null.shape[1]
    fn_1 = np.zeros((L, L))
    fn_2 = np.zeros((L, L))
    fni_scores = np.zeros((L, L))
    print("Size %d" % L)

    for i in range(L - 1):
        for j in range(i + 1, L):
            fn_1[i, j] = np.linalg.norm((j_matrix[i, j]), 'fro')
            fn_2[i, j] = np.linalg.norm((j_matrix_null[i, j]), 'fro')
            # Take the difference between paired and scrambled FN
            fni_scores[i, j] = (fn_1[i, j] - fn_2[i, j])
            fni_scores[j, i] = fni_scores[i, j]

    L = fni_scores.shape[0]
    if L != fni_scores.shape[1]:
        raise ValueError("Input matrix is not symmetric: {}".format(fni_scores.shape))
    col_mean = np.mean(fni_scores, axis=0) * L / (L - 1)
    row_mean = np.mean(fni_scores, axis=1) * L / (L - 1)
    matrix_mean = np.mean(fni_scores) * L / (L - 1)

    # APC correction
    # corrnorms = fni_scores - np.outer(col_mean, row_mean) / matrix_mean
    apc = np.dot(col_mean.reshape(L, 1), col_mean.reshape(1, L)) / matrix_mean
    corrnorms = fni_scores - apc
    corrnorms[np.diag_indices(L)] = 0
    df_scores_apc = cm_make(msa_name, corrnorms)
    df_scores = cm_make(msa_name, fni_scores, apc=False)
    return df_scores, df_scores_apc


def read_matlab_matrix(filein):
    # Function that reads in the couplings and fields Matlab matrix file.
    import h5py
    mat = {}
    f = h5py.File(filein, 'r')
    for k, v in f.items():
        mat[k] = np.array(v)
    h = mat['h']
    J = mat['J']
    q = mat['h'].shape[0]
    N = mat['h'].shape[1]
    couplings = J.T
    fields = h.T
    f.close()
    return fields, couplings


def average_jmatrix(msa_name, nReplicates):
    # Function that averages over coupling matrices.
    # initialize coupling sum to zero
    sum_couplings = 0.0
    resultDir = "scrambled_results\\"
    for i in range(nReplicates):
        dca_matrix = "{}matrix_ising_{}_rep{}_scrambled.fas.mat".format(resultDir, msa_name, i)
        fields, couplings = read_matlab_matrix(dca_matrix)
        sum_couplings += couplings
    average_couplings = sum_couplings / nReplicates
    return average_couplings


# msa = "1GL2_A_1GL2_D"
# vanilla_dca_matrix = "results\\matrix_files\\matrix_ising_{}.fas.mat".format(msa)
# avg_J = average_jmatrix(msa, 5)
# h_paired, J_paired = read_matlab_matrix(vanilla_dca_matrix)
# df_fni, df_fni_apc = calculate_fni_score(msa, J_paired, avg_J)
