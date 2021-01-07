import os
import pandas as pd
import numpy as np
from dca_functions import cm_make


def average_fn(msa_name, nReps):
    fn_matrix = 0.0
    matrixDir = "/scratch/kmm5/scrambled_results/"
    # matrixDir = "scrambled_results\\"
    for i in range(nReps):
        dca_matrix = "{}matrix_ising_{}_rep{}_scrambled.fas.mat".format(matrixDir, msa_name, i)
        j_matrix, j_zero_gauge = load_matrix_matlab(dca_matrix)
        fn_score_matrix = calculate_fn_score(j_matrix)
        fn_matrix += fn_score_matrix
    average_fn_matrix = fn_matrix / nReps
    return average_fn_matrix


# def calculate_fni_score(msa_name, j_matrix, j_matrix_null):
# Function that calculates Frobenius Norm score of coupling matrices.
# q = j_matrix.shape[2]
# L = j_matrix_null.shape[0]
# fn_1 = np.zeros((L, L))
# fn_2 = np.zeros((L, L))
# fni_scores = np.zeros((L, L))
# print("Size %d" % L)
#
# for i in range(L - 1):
#     for j in range(i + 1, L):
#         fn_1[i, j] = np.linalg.norm((j_matrix[i, j]), 'fro')
#         fn_2[i, j] = np.linalg.norm((j_matrix_null[i, j]), 'fro')
#         Take the difference between paired and scrambled FN
# fni_scores[i, j] = (fn_1[i, j] - fn_2[i, j])
# fni_scores[j, i] = fni_scores[i, j]
#
# L = fni_scores.shape[0]
# if L != fni_scores.shape[1]:
#     raise ValueError("Input matrix is not symmetric: {}".format(fni_scores.shape))
# col_mean = np.mean(fni_scores, axis=0) * L / (L - 1)
# row_mean = np.mean(fni_scores, axis=1) * L / (L - 1)
# matrix_mean = np.mean(fni_scores) * L / (L - 1)
#
# APC correction
# corrnorms = fni_scores - np.outer(col_mean, row_mean) / matrix_mean
# apc = np.dot(col_mean.reshape(L, 1), col_mean.reshape(1, L)) / matrix_mean
# corrnorms = fni_scores - apc
# corrnorms[np.diag_indices(L)] = 0
# return fni_scores, corrnorms


def zero_sum_gauge(J_ij, inplace=False):
    """
    Copied from evcoupling.couplings model.py
    Transform coupling matrix into zero-sum gauge
    (i.e., row and column sums of each ij submatrix are 0)

    Parameters
    ----------
    J_ij : np.array
        Coupling matrix of size L x L x num_symbols x num_symbols
        that should be transformed into zero-sum gauge
    inplace : bool, optional (default: False)
        Modify original matrix (True), or return transformed
        matrix in a new matrix

    Returns
    -------
    J_ij_0 : np.array
        J_ij transformed into zero-sum gauge
    """
    L, L2, num_symbols, num_symbols2 = J_ij.shape
    assert L == L2 and num_symbols == num_symbols2

    if inplace:
        J_ij_0 = J_ij
    else:
        J_ij_0 = np.zeros((L, L, num_symbols, num_symbols))

    # go through all pairs of positions
    for i in range(L - 1):
        for j in range(i + 1, L):
            ij_mat = J_ij[i, j]

            # calculate matrix, row and column averages
            avg_ab = np.mean(ij_mat)

            # can't use axis argument of np.mean in numba,
            # so have to calculate rows/cols manually
            avg_a = np.zeros(num_symbols)
            avg_b = np.zeros(num_symbols)
            ij_mat_T = ij_mat.T

            for k in range(num_symbols):
                avg_a[k] = np.mean(ij_mat[k])
                avg_b[k] = np.mean(ij_mat_T[k])

            # subtract correction terms from each entry
            for a in range(num_symbols):
                for b in range(num_symbols):
                    J_ij_0[i, j, a, b] = (
                            ij_mat[a, b] - avg_a[a] - avg_b[b] + avg_ab
                    )
                    J_ij_0[j, i, b, a] = J_ij_0[i, j, a, b]

    return J_ij_0


def calculate_fn_score(j_matrix):
    # Function that calculates Frobenius Norm score of coupling matrices.
    q1 = j_matrix.shape[2]
    q2 = j_matrix.shape[3]
    L1 = j_matrix.shape[0]
    L2 = j_matrix.shape[1]
    scores = np.zeros((L1, L1))
    print("Size %d" % L1)
    assert L1 == L2 and q1 == q2
    for i in range(L1 - 1):
        for j in range(i + 1, L1):
            scores[i, j] = np.linalg.norm((j_matrix[i, j]), 'fro')
            scores[j, i] = scores[i, j]

    return scores


def apc(scores):
    L = scores.shape[0]
    # APC correction
    if L != scores.shape[1]:
        raise ValueError("Input matrix is not symmetric: {}".format(scores.shape))
    col_mean = np.mean(scores, axis=0) * L / (L - 1)
    matrix_mean = np.mean(scores) * L / (L - 1)

    apc_correction = np.dot(col_mean.reshape(L, 1), col_mean.reshape(1, L)) / matrix_mean
    corrnorms = scores - apc_correction
    corrnorms[np.diag_indices(L)] = 0
    return corrnorms


def load_matrix_matlab(filein, freq=False):
    # Function that reads in the couplings and fields Matlab matrix file.
    import h5py
    import scipy.io
    try:
        f = scipy.io.loadmat(filein)
        flag = 1
    except NotImplementedError:
        print("Using h5py...")
        f = h5py.File(filein, 'r')
        flag = 0

    if freq:
        return f['Pi']
    else:
        mat = {}
        for k, v in f.items():
            mat[k] = np.array(v)
        # h = mat['h']
        J = mat['J']
        if flag == 1:
            couplings = J  # use for scipy
        else:
            couplings = J.T  # use for h5py
        # fields = h.T
        # f.close()
        # return fields, couplings
        return couplings


def process_coupling_matrix_output_scores(matrixFile, freqFile):
    """

    :param freqFile:
    :param matrixFile: MATLAB mat file - coupling matrix
    :return: Dataframe with i, j, FN, and FN-apc columns

        1. Read MATLAB mat file then output coupling matrix
        2. Calculate FN score (save FN score matrix to file)
        3. Build contact map from matrix and save to file
        4. Merge FN score column with contact map with FN-apc column
    """

    if os.path.exists(matrixFile):
        # Read MATLAB mat file then output coupling matrix
        J_paired = load_matrix_matlab(matrixFile, freq=False)
        # Calculate FN score-paired and apply ap-correction
        fn_paired = calculate_fn_score(J_paired)
        fn_apc_paired = apc(fn_paired)
        # Build contact maps
        df_fn = cm_make(fn_paired, score='fn')
        df_fnapc = cm_make(fn_apc_paired, score='fn_apc')
        assert len(df_fnapc) == len(df_fn)
        _x = df_fnapc.merge(df_fn, on=['i', 'j'], how='inner')
        if freqFile:
            fi = load_matrix_matlab(freqFile, freq=True)
            DI = direct_information(J_paired, fi)
            df_di = cm_make(DI, score='di')
            _x = _x.merge(df_di, on=['i', 'j'], how='inner')

    else:
        print("Error: {} does not exist.".format(matrixFile))
        _x = 0
    return _x


def rank_hamming(df):
    """
    Rank DCA predictions by designated score and remove pairs with sequence distance < 5
    :param df: Dataframe with four columns i,j, fn_apc, and fn
    :return: Dataframe ranked by score and filtered
    """
    df_sorted = df.sort_values(ascending=False, by=['fn_apc'])
    df_hamming = df_sorted[abs(df_sorted['i'] - df_sorted['j']) > 5].reset_index(drop=True)
    return df_hamming


def test_functions():
    msa = '1EM8_D_1EM8_C'
    matrix_file = "coupling_matrices\\vanilla\\matrix_ising_{}.fas.mat".format(msa)
    freq_file = "frequency_files\\vanilla\\p_dist_{}.fas.mat".format(msa)
    mergedDF = process_coupling_matrix_output_scores(matrix_file, freq_file)
    rankedDF = rank_hamming(mergedDF)
    return len(mergedDF), mergedDF, len(rankedDF), rankedDF


def tilde_fields(J_ij, f_i, f_j):
    """Compute h_tilde fields of the two-site model.

    Parameters
    ----------
    J_ij : np.array
        Matrix of size num_symbols x num_symbols
        containing all coupling strengths of
        position pair (i, j).
    f_i : np.array
        Row i of single-site frequencies.
    f_j : np.array
        Row j of single-site frequencies.

    Returns
    -------
    np.array, np.array
        h_tilde fields of position i and j -
        both arrays of size 1 x num_symbols
    """
    _EPSILON = 1e-4
    diff = 1.0

    num_symbols = f_i.shape[0]

    h_tilde_i = np.full((1, num_symbols), 1 / float(num_symbols))
    h_tilde_j = np.full((1, num_symbols), 1 / float(num_symbols))

    while diff > _EPSILON:
        tmp_1 = np.dot(h_tilde_j, J_ij.T)
        tmp_2 = np.dot(h_tilde_i, J_ij)

        h_tilde_i_updated = f_i / tmp_1
        h_tilde_i_updated /= h_tilde_i_updated.sum()

        h_tilde_j_updated = f_j / tmp_2
        h_tilde_j_updated /= h_tilde_j_updated.sum()

        diff = max(
            np.absolute(h_tilde_i_updated - h_tilde_i).max(),
            np.absolute(h_tilde_j_updated - h_tilde_j).max()
        )

        h_tilde_i = h_tilde_i_updated
        h_tilde_j = h_tilde_j_updated

    return h_tilde_i, h_tilde_j


def direct_information(J_ij, f_i):
    """
    J_ij : np.array
        Matrix of size num_symbols x num_symbols
        containing all coupling strengths of
        position pair (i, j).
    f_i : np.array
        Matrix of size L x num_symbols
        containing column frequencies.
    np.array
        Matrix of size L x L
    """
    L, num_symbols = f_i.shape

    di = np.zeros((L, L))
    for i in range(L):
        for j in range(i + 1, L):
            # extract couplings relevant to
            # position pair (i, j)
            J = np.exp(J_ij[i, j])

            # compute two-site model
            h_tilde_i, h_tilde_j = tilde_fields(J, f_i[i], f_i[j])
            p_di_ij = J * np.dot(h_tilde_i.T, h_tilde_j)
            z = p_di_ij.sum()
            p_di_ij = p_di_ij / z

            # dot product of single-site frequencies
            # of columns i and j
            f_ij = np.dot(
                f_i[i].reshape((1, num_symbols)).T,
                f_i[j].reshape((1, num_symbols))
            )

            # finally, compute direct information as
            # mutual information associated to p_di_ij
            _TINY = 1.0e-100
            di[i, j] = di[j, i] = np.trace(
                np.dot(
                    p_di_ij.T,
                    np.log((p_di_ij + _TINY) / (f_ij + _TINY))
                )
            )

    return di


# l, m, l2, r = test_functions()

# Calculate FN-average of every replicate in scrambled set and then apc
# fn_scrambled_avg = average_fn(msaName, rep)
# fn_apc_scrambled = apc(fn_scrambled_avg)
# df_scrambled = cm_make(msaName, fn_scrambled_avg)
# df_fn_apc_scrambled = cm_make(msaName, fn_apc_scrambled)

# New score based on difference between paired and scrambled scores - with apc applied beforehand
# fni_apc_before = fn_apc_paired - fn_apc_scrambled
# df_fd = cm_make(msaName, fni_apc_before)

# New score based on difference between paired and scrambled scores - with apc applied after the difference
# fni = fn_paired - fn_scrambled_avg
# fni_apc_after = apc(fni)
# df_difference = cm_make(msaName, fn_apc_difference)

# np.save('{}matrix_FNi_{}.npy'.format(outDir, msaName), fni)
# np.save('{}matrix_FNi_apc_before_{}.npy'.format(outDir, msaName), fni_apc_before)
# np.save('{}matrix_FNi_apc_after_{}.npy'.format(outDir, msaName), fni_apc_after)
