def effective_seq_calc(msa_list, eff_cutoff=0.8):
    """compute effective weight for each sequence"""
    from scipy.spatial.distance import pdist, squareform
    import numpy as np
    ncol = msa_list.shape[1]

    # pairwise identity
    msa_sm = 1.0 - squareform(pdist(msa_list, "hamming"))

    # weight for each sequence
    msa_w = (msa_sm >= eff_cutoff).astype(np.float)
    msa_w = 1 / np.sum(msa_w, -1)

    return msa_w
