def filter_gaps(msa_obj, gap_cutoff=0.5):
    """filters alignment to remove gappy positions"""
    # msa must be a list
    import numpy as np
    alphabet = "ARNDCQEGHILKMFPSTWYV-"
    states = len(alphabet)
    tmp = (msa_obj == states - 1).astype(np.float)
    non_gaps = np.where(np.sum(tmp.T, -1).T / msa_obj.shape[0] < gap_cutoff)[0]
    return msa_obj[:, non_gaps], non_gaps
