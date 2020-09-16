def mk_msa(seqs):
    """converts list of sequences to msa"""
    import numpy as np
    from aa2num import aa2num
    from effective_seq_calc import effective_seq_calc
    msa_ori = []
    for seq in seqs:
        msa_ori.append([aa2num(aa) for aa in seq])
    msa_ori = np.array(msa_ori)

    # remove positions with more than > 50% gaps
    # msa, v_idx = self.filt_gaps(msa_ori, 0.5)
    msa = msa_ori

    # compute effective weight for each sequence
    msa_weights = effective_seq_calc(msa, 0.8)

    # compute effective number of sequences
    ncol = msa.shape[1]  # length of sequence
    # w_idx = v_idx[np.stack(np.triu_indices(ncol, 1), -1)]

    return {"msa_ori": msa_ori,
            "msa": msa,
            "weights": msa_weights,
            "neff": np.sum(msa_weights),
            "nrow": msa.shape[0],
            "ncol": ncol,
            "ncol_ori": msa_ori.shape[1]}
    # "v_idx": v_idx,
    # "w_idx": w_idx
    # }
