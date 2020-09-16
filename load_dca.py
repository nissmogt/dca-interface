def load_dca(msa_name, results_dir):
    import os
    import numpy as np
    import pandas as pd
    print("\n\t-- LOADING: {} --".format(msa_name))

    # dca_array = '{}fn_{}_plmdca_rt2.txt'.format(results_dir, msa_name)
    dca_array = "{}FN_{}.txt".format(results_dir, msa_name)
    print("(load_dca)\t{}".format(dca_array))
    df_dca = pd.read_csv(dca_array, delimiter=',', names=['i', 'j', 'score'],
                         usecols=(0, 1, 2))
    df_dca = df_dca.sort_values(ascending=False, by=['score'])
    df_dca = df_dca[abs(df_dca["i"] - df_dca["j"]) > 5]    # ensure sequence separation of 5
    outfile = "{}FN_{}_ranked.txt".format(results_dir, msa_name)
    np.savetxt(outfile, df_dca, fmt='%d\t%d\t%f')
    return df_dca
