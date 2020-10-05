def load_dca_to_df(msa_name, results_dir):
    import numpy as np
    import pandas as pd
    print("\n\t-- LOADING: {} --".format(msa_name))

    dca_raw_output = "{}FN_{}.fas.txt".format(results_dir, msa_name)
    print("(load_dca)\t{}".format(dca_raw_output))

    df_dca = pd.read_csv(dca_raw_output, delimiter=',', names=['i', 'j', 'score'],
                         usecols=(0, 1, 2))
    df_dca = df_dca.sort_values(ascending=False, by=['score'])
    df_dca = df_dca[abs(df_dca["i"] - df_dca["j"]) > 5]    # ensure sequence separation of 5
    outfile = "{}FN_{}_ranked.txt".format(results_dir, msa_name)
    np.savetxt(outfile, df_dca, fmt='%d\t%d\t%f', comments='')
    return df_dca
