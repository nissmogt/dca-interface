import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
# 5 plots


def plot_npairs_vs_zscore(df, labels, msa, i, score):
    plt.figure(100000 * (i+1), figsize=(10, 5), dpi=100)
    plt.scatter(range(len(df[0])), y='zscore', data=df[0], c='orange', label=labels[0])
    plt.scatter(range(len(df[1])), y='zscore', data=df[1], marker='.', label=labels[1], alpha=0.7)
    plt.xlabel("ranked DCA contacts")
    plt.ylabel("z-score")
    plt.title("{}".format(msa))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    # plt.show()
    imgname = "{}{}_{}_zscore_{}_vs_npairs_vanilla_vs_nonbonded_restraint_{}A.png".format(imdir, msa,
                                                                                          pairtype, score, rcutoff)
    plt.savefig(imgname, dpi=900, bbox_inches='tight')
    plt.close(100000*(i+1))


def plot_histogram_scores(dataDf, labelList, msaName, idx, scoreType):
    plt.figure(20000 * (idx + 1), figsize=(10, 5), dpi=100)
    plt.hist(dataDf[0][scoreType], bins=30, edgecolor='black', label=labelList[0], color='orange')
    plt.hist(dataDf[1][scoreType], bins=30, edgecolor='black', label=labelList[1], alpha=0.5)
    plt.yscale('log')
    plt.legend(loc='best')
    plt.xlabel("FN score")
    plt.ylabel("counts")
    plt.title("{}".format(msaName))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    imgname = "{}{}_{}_hist_{}_nonbonded_{}A_vs_vanilla.png".format(imdir, msaName, pairtype, scoreType, rcutoff)
    plt.savefig(imgname, dpi=900, bbox_inches='tight')
    plt.close(20000 * (idx + 1))


def plot_distance_vs_zscore(df, labels, msa, i, score):
    plt.figure(10000 * (i+1), figsize=(10, 5), dpi=100)
    plt.scatter(x='d', y='zscore', data=df[0], c='orange', label=labels[0])
    plt.scatter(x='d', y='zscore', data=df[1], marker='.', label=labels[1], alpha=0.7)
    plt.xlabel("distance")
    plt.ylabel("z-score")
    plt.title("{}".format(msa))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    # plt.show()
    imgname = "{}{}_{}_distance_vs_zscore_{}_vanilla_vs_nonbonded_restraint_{}A.png".format(imdir, msa, pairtype,
                                                                                            score, rcutoff)
    # imgname = "{}{}_noAPC_monomer_distance_vs_zscore_vanilla_vs_nonbonded_restraint12A.png".format(imdir, msa)
    plt.savefig(imgname, dpi=900, bbox_inches='tight')
    plt.close(10000*(i+1))


def plot_std(std_list, labelList, msaName, idx, scoreType):
    plt.figure(100 * (idx + 1), figsize=(10, 5), dpi=100)
    plt.plot(range(len(std_list[0])), std_list[0], c='orange', label=labelList[0])
    plt.plot(range(len(std_list[1])), std_list[1], label=labelList[1])
    plt.xlabel("ranked DCA contacts")
    plt.ylabel("standard deviation (zscore)")
    plt.title("{}".format(msaName))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    imgname = "{}{}_removed_dca_vs_std_{}_{}_nonbonded_{}A_vs_vanilla.png".format(imdir, msaName, pairtype,
                                                                                  scoreType, rcutoff)
    plt.savefig(imgname, dpi=900, bbox_inches='tight')
    plt.close(100 * (idx + 1))


def plot_mean(mean_list, labelList, msaName, idx):
    plt.figure(2000 * (idx + 1), figsize=(10, 5), dpi=100)
    plt.plot(range(len(mean_list[0])), mean_list[0], c='orange', label=labelList[0])
    plt.plot(range(len(mean_list[1])), mean_list[1], label=labelList[1])
    plt.xlabel("ranked DCA contacts")
    plt.ylabel("mean (zscore)")
    plt.title("{}".format(msaName))
    plt.legend(loc='best')
    plt.grid(axis='both', alpha=0.5)
    imgname = "{}{}_removed_dca_vs_mean_{}_nonbonded_{}A_vs_vanilla.png".format(imdir, msaName, pairtype, rcutoff)
    plt.savefig(imgname, dpi=900, bbox_inches='tight')
    plt.close(2000 * (idx + 1))


def stats_removed_pairs(df):
    _std = np.zeros(len(df))
    _mean = np.zeros(len(df))
    for p in range(len(_std)):
        _std[p] = np.std(df['zscore'][p:])
        _mean[p] = np.mean(df['zscore'][p:])
    return _std, _mean


def zscore_calc(data, n_pairs, score='fn_apc', which='monomer', stats=False):
    if which == 'monomer':
        _df = data[data['chain_1'] == data['chain_2']].reset_index()
    if which == 'interface':
        _df = data[data['chain_1'] != data['chain_2']].reset_index()
    else:
        _df = data
    mean = np.mean(data[score])
    std = np.std(data[score])
    _df['zscore'] = (_df[score] - mean) / std
    if stats:
        _s, _m = stats_removed_pairs(_df)
    else:
        _s = 0
        _m = 0
    return _df, _s, _m


results_directory = ["vanilla_results\\", "nonbonded_restraints_results\\12A\\"]
# results_directory = ["vanilla_results\\", "sasa_restraints_results\\sasa_5\\"]
imdir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2021\\JAN_2021\\"

dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
          "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
          '5MU7_B_5MU7_A', '5M72_A_5M72_B', 'HKRR_C_HKRR_A']
label1 = results_directory[0].split('\\')[0]
label2 = results_directory[1].split('\\')[0]
labels = [label1, label2]

# pairtype = 'interface'
rcutoff = 12
score = 'fn'

# i = 1
for pairtype in ['interface', 'monomer']:
    for i in range(11):
        nPairs = 50
        msa = dimers[i]
        data1 = "{}FN_all_{}_mapped_aa_dist.txt".format(results_directory[0], msa)
        data2 = "{}FN_all_{}_mapped_aa_dist.txt".format(results_directory[1], msa)
        df1_all = pd.read_csv(data1, delimiter='\t')
        df2_all = pd.read_csv(data2, delimiter='\t')
        df1, s1, m1 = zscore_calc(df1_all, nPairs, which=pairtype, stats=True, score=score)
        df2, s2, m2 = zscore_calc(df2_all, nPairs, which=pairtype, stats=True, score=score)
        df = [df1, df2]
        s = [s1, s2]
        m = [m1, m2]
        plot_std(s, labels, msa, i, score)
        plot_mean(m, labels, msa, i)
        plot_distance_vs_zscore(df, labels, msa, i, score)
        plot_npairs_vs_zscore(df, labels, msa, i, score)
        plot_histogram_scores(df, labels, msa, i, score)


# header = "i\tj\tscore\tdist_aa\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id\tratio_i\tratio_j\tzscore"
    # out1 = "{}FN_inter_{}_mapped_aa_dist_sasa_zscore.txt".format(results_directory[0], msa)
    # out2 = "{}FN_inter_{}_mapped_aa_dist_sasa_zscore.txt".format(results_directory[1], msa)
    # df1.to_csv(out1, sep='\t', index=False, header=header, float_format='%.5f')
    # df2.to_csv(out2, sep='\t', index=False, header=header, float_format='%.5f')
