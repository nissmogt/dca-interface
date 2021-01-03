import matplotlib.pyplot as plt
"""
Calculates true positives based on pair distance cutoff.
Plots true positive rate as a function of ranked DCA score.
Plots true positive rate (true positives normalized by total predicted pairs) as a function of FN-apc score.
"""
import glob
import os
import pandas as pd
import numpy as np
from analysis_functions import (plot_npairs_as_function_of_scores, plot_tpr_as_function_of_pairs,
                                compare_two_tpr, plot_tp_as_function_of_scores, plot_fp_as_function_of_scores)


def compare_fni_fn(n_pairs, threshold, width):
    tpr_n, tpr, normalize, binned_scores, fn_system_list = calculate_true_positives(n_pairs, threshold, width, False)
    fni_tpr_n, fni_tpr, fni_norm, fni_bin_scores, fni_system_list = calculate_true_positives(n_pairs, threshold, width,
                                                                                             True)
    n_sys = len(fni_system_list)
    compare_two_tpr(n_pairs, fni_tpr, tpr, threshold, n_sys)


def normalize_tpr(true_positive_array, norm, nBins, binwidth):
    TPR_norm = np.zeros(nBins)
    sum_norm = sum(norm)

    return TPR_norm / sum_norm*binwidth


def make_bins(nSystems, systemsList, dx):
    """
    Define max and min of bin from all system score extrema

    :param nSystems:
    :param systemsList:
    :param dx: string; Bin width
    :return: np.array binned scores from min to max in steps of dx
    """
    maxBinvalue = np.zeros(nSystems)
    minBinvalue = np.zeros(nSystems)
    for i in range(nSystems):
        maxBinvalue[i] = systemsList[i]["score"].max()
        minBinvalue[i] = systemsList[i]["score"].min()

    minBinvalue = min(minBinvalue)
    maxBinvalue = max(maxBinvalue)
    binned_scores = np.arange(0, maxBinvalue + dx, dx)
    print("min bin: {}, max bin: {}".format(minBinvalue, maxBinvalue))
    return binned_scores


def make_system_list(n_pairs, fni=False, from_list=None):
    scrambledDir = "scrambled_results\\fni_matrices\\"
    # scrambledDir = "scrambled_results\\from_averaging_jmatrices\\apc\\"
    if from_list:
        systemFiles = []
        for idx, msa in enumerate(from_list):
            if fni:
                systemFiles.append("{}{}_FNi_mapped_top50.txt".format(scrambledDir, msa))
            else:
                systemFiles.append("vanilla_results\\FN_all_{}_mapped_aa_dis_sasa.txt".format(msa))
                # systemFiles.append("restraint_results\\FN_{}_inter_mapped_aa_dist_top50.txt".format(msa))
                # systemFiles.append("results\\FN_{}_inter_mapped_aa_dist_top50.txt".format(msa))
    else:
        if fni:
            systemFiles = glob.glob("{}*_FNi_mapped_top50.txt".format(scrambledDir))
            # systemFiles = glob.glob("{}FNi_apc_*_inter_mapped_aa_dist_top{}.txt".format(scrambledDir, n_pairs))
        else:
            # systemFiles = glob.glob("restraint_results\\FN_*_inter_mapped_aa_dist_top50.txt")
            systemFiles = glob.glob("vanilla_results\\FN_all_*_mapped_aa_dis_sasa.txt")
    systemsList = []
    systemNames = []
    count = 0
    for idx, system in enumerate(systemFiles):
        systemName = os.path.basename(system).split("_")
        if fni:
            # systemName = '_'.join(systemName[2:6])
            systemName = '_'.join(systemName[:4])
            dist_file = "{}{}_FNi_mapped_top50.txt".format(scrambledDir, systemName)
            # dist_file = "{}FNi_apc_{}_inter_mapped_aa_dist_top{}.txt".format(scrambledDir, systemName, n_pairs)
        else:
            systemName = '_'.join(systemName[2:6])
            dist_file = "vanilla_results\\FN_all_{}_mapped_aa_dis_sasa.txt".format(systemName)
        if os.path.exists(dist_file):
            count += 1
            df_system = pd.read_csv(system, delimiter='\t')
            # systemName = '_'.join(systemName[2:6])
            systemNames.append("{}.fas\n".format(systemName))
            systemsList.append(df_system[:n_pairs])
            # print(systemName)

    # write list of systems analyzed into a file
    if fni:
        f = open("benchmark_systems_FNi_count{}.txt".format(count), 'w')
    else:
        f = open("benchmark_systems_FN_count{}.txt".format(count), 'w')
    f.writelines(systemNames)
    f.close()
    return systemsList, systemNames


def batch_calculate_true_positives(n_pairs, cutoff, dx, fni=False, msalist=False):
    if msalist:
        systemsList, systemNames = make_system_list(n_pairs, fni=fni, from_list=msalist)
    else:
        systemsList, systemNames = make_system_list(n_pairs, fni=fni)
    nSystems = len(systemsList)
    # add all systems to a list
    binned_scores = make_bins(nSystems, systemsList, dx)
    nBins = len(binned_scores)
    tp = np.zeros((nSystems, nBins))
    fp = np.zeros((nSystems, nBins))
    systems = []
    for i in range(nSystems):
        for pairs in range(len(systemsList[i])):
            FN = systemsList[i]["score"][pairs]
            distance = systemsList[i]["dist_aa"][pairs]
            for binIndex in range(nBins - 1):
                if binned_scores[binIndex] <= FN < binned_scores[binIndex + 1]:
                    if distance <= cutoff:
                        tp[i, binIndex] += 1  # true positive
                    else:
                        fp[i, binIndex] += 1

        systems.append([systemNames[i].strip(".fas\n"), tp + fp])

    # if fni:
    #     np.savetxt('table_{}systems_fni_scores_top{}.txt'.format(nSystems, n_pairs), systems, fmt='%s',
    #                delimiter='\t', comments='')
    # else:
    #     np.savetxt('table_{}systems_fn_scores_top{}.txt'.format(nSystems, n_pairs), systems, fmt='%s',
    #                delimiter='\t', comments='')

    # normalized_tp = tp[0]
    # normalized_fp = fp[0]
    # normalized_tp = np.nan_to_num(tp[0] / (sum(tp[0])))
    # normalized_fp = np.nan_to_num(fp[0] / (sum(fp[0])))

    # normalized_counts = counts / (sum(counts)*dx)
    return tp, fp, binned_scores, systems


def calculate_true_positives(msa, n_pairs, cutoff, dx, fni=False):
    if fni:
        dist_file = "scrambled_results\\fni_matrices\\{}_FNi_mapped_top50.txt".format(msa)
    else:
        # dist_file = "restraint_results\\FN_{}_inter_mapped_aa_dist_top50.txt".format(msa)
        dist_file = "results\\FN_{}_inter_mapped_aa_dist_top50.txt".format(msa)

    df_system = pd.read_csv(dist_file, delimiter='\t', header=0)[:n_pairs]
    # add all systems to a list
    maxBinvalue = df_system["score"].max()
    minBinvalue = df_system["score"].min()

    binned_scores = np.arange(0, 1.4, dx)
    print("min bin: {}, max bin: {}".format(minBinvalue, maxBinvalue))

    nBins = len(binned_scores)
    TPR = np.zeros(nBins)
    counts = np.zeros(nBins)
    tpr_unbinned = np.zeros(n_pairs)
    high_scores = np.zeros(nBins)
    systems = []
    # for i in range(nSystems):
    i = 0
    i_counts = np.zeros(nBins)
    tp_count = 0
    for pairs in range(len(df_system)):
        FN = df_system["score"][pairs]
        distance = df_system["dist_aa"][pairs]
        # if distance <= cutoff:
        #     tp_count += 1
        #     tpr_unbinned[pairs] += 1
        for binIndex in range(nBins - 1):
            if binned_scores[binIndex] <= FN < binned_scores[binIndex + 1]:
                counts[binIndex] += 1  # count all pairs that lie in this bin
                if distance <= cutoff:
                    TPR[binIndex] += 1  # add count to tpr bin

    normalized_tpr = TPR / sum(counts)*dx
    return normalized_tpr, tpr_unbinned, counts, binned_scores, systems


def ppv(df_system, n_pairs, cutoff, fni=False):

    tp = 0
    fp = 0
    ppv_array = np.zeros(n_pairs)
    for pairs in range(len(df_system)):
        FN = df_system["score"][pairs]
        distance = df_system["dist_aa"][pairs]
        if distance <= cutoff:
            tp += 1
        else:
            fp += 1
        ppv_array[pairs] = tp / (tp + fp)

    return ppv_array


nPairs = 20000
contactCutoff = 12
binWidth = 0.05

dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
          "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
          '5MU7_B_5MU7_A']


# for i in range(len(dimers)):
# Plot true positive histograms
m = dimers[-1]
dist_fni = "scrambled_results\\fni_matrices\\{}_FNi_mapped_top10000.txt".format(m)
df_dist_fni = pd.read_csv(dist_fni, delimiter='\t', header=0)[:nPairs]

dist_no_rest = "vanilla_results\\FN_all_{}_mapped_aa_dist.txt".format(m)
# df_dist_no_rest = pd.read_csv(dist_no_rest, delimiter='\t', header=0)[:nPairs]
df_dist_no_rest = pd.read_csv(dist_no_rest, delimiter='\t', header=0)
df_dist_no_rest = df_dist_no_rest[df_dist_no_rest["chain_1"] != df_dist_no_rest["chain_2"]][:nPairs].reset_index(drop=True)
dist_nonbonded_APC = "nonbonded_restraints_results\\APC\\8A\\FN_all_{}_mapped_aa_dist.txt".format(m)
# df_dist_nonbonded_APC = pd.read_csv(dist_nonbonded_APC, delimiter='\t', header=0)[:nPairs]
df_dist_nonbonded_APC = pd.read_csv(dist_nonbonded_APC, delimiter='\t', header=0)
df_dist_nonbonded_APC = df_dist_nonbonded_APC[df_dist_nonbonded_APC["chain_1"] != df_dist_nonbonded_APC["chain_2"]][:nPairs].reset_index(drop=True)
dist_nonbonded_12APC = "nonbonded_restraints_results\\APC\\12A\\FN_all_{}_mapped_aa_dist.txt".format(m)
# df_dist_nonbonded_12APC = pd.read_csv(dist_nonbonded_12APC, delimiter='\t', header=0)[:nPairs]
df_dist_nonbonded_12APC = pd.read_csv(dist_nonbonded_12APC, delimiter='\t', header=0)
df_dist_nonbonded_12APC = df_dist_nonbonded_12APC[df_dist_nonbonded_12APC["chain_1"] != df_dist_nonbonded_12APC["chain_2"]][:nPairs].reset_index(drop=True)
dist_nonbonded_15APC = "nonbonded_restraints_results\\APC\\15A\\FN_all_{}_mapped_aa_dist.txt".format(m)
df_dist_nonbonded_15APC = pd.read_csv(dist_nonbonded_15APC, delimiter='\t', header=0)[:nPairs]
# df_dist_nonbonded_15APC = pd.read_csv(dist_nonbonded_15APC, delimiter='\t', header=0)
df_dist_nonbonded_15APC = df_dist_nonbonded_15APC[df_dist_nonbonded_15APC["chain_1"] != df_dist_nonbonded_15APC["chain_2"]][:nPairs].reset_index(drop=True)
dist_nonbonded_20APC = "nonbonded_restraints_results\\APC\\20A\\FN_all_{}_mapped_aa_dist.txt".format(m)
# df_dist_nonbonded_20APC = pd.read_csv(dist_nonbonded_20APC, delimiter='\t', header=0)[:nPairs]
df_dist_nonbonded_20APC = pd.read_csv(dist_nonbonded_20APC, delimiter='\t', header=0)
df_dist_nonbonded_20APC = df_dist_nonbonded_20APC[df_dist_nonbonded_20APC["chain_1"] != df_dist_nonbonded_20APC["chain_2"]][:nPairs].reset_index(drop=True)
# ppv_no_restraint = ppv(df_dist_no_rest, nPairs, contactCutoff)
# ppv_fni = ppv(df_dist_fni, nPairs, contactCutoff)
# ppv_nonbonded_APC = ppv(df_dist_nonbonded_APC, nPairs, contactCutoff)
# ppv_nonbonded_12APC = ppv(df_dist_nonbonded_12APC, nPairs, contactCutoff)
# ppv_nonbonded_15APC = ppv(df_dist_nonbonded_15APC, nPairs, contactCutoff)
# ppv_nonbonded_20APC = ppv(df_dist_nonbonded_20APC, nPairs, contactCutoff)
ss = 5
# plt.figure(figsize=(10, 5), dpi=90)
# x = range(1, nPairs+1)
# plt.scatter(x='dist_aa', y='score', data=df_dist_no_rest, label='no restraints')
# plt.scatter(x='dist_aa', y='score', data=df_dist_fni, label='fn-fn_scrambled')
# plt.scatter(x='dist_aa', y='score', data=df_dist_nonbonded_APC, label='nonbonded 8APC', marker='s')
# plt.scatter(x='dist_aa', y='score', data=df_dist_nonbonded_12APC, label='nonbonded 12APC', marker="*")
# plt.scatter(x='dist_aa', y='score', data=df_dist_nonbonded_15APC, label='nonbonded 15APC', marker="^")
# plt.scatter(x='dist_aa', y='score', data=df_dist_nonbonded_20APC, label='nonbonded 20APC', marker="x")
# plt.plot(x, ppv_no_restraint, label="no restraints", marker='s')
# plt.plot(x, ppv_fni, label="FN - FN_scrambled", marker='d', markersize=ss)
# plt.plot(x, ppv_nonbonded_APC, label="nonbonded restraints APC 8A", marker='^', markersize=ss)
# plt.plot(x, ppv_nonbonded_12APC, label="nonbonded restraints APC 12A", marker='*', markersize=ss)
# plt.plot(x, ppv_nonbonded_15APC, label="nonbonded restraints APC 15A", marker='x', markersize=ss)
# plt.plot(x, ppv_nonbonded_20APC, label="nonbonded restraints APC 20A", markersize=ss)
# plt.legend(loc='best')
# plt.grid(axis='both', alpha=0.3)
# plt.yticks(np.arange(0, 3.0, 0.1))
# plt.title(m)
# plt.yscale('log')
# plt.xlabel('rank ordered DCA contacts')
# plt.ylabel('PPV')
# plt.xlabel('distance')
# plt.ylabel('score')
# plt.xlim(0, 60)
# plt.ylim(0, 3.0)
# plt.show()


msaFile = dimers
fni_flag = False
tp_norm, fp_norm, scores_binned, s = batch_calculate_true_positives(nPairs, contactCutoff, binWidth, fni_flag,
                                                                    msalist=msaFile)
# nSys = len(s)
# if nSys != 0:
#     tp_norm = sum(tp_norm) / nSys
#     fp_norm = sum(fp_norm) / nSys
# for i in range(nSys):
#     plot_tp_as_function_of_scores(tp_norm[i], scores_binned[i], binWidth, nPairs, contactCutoff, nSys, fni_flag)
# plot_fp_as_function_of_scores(fp_norm, scores_binned, binWidth, nPairs, contactCutoff, nSys, fni_flag)
# plot_npairs_as_function_of_scores(tp_norm + fp_norm, scores_binned, binWidth, nPairs, contactCutoff, nSys, fni_flag)
# plot_tpr_as_function_of_pairs(nPairs, unbinned_tpr, contactCutoff, nSys, fni_flag)

