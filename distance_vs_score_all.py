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
                                plot_tpr_as_function_of_scores)

img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"

n_pairs = 50
cutoff = 15
systemFiles = glob.glob("results\\*_inter_mapped_aa_dist_top{}_tp.txt".format(n_pairs))
nSystems = len(systemFiles)
systemsList = []
systemNames = []

# add all systems to a list
for idx, system in enumerate(systemFiles):
    df_system = pd.read_csv(system, delimiter='\t')
    systemName = os.path.basename(system).split("_")
    systemName = '_'.join(systemName[1:5])
    systemNames.append("{}.fas\n".format(systemName))
    systemsList.append(df_system)
    print(systemName)

f = open("benchmark_systems_count{}.txt".format(idx+1), 'w')
f.writelines(systemNames)
f.close()
# define max and min of bin from all system score extrema
maxBinvalue = np.zeros(nSystems)
minBinvalue = np.zeros(nSystems)
for i in range(nSystems):
    maxBinvalue[i] = systemsList[i]["score"].max()
    minBinvalue[i] = systemsList[i]["score"].min()

minBinvalue = min(minBinvalue)
maxBinvalue = max(maxBinvalue)
dx = 0.025  # width of bin
binned_scores = np.arange(minBinvalue, maxBinvalue + dx, dx)
nBins = len(binned_scores)
TPR = np.zeros(nBins)
norm = np.zeros(nBins)
tpr_unbinned = np.zeros(n_pairs)
for i in range(nSystems):
    for pairs in range(len(systemsList[i])):
        FN = systemsList[i]["score"][pairs]
        distance = systemsList[i]["dist"][pairs]
        if distance < cutoff:
            tpr_unbinned[pairs] += 1.0 / nSystems
        for binIndex in range(nBins):
            if binned_scores[binIndex] <= FN < binned_scores[binIndex + 1]:
                norm[binIndex] += 1  # count all pairs that lie in this bin
                if distance < cutoff:
                    TPR[binIndex] += 1  # add count to tpr bin
# divide TPR by total number of pairs for each system
TPR_norm = np.divide(TPR, norm, where=norm != 0)

# plot_tpr_as_function_of_pairs(n_pairs, tpr_unbinned, cutoff, nSystems)
# plot_tpr_as_function_of_scores(TPR_norm, binned_scores, dx, n_pairs, cutoff, nSystems)
# plot_npairs_as_function_of_scores(norm, binned_scores, dx, n_pairs, cutoff, nSystems)
