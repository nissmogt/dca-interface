import glob
import os
import numpy as np
import matplotlib.pyplot as plt

img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\AUG_2020\\"

results_dir = 'results\\'
high_tp = ["1EFP_A_1EFP_B", "1QOP_A_1QOP_B", "1I1Q_A_1I1Q_B", "2ONK_A_2ONK_C", "2NU9_A_2NU9_B", "3G5O_A_3G5O_B"]
low_tp = ["3A0R_A_3A0R_B", "3RPF_A_3RPF_C"]
l_length = len(high_tp) + len(low_tp)

results_list = []
tp_list = []
for name in high_tp+low_tp:
    results_list.append('{}mapped_cm_FNi_apc_{}.txt'.format(results_dir, name))
    tp_list.append('{}TP_{}.txt'.format(results_dir, name))

n_pairs = 50
scores = []
tp = []
for i in range(len(results_list)):
    dca_pairs = np.loadtxt(results_list[i])
    tp_i = np.loadtxt(tp_list[i])
    scores.append(dca_pairs.transpose()[2][:n_pairs])
    tp.append(100 * tp_i[:n_pairs] / n_pairs)
avg_corrected = []
fig, ax = plt.subplots(len(scores), sharex='col', sharey='row')
for j in range(len(scores)):
    name = os.path.basename(results_list[j]).strip("mapped_cm_FNi_apc_")
    avg_corrected.append(scores[j] - np.average(scores[j]))
    # std_dev = np.std(scores[j])
    ax[j].hist(avg_corrected[j], bins=10, label="{}, AVG_FN= {:.3f}, TP= {}%".format(name, np.average(scores[j]), tp[j][-1]),
               edgecolor='black')
    ax[j].legend(loc='best')
    # np.savetxt("avg_scores_TP_{}".format(name), [avg_corrected, tp[j][-1]], fmt="%.3f")
plt.xlabel("Average-corrected score")
plt.show()
