import numpy as np
import matplotlib.pyplot as plt

from plot_cm import plot_cm
from sasa_calc import tpr_dca_sasa, sasa_filter, draw_sasa_res, draw_dca_sasa, sasa_affect_plot
from make_distance_matrix_dca_map import vectorize_pdb_contacts, vectorize_dca_contacts, tpr_top_pairs

pdb_file = 'files\\XDHB_XDHC\\pdb\\3hrd_CD_pdbfixed.pdb'
dca_file = 'files\\XDHB_XDHC\\map_cm_XDHB_XDHC.txt'
monomer_file = 'files\\XDHB_XDHC\\sasa\\sasa_ps14_mon_3hrd_CD_pdbfixed.txt'

# Filter DCA predictions using sasa dictionary
n_pairs = 15
threshold_sasa = [0, 1, 10, 50]
cutoff = 10
# XDHB-XDHC Parameters
length_a = 291
dimer_length = 452
chains = ['A', 'C']




# tpr = []
# pdb_flat_matrix = vectorize_pdb_contacts(pdb_file, dimer_length,
#                                         chain, calpha_cutoff)
# Plot the effect of SASA filter on number of DCA pairs
# sasa_affect_plot(dca_file, n_contacts, monomer_file, threshold_sasa)

# Plot True positive rate as function of number of SASA-filtered DCA pairs
# tpr = tpr_dca_sasa(pdb_flat_matrix, dca_file, n_contacts, monomer_file,
#                   threshold_sasa, dimer_length, chain, calpha_cutoff)

# Make draw tcl script for sasa threshold
# draw_sasa_res(monomer_file, 1)
# draw_dca_sasa(dca_file, n_contacts, length_a, chain, monomer_file, 0)
# draw_dca_sasa(dca_file, n_contacts, length_a, chain, monomer_file, 1)
# draw_dca_sasa(dca_file, n_contacts, length_a, chain, monomer_file, 10)
# df_dca = pd.DataFrame(dca_sasa_ratio, columns=index)
# plot_dca_sasa_filter(pdb_file, df_dca, chain, 10, cutoff)
