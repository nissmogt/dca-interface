import numpy as np
import matplotlib.pyplot as plt

from sasa_calc import tpr_dca_sasa, sasa_filter, draw_sasa_res, draw_dca_sasa, sasa_affect_plot
from dca_performance import vectorize_pdb_contacts, vectorize_dca_contacts, loop_tpr

pdb_file = 'files\\YBGJ_YBGK\\pdb\\pdbfixed_corrected.pdb'
dca_file = 'files\\YBGJ_YBGK\\map_cm_YBGJ_YBGK.txt'
monomer_file = 'files\\YBGJ_YBGK\\pdb\\sasa_pdbfixed_corrected_BA.txt'

# Filter DCA predictions using sasa dictionary
n_contacts = 15
threshold_sasa = [0, 0.5, 1, 2]
calpha_cutoff = 15
# XDHB-XDHC Parameters
# length_a = 291
# dimer_length = 452
# chain = ['C', 'D']
# FDNI-YSAA Parameters
# length_a = 217
# dimer_length = 511
# chain = ['C', 'B']
# YBGJ-YBGK Parameters
length_a = 211
dimer_length = 501
chain = ['B', 'A']
tpr = []
pdb_flat_matrix = vectorize_pdb_contacts(pdb_file, dimer_length,
                                         chain, calpha_cutoff)
# Plot the effect of SASA filter on number of DCA pairs
sasa_affect_plot(dca_file, n_contacts, monomer_file, threshold_sasa)

# Plot True positive rate as function of number of SASA-filtered DCA pairs
tpr = tpr_dca_sasa(pdb_flat_matrix, dca_file, n_contacts, monomer_file,
                   threshold_sasa, dimer_length, chain, calpha_cutoff)

# Make draw tcl script for sasa threshold
draw_sasa_res(monomer_file, 1)
draw_dca_sasa(dca_file, n_contacts, length_a, chain, monomer_file, 0)
draw_dca_sasa(dca_file, n_contacts, length_a, chain, monomer_file, 1)
draw_dca_sasa(dca_file, n_contacts, length_a, chain, monomer_file, 10)
# df_dca = pd.DataFrame(dca_sasa_ratio, columns=index)
# plot_dca_sasa_filter(pdb_file, df_dca, chain, 10, cutoff)
