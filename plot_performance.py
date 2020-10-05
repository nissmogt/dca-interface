import numpy as np
from make_distance_matrix_dca_map import plot_performance

# Static parameters
n_contacts = 100
L = 452
chain = ['C', 'D']
cutoff = 15

# Load files into array
dcafile = "files\\map_scan_XDHB_XDHC_CD.txt"
# gremlinfile = "files\\gremlin_map_interface_xdhb.txt"
pdbfile = "files\\XDHB_XDHC\\3hrd_CD_pdbfixed.pdb"
dca_array = np.loadtxt(dcafile, usecols=(0, 1), dtype=int)

plot_performance(pdbfile, dca_array, n_contacts, L, chain)
