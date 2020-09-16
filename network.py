import networkx as nx
import numpy as np

results_dir = "results\\"
resi, resj, scores = np.loadtxt(results_dir+"mapped_cm_FNi_apc_1EFP_A_1EFP_B.txt", unpack=True)
G = nx.Graph()


