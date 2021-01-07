import matplotlib.pyplot as plt
import collections
import numpy as np
import networkx as nx
import pandas as pd

dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
          "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
          '5MU7_B_5MU7_A']
n=5
data = "vanilla_results//FN_all_{}_mapped_aa_dist.txt".format(dimers[n])
df_dca = pd.read_csv(data, header=0, delimiter='\t')
all_residues = df_dca.iloc[:, :2].to_numpy()
monomer_residues = df_dca[df_dca["chain_1"] == df_dca["chain_2"]].iloc[:, :2].to_numpy()
interface_residues = df_dca[df_dca["chain_1"] != df_dca["chain_2"]].iloc[:, :2].to_numpy()

nPairs = 50
G = nx.Graph()
G_monomer = nx.Graph()
G_interface = nx.Graph()
G.add_edges_from(all_residues[:nPairs])
G_monomer.add_edges_from(monomer_residues[:nPairs])
G_interface.add_edges_from(interface_residues[:nPairs])
# nx.draw_networkx(G, pos=nx.planar_layout(G))
# nx.draw_networkx(G_interface, pos=nx.planar_layout(G_interface), node_color='red')

degree_sequencem = sorted([dm for n, dm in G_monomer.degree()], reverse=True)
degree_sequence = sorted([d for n, d in G_interface.degree()], reverse=True)
degreeCountm = collections.Counter(degree_sequencem)
degreeCount = collections.Counter(degree_sequence)
degm, cntm = zip(*degreeCountm.items())
deg, cnt = zip(*degreeCount.items())

# fig, ax = plt.subplots()
# plt.bar(degm, cntm, width=0.8)
# plt.bar(deg, cnt, width=0.8, alpha=0.5, edgecolor='black')
# ax.set_ylim(0,100)
# ax.set_xticks([d + 0.4 for d in deg])
# ax.set_xticklabels(deg)
# plt.axes([0.4, 0.4, 0.5, 0.5])
Gcc = G_interface.subgraph(sorted(nx.connected_components(G_interface), key=len, reverse=True)[0])
posm = nx.spring_layout(G_monomer, seed=np.random.seed(1))
posi = nx.spring_layout(G_interface, seed=np.random.seed(1))
plt.axis('off')
nx.draw_networkx_nodes(G_interface, posi, node_size=20, with_labels=True, node_color="red")
nx.draw_networkx_edges(G_interface, posi, alpha=0.4, edge_color='red')
nx.draw_networkx_nodes(G_monomer, posm, node_size=20, with_labels=True)
nx.draw_networkx_edges(G_monomer, posm, alpha=0.4)
plt.show()