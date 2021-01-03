import pandas as pd
import numpy as np


def monomer_restraint1(sysName, df, cutoff):
    from mapping_functions import apply_map

    referenceMap = "results\\reference_maps\\ref_map_{}.txt".format(sysName)
    dfMonomer = df[df["chain_1"] == df["chain_2"]]
    dfMonomer = dfMonomer[dfMonomer["d"] > cutoff].reset_index(drop=True)
    monomer_array = dfMonomer.iloc[:, :3].to_numpy()

    map_pdb_dca = pd.read_csv(referenceMap, delimiter="\t", header=0, dtype=str)
    map_pdb_dca = map_pdb_dca.replace("?", np.nan).dropna()  # some pdbs have unknown seq res UNK
    map_to_dca = dict(zip(map_pdb_dca["pdb_i"], map_pdb_dca["dca_i"]))

    mappedArray = apply_map(monomer_array, map_to_dca)
    r = mappedArray[:, :2]
    exclusions = r.astype(dtype='int')


def monomer_restraint(sysName, df, cutoff):
    from read_db import get_lengths
    from get_region import get_dca_indices
    from mapping_functions import apply_map

    ch = get_lengths(sysName)
    _, dca_ch, _ = get_dca_indices(sysName, ch[0])
    print(dca_ch)
    msa_pairs = []
    protein1 = []
    protein2 = []
    for i in range(1, dca_ch[0]):
        for j in range(i+1, dca_ch[0] + 1):
            protein1.append([i, j])
    for i in range(1, dca_ch[1]):
        for j in range(i+1, dca_ch[1] + 1):
            protein2.append([i+dca_ch[0], j+dca_ch[0]])

    msa_pairs = np.array(protein1 + protein2)

    referenceMap = "results\\reference_maps\\ref_map_{}.txt".format(sysName)
    dfMonomer = df[df["chain_1"] == df["chain_2"]]
    dfMonomer = dfMonomer[dfMonomer["d"] <= cutoff].reset_index(drop=True)
    monomer_array = dfMonomer.iloc[:, :3].to_numpy()

    map_pdb_dca = pd.read_csv(referenceMap, delimiter="\t", header=0, dtype=str)
    map_pdb_dca = map_pdb_dca.replace("?", np.nan).dropna()  # some pdbs have unknown seq res UNK
    map_to_dca = dict(zip(map_pdb_dca["pdb_i"], map_pdb_dca["dca_i"]))

    mappedArray = apply_map(monomer_array, map_to_dca)
    r = mappedArray[:, :2]
    r = r.astype(dtype='int')

    msa_row = msa_pairs.view([('', msa_pairs.dtype)] * msa_pairs.shape[1])
    pdb_row = r.view([('', r.dtype)] * r.shape[1])
    exclusions_list = np.setdiff1d(msa_row, pdb_row).view(msa_pairs.dtype).reshape(-1, msa_pairs.shape[-1])
    assert len(msa_pairs) - len(r) == len(exclusions_list)
    return exclusions_list


dimers = ['5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A']
# dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
#       "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
#       '5MU7_B_5MU7_A']
# dimers = "HKRR_C_HKRR_A"
# for i in range(len(dimers)):
i = -1
sysName = dimers[i]
outfile = "nonbonded_restraints\\{}.fas_restraints.txt".format(sysName)
distanceMatrix = "distance_matrix\\heavy_atom_distance_matrix_{}.txt".format(sysName)
df_matrix = pd.read_csv(distanceMatrix, delimiter="\t", header=0)
res = monomer_restraint(sysName, df=df_matrix)
np.savetxt(outfile, res, fmt='%d\t%d')
