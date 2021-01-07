import pandas as pd
import numpy as np


def mapIt(array, mapList):
    mappedList = []
    for i in array:
        if str(i) in mapList.keys():
            mappedList.append(int(mapList[str(i)]))

    return np.array(mappedList)


def buried_restraints(msa, buriedCutoff, sasaResidues):
    referenceMap = "results\\reference_maps\\ref_map_{}.txt".format(msa)


    chains = [msa.split('_')[1], msa.split('_')[3]]
    monomer1 = sasaResidues[sasaResidues["Chain"] == chains[0]]
    monomer2 = sasaResidues[sasaResidues["Chain"] == chains[1]]

    filterMonomer1 = monomer1[monomer1["Ratio"] <= buriedCutoff].reset_index()
    filterMonomer2 = monomer2[monomer2["Ratio"] <= buriedCutoff].reset_index()

    array1 = filterMonomer1["index"].to_numpy()
    array2 = filterMonomer2["index"].to_numpy()

    map_pdb_dca = pd.read_csv(referenceMap, delimiter="\t", header=0, dtype=str)
    map_pdb_dca = map_pdb_dca.replace("?", np.nan).dropna()  # some pdbs have unknown seq res UNK
    map_to_dca = dict(zip(map_pdb_dca["pdb_i"], map_pdb_dca["dca_i"]))

    mappedArray1 = mapIt(array1, map_to_dca)
    mappedArray2 = mapIt(array2, map_to_dca)

    r = np.array([[i, j] for i in mappedArray1
                         for j in mappedArray2])

    return r


bc = 5
# dimers = ['5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A']
msaName = '1EM8_D_1EM8_C'
outfile = "sasa_restraints\\{}_{}_buried_restraints.txt".format(msaName, bc)
sasaFile = "sasa\\sidechain_sasa\\{}_freesasa.txt".format(msaName)
resList = pd.read_csv(sasaFile, header=0, delimiter='\t')
#
res = buried_restraints(msaName, buriedCutoff=bc, sasaResidues=resList)
# np.savetxt(outfile, res, fmt='%d\t%d')

