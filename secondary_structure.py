import biotite.structure.io as strucio
import biotite.structure as struc
import matplotlib.pyplot as plt
import numpy as np

def pdb2biotite(msaName):
    pdbDir = "PDB_benchmark_structures\\"
    # msaName = "1FM0_E_1FM0_D"
    array = strucio.load_structure("{}.pdb".format(pdbDir + msaName[:4]))
    return array


def calculate_ss(msaName, chain):
    # msaName = "1FM0_E_1FM0_D"
    array = pdb2biotite(msaName)
    array = array[array.hetero == False]    # filters out hetatoms
    # Estimate secondary structure
    if len(chain) > 1:
        sse = []
        for ch in chain:
            sse.append(struc.annotate_sse(array, chain_id=ch))
        return np.append(sse[0], sse[1])
    else:
        return struc.annotate_sse(array, chain_id=chain)


def calculate_dihedral(msaName):
    array = pdb2biotite(msaName)
    # plt.figure(0)
    # plt.plot(phi * 360/(2*np.pi), psi * 360/(2*np.pi),
    #          marker="o", linestyle="None")
    # plt.xlim(-180,180)
    # plt.ylim(-180,180)
    # plt.xlabel("$\phi$")
    # plt.ylabel("$\psi$")
    # plt.close()
    return struc.dihedral_backbone(array)


def _add_sse_to_sasa_file(msaName):
    """
    function to add calculated secondary structure values to sasa file.
    :param msaName:
    :return:
    """
    import os
    import pandas as pd

    # msaName = "1FM0_E_1FM0_D"
    sasa_file = "sasa\\total_{}_freesasa.txt".format(msaName)
    dfSasa = pd.read_csv(sasa_file, delimiter='\t', header=0)

    # Calculate ss
    chains = [msaName.split("_")[1], msaName.split("_")[3]]
    sse_chains = calculate_ss(msaName, chains)

    dfSasa["ss"] = sse_chains
    outfile = "{}_ss.txt".format(sasa_file.strip(".txt"))
    header = 'Residue\tResidue_num\tChain\tTotal\tApolar\tBackbone\tmainchain\tSidechain\tRatio\tss'
    dfSasa.to_csv(outfile, sep='\t', index=False, header=header)


# msaName = "1FM0_E_1FM0_D"
# msaName = "1SXJ_A_1SXJ_G"
msaName = "1WA5_A_1WA5_C"
_add_sse_to_sasa_file(msaName)
