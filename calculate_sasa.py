import time
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def compare_to_getarea(pdbid, chain='A'):
    getarea_file = "sasa\\PDB\\{}_{}_SASA.csv".format(pdbid.upper(), chain)
    df_getarea = pd.read_csv(getarea_file, delimiter=",", usecols=(1, 6))
    df_getarea = df_getarea.dropna().reset_index(drop=True)

    freesasa_file = "sasa\\PDB\\{}_freesasa.csv".format(pdbid.upper())
    df_freesasa = pd.read_csv(freesasa_file, delimiter=",", usecols=(1, 6))

    output_dir = os.path.dirname(getarea_file)
    a = np.array([df_getarea["Ratio"].to_numpy(), df_freesasa["Ratio"].to_numpy()])
    a = a.transpose()
    np.savetxt("{}\\{}_compare.csv".format(output_dir, pdbid), a, fmt="%f",
               header="GetArea, FreeSasa", comments='', delimiter=',')

    maxGA = df_getarea["Ratio"].max()
    maxFS = df_freesasa["Ratio"].max()
    max_range = int(max(maxGA, maxFS)) + 1
    plt.figure(0)
    plt.plot(range(max_range), range(max_range), color='black')
    plt.scatter(df_getarea["Ratio"], df_freesasa["Ratio"], s=30, color="xkcd:blue")
    plt.xlabel("GetArea Ratio")
    plt.ylabel("freesasa Ratio SR-np:1k (uses GetArea randomcoil)")
    plt.savefig("{}\\{}_compare_to_getarea.png".format(output_dir, pdbid), dpi=900, bbox_inches="tight")
    plt.close()
    # plt.show()
    return df_freesasa, df_getarea


def _test():
    msaName = "1FM0_E_1FM0_D"
    sasa_file = "sasa\\total_{}_freesasa.txt".format(msaName)
    pdbMatrix = "distance_matrix\\"
    pdbFile = "{}heavy_atom_distance_matrix_{}.txt".format(pdbMatrix, msaName)
    # pdbFile = "PDB_benchmark_structures\\1em8.1891.pdb.contacts"
    dfPDB = pd.read_csv(pdbFile, delimiter='\t', header=0)
    # dfPDB = pd.read_csv(pdbFile, delim_whitespace=True, header=0, names=["chain_2", "j", "chain_1", "i"])
    # Load ratios into pandas dataframe
    dfSasa = pd.read_csv(sasa_file, delimiter='\t', header=0)
    # dfSasa = calculate_sasa(pdbfile=pdbFile, chain='D', multichain=False)
    # dfSasa = dfSasa.dropna()
    dfPDB = get_system_sasa(dfPDB, dfSasa)

    assert len(dfPDB) > 0, "Length of dfPDB is 0! Check get_system_sasa"
    outfile = "sasa\\test_{}_sasa.txt".format(os.path.basename(pdbFile).strip(".txt"))
    header = "i\tj\td\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id\tratio_i\tratio_j"
    dfPDB.to_csv(outfile, sep='\t', index=False, header=header)

    return dfPDB


def _test_classify(msaName, width=1, cutoff=6.0):
    file_sasa = "sasa\\test_heavy_atom_distance_matrix_{}_sasa.txt".format(msaName)
    # pattern = "sasa\\heavy_atom_distance_matrix_*_sasa.txt"
    # pdbSasaList = glob.glob("{}".format(pattern))
    binned_ratios = np.arange(0, 102, width)
    nBins = len(binned_ratios)
    TPR = np.zeros(nBins)
    pair_count = np.zeros(nBins)

    sysInterface_list = []
    df_systems = pd.read_csv(file_sasa, delim_whitespace=True, header=0,
                             names=["chain_2", "j", "chain_1", "i", "ratio_j", "ratio_i"])
    # df_systems = pd.read_csv(file_sasa, delimiter='\t', header=0)
    df_pdb_sasa = df_systems.dropna()
    df_interface_sasa = df_pdb_sasa[df_pdb_sasa["chain_1"] != df_pdb_sasa["chain_2"]].reset_index(drop=True)
    # df_interface_sasa = df_interface_sasa[df_interface_sasa["d"] <= cutoff].reset_index(drop=True)
    total_pairs = len(df_pdb_sasa)
    # TODO: save list to file with sysname, fraction_interface
    for idx in range(len(df_interface_sasa)):
        ratio_i = df_interface_sasa["ratio_i"][idx]
        ratio_j = df_interface_sasa["ratio_j"][idx]
        min_ratio = min(ratio_i, ratio_j)
        # distance = df_interface_sasa["d"][idx]
        for binIndex in range(nBins - 1):
            if binned_ratios[binIndex] <= min_ratio < binned_ratios[binIndex + 1]:
                pair_count[binIndex] += 1  # count all pairs that lie in this bin
                # if distance <= cutoff:
                TPR[binIndex] += 1  # add count to tpr bin
    return TPR, pair_count


def make_bins(nSystems, systemsList, dx, cutoff=6.0, types='interface'):
    """
    Define max and min of bin from all system score extrema

    :param cutoff:
    :param nSystems:
    :param systemsList:
    :param dx: string; Bin width
    :param types:  NOTE: make_bins has an option called 'types' which is used to distinguish
    between 'interface' and 'monomer' for a system. I.e., it makes bins depending on types.
    :return: np.array binned scores from min to max in steps of dx
    """
    assert len(systemsList) > 0, "Check systemsList! Length is zero."
    maxBinvalue = np.zeros(nSystems)
    minBinvalue = np.zeros(nSystems)
    for i in range(nSystems):
        df_pdb_sasa = pd.read_csv(systemsList[i], delimiter='\t', header=0)
        if types == 'interface':
            df_type = df_pdb_sasa[df_pdb_sasa["chain_1"] != df_pdb_sasa["chain_2"]]
        else:
            df_type = df_pdb_sasa[df_pdb_sasa["chain_1"] == df_pdb_sasa["chain_2"]]
        df_type = df_type[df_type["d"] <= cutoff].reset_index(drop=True)
        maxBinvalue[i] = max(df_type["ratio_i"].max(), df_type["ratio_j"].max())
        minBinvalue[i] = min(df_type["ratio_i"].min(), df_type["ratio_j"].min())
    maxBinvalue = max(maxBinvalue)
    minBinvalue = min(minBinvalue)
    binned_scores = np.arange(0, maxBinvalue+2+dx, dx)
    print("min bin: {}, max bin: {}, width: {}".format(minBinvalue, maxBinvalue, dx))
    return binned_scores


def get_system_sasa(df_system, df_sasa):
    """
    This function adds SASA ratios to the corresponding residue pairs in a pdb distance file. Used in 'calculate_sasa'.
    :param df_system: Pandas Dataframe composed of residue pairs and other information about the PDB.
    :param df_sasa: Pandas Dataframe composed of residue number and SASA values.
    :return: Pandas Dataframe of df_system with two added columns of SASA ratios for each residue in pairs.
    """
    sasa_system = []
    assert len(df_system) > 0, "(get_system_sasa) Check df_system, it's empty!"
    for i in range(len(df_system)):
        resi = df_system["i"][i]
        resj = df_system["j"][i]
        if resi <= len(df_sasa) and resj <= len(df_sasa):
            sasa_system.append([df_sasa["Ratio"][resi - 1],
                                df_sasa["Ratio"][resj - 1]])
        else:
            sasa_system.append([np.nan, np.nan])
    df_system["ratio_i"] = np.transpose(sasa_system)[0]
    df_system["ratio_j"] = np.transpose(sasa_system)[1]
    return df_system


def calculate_sasa(pdbfile, chain, multichain=True, relative_type='sidechain'):
    """

    :param pdbfile: String of PDB file name.
    :param chain: String or List of chain identifiers.
    :param multichain: Boolean. True to separate chains. This allows SASA calculation for a single unattached monomer.
    False if you want to calculate SASA for the structure 'as-is'.
    :return: Pandas Dataframe of residue number, types, and sasa values as columns.
    """
    import freesasa as fs
    dict_max_acc = {
        # Miller max acc: Miller et al. 1987 https://doi.org/10.1016/0022-2836(87)90038-6
        # Wilke: Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635
        # Sander: Sander & Rost 1994 https://doi.org/10.1002/prot.340200303
        "Miller": {
            "ALA": 113.0,
            "ARG": 241.0,
            "ASN": 158.0,
            "ASP": 151.0,
            "CYS": 140.0,
            "GLN": 189.0,
            "GLU": 183.0,
            "GLY": 85.0,
            "HIS": 194.0,
            "ILE": 182.0,
            "LEU": 180.0,
            "LYS": 211.0,
            "MET": 204.0,
            "PHE": 218.0,
            "PRO": 143.0,
            "SER": 122.0,
            "THR": 146.0,
            "TRP": 259.0,
            "TYR": 229.0,
            "VAL": 160.0,
        },
        "Wilke": {
            "ALA": 129.0,
            "ARG": 274.0,
            "ASN": 195.0,
            "ASP": 193.0,
            "CYS": 167.0,
            "GLN": 225.0,
            "GLU": 223.0,
            "GLY": 104.0,
            "HIS": 224.0,
            "ILE": 197.0,
            "LEU": 201.0,
            "LYS": 236.0,
            "MET": 224.0,
            "PHE": 240.0,
            "PRO": 159.0,
            "SER": 155.0,
            "THR": 172.0,
            "TRP": 285.0,
            "TYR": 263.0,
            "VAL": 174.0,
            "MSE": 224.0,
            "SEC": 167.0,
        },
        "Sander": {
            "ALA": 106.0,
            "ARG": 248.0,
            "ASN": 157.0,
            "ASP": 163.0,
            "CYS": 135.0,
            "GLN": 198.0,
            "GLU": 194.0,
            "GLY": 84.0,
            "HIS": 184.0,
            "ILE": 169.0,
            "LEU": 164.0,
            "LYS": 205.0,
            "MET": 188.0,
            "PHE": 197.0,
            "PRO": 136.0,
            "SER": 130.0,
            "THR": 142.0,
            "TRP": 227.0,
            "TYR": 222.0,
            "VAL": 142.0,
        },
    }
    theoreticalMaxASA = dict_max_acc["Wilke"]

    # Calculates SASA for unseparated chains.
    if not multichain:
        structure = fs.Structure(pdbfile)
    else:
        # Separate chains if multichain structure. This allows SASA calculation for a single unattached monomer.
        structures = fs.structureArray(pdbfile, options={"separate-chains": True})
        chains = []
        for c in range(len(structures)):
            chains.append(structures[c].chainLabel(1))
        structure = structures[chains.index(chain)]
        print("using {} separating chains {}".format(chains.index(chain), chains))

    print("Number of atoms of {}: {}".format(pdbfile, structure.nAtoms()))
    result = fs.calc(structure, fs.Parameters({'algorithm': fs.ShrakeRupley, 'n-points': 10000}))
    res = result.residueAreas()
    residue = []
    resnum = []
    total = []
    apolar = []
    mainchain = []
    sidechain = []
    ratio = []

    for idx, v in res[chain].items():
        residue.append(v.residueType)
        resnum.append(v.residueNumber)
        total.append(v.total)
        apolar.append(v.apolar)
        mainchain.append(v.mainChain)
        sidechain.append(v.sideChain)
        if v.residueType == 'GLY':
            ratio.append(100 * v.mainChain / theoreticalMaxASA[v.residueType])
        elif v.residueType not in theoreticalMaxASA.keys():
            possibleSASA = []
            for i, maxSASA in enumerate(theoreticalMaxASA.values()):
                # If the residue is unknown but has a SASA,
                # calculate the rSASA dividing by theoretical maxSASA and then use the average of that value
                possibleSASA.append(100 * v.sideChain / maxSASA)
            ratio.append(np.average(possibleSASA))
        else:
            if relative_type == 'sidechain':
                ratio.append(100 * v.sideChain / theoreticalMaxASA[v.residueType])
            else:
                ratio.append(100 * v.total / theoreticalMaxASA[v.residueType])

        # if v.hasRelativeAreas:
        #     ratio.append(v.relativeSideChain)
        # else:
        #     ratio.append(np.nan)

    df_sasa = pd.DataFrame({'Residue': residue, 'Residue_num': resnum, 'Chain': chain, 'Total': total, 'Apolar': apolar,
                            'Backbone': mainchain, 'Sidechain': sidechain, 'Ratio': ratio})
    area_class = fs.classifyResults(result, structure)
    print("Total : %.2f A2" % result.totalArea())
    for key in area_class:
        print(key, ": %.2f A2" % area_class[key])

    return df_sasa


def batch_calculate(calc_type, result_dir, inputList=None):
    sysDir = "PDB_benchmark_alignments\\"
    pdbDir = "PDB_benchmark_structures\\"
    dcaMatrixDir = "scrambled_results\\fni_matrices\\"
    if inputList:
        sysList = inputList
    else:
        sysList = glob.glob("{}*.fas".format(sysDir))

    count = 0
    start_time = time.time()
    for sysFile in sysList:
        # TODO: I want to be able to add a list of systems to NOT run
        msaName = os.path.basename(sysFile).strip(".fas")
        if len(sysList) > 0:
        # if os.path.exists("{}matrix_FNi_{}.npy".format(dcaMatrixDir, msaName)):
            # if os.path.exists("sasa\\{}_freesasa.txt".format(msaName)):
            print(msaName)
            count += 1
            # sasa_file = "sasa\\total_{}_freesasa.txt".format(msaName)
            sasa_file = "sasa\\sidechain_sasa\\{}_freesasa.txt".format(msaName)

            if calc_type == 'sasa':
                dfSystems = []
                chains = [msaName.split("_")[1], msaName.split("_")[3]]
                for chain in chains:
                    pdbName = "{}{}.pdb".format(pdbDir, msaName[:4])
                    dfSystems.append(calculate_sasa(pdbName, chain, multichain=True))
                df_concatenated_chains = pd.concat(dfSystems)
                header = "Residue\tResidue_num\tChain\tTotal\tApolar\tBackbone\tSidechain\tRatio"
                format_string = "%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f"
                np.savetxt(sasa_file, df_concatenated_chains, header=header, fmt=format_string, comments='')

            elif calc_type == 'pdb':
                pdbMatrix = "distance_matrix\\"
                pdbFile = "{}heavy_atom_distance_matrix_{}.txt".format(pdbMatrix, msaName)
                # pdbFile = "PDB_benchmark_structures\\1em8.1890.pdb.contacts"
                dfPDB = pd.read_csv(pdbFile, delimiter='\t', header=0)
                # Load ratios into pandas dataframe
                dfSasa = pd.read_csv(sasa_file, delimiter='\t', header=0)
                # dfSasa = dfSasa.dropna()
                dfPDB = get_system_sasa(dfPDB, dfSasa)

            elif calc_type == 'dca':
                dist = "{}FN_inter_{}_mapped_aa_dist.txt".format(result_dir, msaName)
                # pdbFile = "PDB_benchmark_structures\\1em8.1890.pdb.contacts"
                dfDCA = pd.read_csv(dist, delimiter='\t', header=0)
                # Load ratios into pandas dataframe
                dfSasa = pd.read_csv(sasa_file, delimiter='\t', header=0)
                # dfSasa = dfSasa.dropna()
                dfDCA = get_system_sasa(dfDCA, dfSasa)
                assert len(dfDCA) > 0, "Length of dfPDB is 0! Check get_system_sasa"
                outfile = "{}_sasa.txt".format(dist.strip(".txt"))
                header = "i\tj\tscore\tdist_aa\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id\tratio_i\tratio_j"
                dfDCA.to_csv(outfile, sep='\t', index=False, header=header, float_format='%.5f')

            elif calc_type == 'draw':
                from draw_sasa import draw_sasa_res
                draw_sasa_res(msaName)

            elif calc_type == 'ss':
                from secondary_structure import _add_sse_to_sasa_file
                _add_sse_to_sasa_file(msaName)


    print("-- Total time to run {} systems: {} --".format(count, time.time() - start_time))
    print("=" * 72)


def classify_sasa(cutoff=6.0, inputList=None, types='interface'):
    """

    :param cutoff:
    :param inputList: Optional A list of msa names
    :return:
    """
    # This part is used if you only want a set list of systems to classify.
    if inputList:
        pdbSasaList = []
        for msa in inputList:
            pattern = "sasa\\sidechain_sasa\\heavy_atom_distance_matrix_{}_sasa.txt".format(msa)
            pdbSasaList.append(pattern)
    # Else we use all systems we've calculated sasa for
    else:
        pattern = "sasa\\sidechain_sasa\\heavy_atom_distance_matrix_*_sasa.txt"
        pdbSasaList = glob.glob("{}".format(pattern))

    # This next part consists of calculating bins for our system/s
    # NOTE: make_bins has an option called 'types' which is used to distinguish
    # between 'interface' and 'monomer' for a system. I.e., it makes bins depending on types.
    width = 1
    print("Making bins")
    binned_ratios = make_bins(len(pdbSasaList), pdbSasaList, width, cutoff, types=types)
    nBins = len(binned_ratios)
    TPR = np.zeros(nBins)
    type_1 = np.zeros(nBins)
    type_2 = np.zeros(nBins)
    type_3 = np.zeros(nBins)
    type_4 = np.zeros(nBins)
    type_5 = np.zeros(nBins)
    type_6 = np.zeros(nBins)

    sysInterface_list = []
    for k, sys in enumerate(pdbSasaList):
        df_systems = pd.read_csv(sys, delimiter='\t', header=0)
        df_pdb_sasa = df_systems.dropna()
        if types == 'interface':
            df_interface_sasa = df_pdb_sasa[df_pdb_sasa["chain_1"] != df_pdb_sasa["chain_2"]]
            df_interface_sasa = df_interface_sasa[df_interface_sasa["d"] <= cutoff].reset_index(drop=True)
            # df_interface_sasa = df_interface_sasa[df_interface_sasa["d"] > 12.0].reset_index(drop=True)
            df_to_classify = df_interface_sasa
            total_pairs = len(df_to_classify)
        else:
            df_pdb_sasa = df_pdb_sasa[df_pdb_sasa["chain_1"] == df_pdb_sasa["chain_2"]]
            df_pdb_sasa = df_pdb_sasa[df_pdb_sasa["d"] <= cutoff].reset_index(drop=True)
            df_to_classify = df_pdb_sasa
            total_pairs = len(df_to_classify)

        sysName = "_".join(pdbSasaList[k].split("_")[4:7])
        print("System: {}\tinterface pairs < {}A: {}/{}".format(sysName, cutoff,
                                                                len(df_to_classify), total_pairs))
        sysInterface_list.append([sysName, len(df_to_classify) / total_pairs])
        # TODO: save list to file with sysname, fraction_interface
        for idx in range(len(df_to_classify)):
            ratio_i = df_to_classify["ratio_i"][idx]
            ratio_j = df_to_classify["ratio_j"][idx]
            min_ratio = min(ratio_i, ratio_j)
            # distance = df_to_classify["d"][idx]
            for binIndex in range(nBins - 1):
                if binned_ratios[binIndex] <= min_ratio < binned_ratios[binIndex + 1]:
                    # if distance <= cutoff:
                    TPR[binIndex] += 1  # add count to tpr bin
                    # Type 1: Exposed+Exposed
                    if ratio_i >= 50 and ratio_j >= 50:
                        type_1[binIndex] += 1
                    # Type 2: Exposed+LessExposed
                    if (ratio_i >= 50 and 20 <= ratio_j < 50) or (ratio_j >= 50 and 20 <= ratio_i < 50):
                        type_2[binIndex] += 1
                    # Type 3: Exposed+Buried
                    if (ratio_i >= 50 and ratio_j < 20) or (ratio_j >= 50 and ratio_i < 20):
                        type_3[binIndex] += 1
                    # Type 4: LessExposed+Buried
                    if (20 <= ratio_i < 50 and ratio_j < 20) or (20 <= ratio_j < 50 and ratio_i < 20):
                        type_4[binIndex] += 1
                    # Type 5: Buried+Buried
                    if ratio_i < 20 and ratio_j < 20:
                        type_5[binIndex] += 1
                    # Type 6: LessExposed+LessExposed
                    if 20 <= ratio_i < 50 and 20 <= ratio_j < 50:
                        type_6[binIndex] += 1

    nSystems = k
    interface_type = [type_1, type_2, type_3, type_4, type_5, type_6]
    print("Total systems: ", nSystems)
    return TPR, interface_type, sysInterface_list


def plot_histograms(tpr, types, cutoff=6.0):
    tpr_norm = tpr / sum(tpr)
    # Plot histogram of number of dca predictions in each bin
    plt.figure(0)
    plt.bar(range(len(tpr)), tpr_norm, width=1, align='edge', color="xkcd:gold", edgecolor="black")
    # plt.vlines(50, ymin=0, ymax=max(tpr_norm), color="black", linestyles="dashed")
    # plt.vlines(20, ymin=0, ymax=max(tpr_norm), color="crimson", linestyles="dashed")
    # plt.yscale(value="log")
    plt.ylabel(
        "Probability of seeing an interface pair with at least one residue at this relative SASA value (at distance <= {}$\AA$)".format(cutoff))
    plt.xlabel("relative SASA")

    plt.figure(1)
    labels = ["Exposed+Exposed", "Exposed+ModeratelyExposed", "Exposed+Buried", "ModeratelyExposed+Buried",
              "Buried+Buried", "ModeratelyExposed+ModeratelyExposed"]
    colors = cm.Set1(np.linspace(0, 1, 10))
    for i in range(len(types)):
        plt.bar(range(len(types[i])), types[i] / sum(tpr), width=1, align='edge', alpha=0.7,
                color=colors[i], edgecolor='black', label=labels[i])
    plt.ylabel(
        "Probability of seeing an interface pair with at least one residue at this relative SASA value (at distance <= {}$\AA$)".format(cutoff))
    plt.xlabel("relative SASA")
    plt.legend()
    plt.show()


# dimers = ["1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B', '5MU7_B_5MU7_A']
# dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "2OXG_Z_2OXG_Y", "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A']
dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
          "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5M72_A_5M72_B', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
          '5MU7_B_5MU7_A']
# r = 'nonbonded_restraints_results\\APC\\20A\\'
r = "vanilla_results//"
# batch_calculate(calc_type='dca', result_dir=r, inputList=dimers)
# for d in range(len(dimers)):
#     msa = dimers[d]
#     data = "results\\FN_{}_inter_mapped_aa_dist_top10000.txt".format(msa)
#     sasaData = "sasa\\sidechain_sasa\\{}_freesasa.txt".format(msa)
#
#     dfDCA = pd.read_csv(data, delimiter='\t', header=0)
#     dfSASA = pd.read_csv(sasaData, delimiter='\t', header=0)
#     df = get_system_sasa(dfDCA, dfSASA)
#     sThreshold = [0, 5, 15, 25, 35, 45, 55]
#     dThreshold = 12
#     dx = 0.025
#     binned_scores = np.arange(0, max(df["score"]), dx)
#     nBins = len(binned_scores)
#     fp = np.zeros((len(sThreshold), nBins))
#     tp = np.zeros((len(sThreshold), nBins))
#     fig, ax = plt.subplots(nrows=len(sThreshold), ncols=2, sharex='all', sharey=False)
#     for si, s in enumerate(sThreshold):
#         for i in range(len(df)):
#             distance = df["dist_aa"][i]
#             score = df["score"][i]
#             ratio_i = df["ratio_i"][i]
#             ratio_j = df["ratio_j"][i]
#             min_ratio = min(ratio_i, ratio_j)
#             for binIndex in range(nBins - 1):
#                 if binned_scores[binIndex] <= score < binned_scores[binIndex + 1]:
#                     if distance > dThreshold and min_ratio > s:
#                         fp[si, binIndex] += 1
#                     elif distance < dThreshold and min_ratio > s:
#                         tp[si, binIndex] += 1
#
#         ax[si, 0].bar(binned_scores, fp[si], width=dx, align='edge', edgecolor="black", color="purple", label='sasa > {}'.format(s))
#         ax[si, 1].bar(binned_scores, tp[si], width=dx, align='edge', edgecolor="black", alpha=0.7, label='sasa > {}'.format(s))
#         ax[si, 0].legend(loc='best')
#         ax[si, 1].legend(loc='best')
#         ax[si, 0].set_xlabel("FN-apc score")
#         ax[si, 0].set_ylabel("FP counts")
#         ax[si, 1].set_xlabel("FN-apc score")
#         ax[si, 1].set_ylabel("TP counts")
#         ax[si, 0].set_yscale("log")
#         ax[si, 1].set_yscale("log")
#         ax[si, 1].set_ylim(0, 100)
#     plt.title(msa)
#     plt.show()


#     if distance > dThreshold:




# batch_calculate(calc_type='pdb')
# batch_calculate(calc_type='sasa', inputList=dimers)
# _test()
# threshold = 6.0
# true_interface, interface_types, slist_ = classify_sasa(cutoff=threshold, inputList=["1FM0_E_1FM0_D"])
# plot_histograms(true_interface, interface_types, cutoff=threshold)
# tpr2, types2, slist2 = classify_sasa2(cutoff=threshold)
# tpr_norm = tpr / sum(tpr)
# tpr_norm2 = tpr2 / sum(tpr2)
# Plot histogram of number of dca predictions in each bin
# plt.figure(0)
# plt.bar(range(len(tpr)), tpr_norm, width=1, align='edge', color="xkcd:gold", edgecolor="black", label='monomer')
# plt.bar(range(len(tpr2)), tpr_norm2, width=1, align='edge', color="xkcd:cyan", edgecolor="black", label='dimer')
# plt.ylabel(
#     "Probability of seeing an interface pair with at least one residue at this relative SASA value (at distance <= {}$\AA$)".format(threshold))
# plt.xlabel("relative SASA")