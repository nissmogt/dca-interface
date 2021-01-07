def interface_density(msaName, cutoff=6):
    import pandas as pd

    # cutoff = 6.0
    # msaName = "1FM0_E_1FM0_D"
    chains = [msaName.split("_")[1], msaName.split("_")[3]]
    pdbMatrix = "distance_matrix\\"
    pdbFile = "{}heavy_atom_distance_matrix_{}.txt".format(pdbMatrix, msaName)
    dfPDB = pd.read_csv(pdbFile, delimiter='\t', header=0)
    dfInterface = dfPDB[dfPDB["chain_1"] != dfPDB["chain_2"]]
    dfInterface_cutoff = dfInterface[dfInterface["d"] <= cutoff]

    chain1_length = dfPDB[dfPDB["chain_2"] == chains[0]]["j"].max()
    chain2_length = dfPDB[dfPDB["chain_2"] == chains[1]]["j"].max() - chain1_length
    contact_density = len(dfInterface_cutoff)

    output_list = [msaName, contact_density, chain1_length, chain2_length]
    return output_list


def batch_density(cutoff=6):
    import csv
    import glob
    import os
    sysDir = "PDB_benchmark_alignments\\"
    dcaMatrixDir = "scrambled_results\\fni_matrices\\"
    sysList = glob.glob("{}*.fas".format(sysDir))
    density_list = []
    count = 0
    for sys in sysList:
        msaName = os.path.basename(sys).strip(".fas")
        if os.path.exists("{}matrix_FNi_{}.npy".format(dcaMatrixDir, msaName)):
            count += 1
            print("{} {}".format(count, msaName))
            density_list.append(interface_density(msaName, cutoff))

    outfile = "benchmark{}_interface_contact_per_L_cutoff{}.csv".format(count, cutoff)
    with open(outfile, "w", encoding='utf-8', newline='') as output:
        writer = csv.writer(output)
        writer.writerow(["MSA", "Ninterface/L1*L2", "L1", "L2"])
        writer.writerows(density_list)


batch_density(6)