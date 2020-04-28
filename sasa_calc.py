def delta_sasa(pdbfile_1, pdbfile_2):
    import numpy as np
    sasa_res, ratio_1 = np.genfromtxt(pdbfile_1, unpack=True, usecols=(1, 6), skip_header=1)
    chainb_start = np.where(sasa_res == 1.0)[0][-1]
    sasa_res[chainb_start:] = sasa_res[chainb_start:] + sasa_res[chainb_start - 1]

    sasa_res, ratio_2 = np.genfromtxt(pdbfile_2, unpack=True, usecols=(1, 6), skip_header=1)
    chainb_start = np.where(sasa_res == 1.0)[0][-1]
    sasa_res[chainb_start:] = sasa_res[chainb_start:] + sasa_res[chainb_start - 1]

    delta_sasa = np.abs(ratio_1 - ratio_2)
    return dict(zip(sasa_res, delta_sasa))


def process_sasa(sasa_file, *args):
    """
    Processes GetArea Server file and returns dictionary of SASA ratio and total for a given residue
    :param sasa_file: GetArea Server output with footer removed
    :param args: 'dimer' if sasa_file is a dimer
    :return: Dictionary of residue IDs and sasa ratio, and un-normalized sasa
    :bugs: No bugs as of 4/14/2020
    """
    import numpy as np
    sasa_res, total_sasa, ratio = np.genfromtxt(sasa_file, unpack=True, usecols=(1, 2, 6), skip_header=1)
    # Finds where the next chain begins if a dimer
    if 'dimer':
        print("Processing dimer...")
        chain_b_start = np.where(sasa_res == 1.0)[0][-1]
        sasa_res[chain_b_start:] = sasa_res[chain_b_start:] + sasa_res[chain_b_start - 1]
    return dict(zip(sasa_res, ratio)), dict(zip(sasa_res, total_sasa))


def sasa_filter(contact_file, number_of_contacts, sasa_file, sasa_cutoff):
    """
    SASA cutoff filter function that takes inputs and uses the cutoff
    to filter out DCA pairs that have a SASA ratio < cutoff.

    :param contact_file: file in i,j,k column format
    :param number_of_contacts: Integer - number of DCA pairs to include
    :param sasa_file: GetArea Server output with footer removed
    :param sasa_cutoff: Cutoff for SASA ratio or total
    :return: Numpy array of SASA-filtered DCA pairs
    """
    import numpy as np
    res_i, res_j, dca_score = np.loadtxt(contact_file, unpack=True, usecols=(0, 1, 2))
    sasa_ratio, sasa_total = process_sasa(sasa_file, 'dimer')

    dca_filtered_array = []
    for i in range(number_of_contacts):
        try:
            sasa_res_i = sasa_ratio[int(res_i[i])]
            sasa_res_j = sasa_ratio[int(res_j[i])]
            if sasa_res_i >= sasa_cutoff and sasa_res_j >= sasa_cutoff:
                dca_filtered_array.append([res_i[i], res_j[i], dca_score[i],
                                           sasa_ratio[res_i[i]], sasa_ratio[res_j[i]]])
        except KeyError as error:
            print("ERROR: check res: %s" % error.args[0])

    print("Number of DCA pairs above SASA cutoff of %s: %s%%\n" %
          (sasa_cutoff, (len(dca_filtered_array) / float(number_of_contacts) * 100)))
    return np.array(dca_filtered_array)


def sasa_affect_plot(dca_file, number_of_contacts, sasa_file, threshold_list):
    """

    :param dca_file: DCA contact file format i, j, k
    :param number_of_contacts: Number of contacts to use (Integer)
    :param sasa_file: GetArea Server file with footer removed
    :param threshold_list: A list of SASA thresholds to try
    :return: A plot of Total DCA pairs vs SASA-filtered pairs
    """
    import numpy as np
    import matplotlib.pyplot as plt

    filtered_list = np.zeros((len(threshold_list), number_of_contacts))
    thresh_index = 0
    for i in threshold_list:
        for j in range(1, number_of_contacts):
            sasa_filtered_dca_array = sasa_filter(dca_file, j, sasa_file, i)
            filtered_list[thresh_index, j] = (len(sasa_filtered_dca_array))
        thresh_index += 1

    plt.figure(1)
    for i in range(thresh_index):
        plt.scatter(range(number_of_contacts), filtered_list[i, :], label='cutoff: %s' % threshold_list[i])
        plt.plot(range(number_of_contacts), filtered_list[i, :])
    plt.xlabel("number of DCA pairs")
    plt.ylabel("number of SASA-filtered pairs")
    plt.legend(loc='best')
    plt.grid(axis='both', c='black', alpha=0.5)
    plt.show()


def calc_sasa(pdbfile):
    import freesasa as fs
    import numpy as np
    structure = fs.Structure(pdbfile)
    # result = fs.calc(structure)
    result = fs.calc(structure, fs.Parameters({'algorithm': fs.LeeRichards, 'n-slices': 150}))
    area_class = fs.classifyResults(result, structure)

    out = []
    res_name = []
    sasa = []
    res_sasa = []
    atomcount = 0
    for i in range(structure.nAtoms()):
        atom_name = structure.atomName(i)
        res_name.append(structure.residueName(i))

        if i != 0:
            if res_name[i] == res_name[i - 1]:
                atomcount += 1
            else:
                res_sasa.append(sum(sasa))
                #             print("Total SASA: %f for %d atoms for %s" % (sasa_avg, atomcount, res_name[i-1]))
                sasa_avg = 0
                atomcount = 0
                sasa = []
        sasa.append(result.atomArea(i))
    #     print(structure.chainLabel(i), atom_name)
    res_sasa_array = np.array(res_sasa)

    print("Total : %.2f A2" % result.totalArea())
    for key in area_class:
        print(key, ": %.2f A2" % area_class[key])

    return res_sasa_array


def tpr_vs_sasa(pdb_file, dca_file, n_contacts, sasa_ratio_dict, threshold_sasa, dimer_length, chain, calpha_cutoff):
    from dca_performance import vectorize_pdb_contacts, vectorize_dca_contacts
    from sklearn.metrics import precision_score
    import matplotlib.pyplot as plt
    import numpy as np
    import time
    ds = 10
    tpr = []
    start_time = time.time()
    pdb_flat_matrix = vectorize_pdb_contacts(pdb_file, dimer_length)
    pair_per_threshold = []
    for i in np.arange(0, threshold_sasa, ds):

        # dca_sasa_total = sasa_filter(dca_file, n_contacts, sasa_total_dict, i)
        # pair_per_threshold.append(len(dca_sasa_total)/n_contacts)
        dca_sasa_ratio = sasa_filter(dca_file, n_contacts, sasa_ratio_dict, i)
        pair_per_threshold.append(len(dca_sasa_ratio) / n_contacts)

        # Writes a VMD tcl script that draws the DCA-SASA-filtered pairs
        # draw_dca_sasa(dca_sasa_total, n_contacts, length_a, chain)

        # df_dca = pd.DataFrame(dca_sasa_ratio, columns=index)
        # plot_dca_sasa_filter(pdb_file, df_dca, chain, 10, cutoff)
        # tpr_list = loop_tpr(pdb_file, dca_sasa_total, n_contacts, dimer_length, chain, cutoff)
        if len(dca_sasa_ratio) > 0:
            transpose_array = dca_sasa_ratio.transpose()[:2]
            dca_ij_array = transpose_array.transpose()
            dca_flat_matrix = vectorize_dca_contacts(dca_ij_array, dimer_length, n_contacts)

            tpr.append(precision_score(pdb_flat_matrix, dca_flat_matrix, zero_division=1))
        else:
            tpr.append(0)

    print("FOR LOOP: --- %s seconds ---" % (time.time() - start_time))

    plt.title("True positive rate wrt SASA cutoff")
    plt.xlabel("SASA cutoff")
    plt.ylabel("TPR")
    x_range = np.arange(0, threshold_sasa, ds)
    plt.plot(x_range, pair_per_threshold, label='Fraction of DCA pairs')
    plt.plot(x_range, tpr, c='black')
    plt.scatter(x_range, tpr, c='xkcd:red', edgecolors='black', label='TPR')
    plt.grid(which='both', alpha=0.5)
    plt.legend(loc='best')
    plt.show()


def tpr_dca_sasa(pdb_flat_matrix, dca_file, number_of_contacts, sasa_file, threshold_sasa, dimer_length, chain, calpha_cutoff):
    from dca_performance import tpr_top_pairs
    import time
    import numpy as np
    import matplotlib.pyplot as plt
    tpr = []
    start_time = time.time()
    for i in range(len(threshold_sasa)):
        dca_sasa_ratio = sasa_filter(dca_file, number_of_contacts, sasa_file, threshold_sasa[i])

        if len(dca_sasa_ratio) > 0:
            transpose_array = dca_sasa_ratio.transpose()[:2]
            new_dca_array = transpose_array.transpose()
            print('len dca_sasa_ratio: ', len(dca_sasa_ratio))

            start_time2 = time.time()
            tpr.append(tpr_top_pairs(pdb_flat_matrix, new_dca_array, number_of_contacts,
                                     dimer_length, chain, calpha_cutoff))
            print("TPR LOOP: --- %s seconds ---" % (time.time() - start_time2))
        else:
            tpr.append(0)

    print("\tTotal Time: --- %s seconds ---" % (time.time() - start_time))
    plt.figure(0)
    for i in range(len(tpr)):
        top_pairs = range(1, len(tpr[i]) + 1)
        # multiply i and ds to get correct label number
        plt.plot(top_pairs, tpr[i], lw=2, label=("rSASA cutoff: %s" % (threshold_sasa[i])))
        plt.scatter(top_pairs, tpr[i], edgecolors='black')
    plt.xlabel("number of pre-filtered DCA pairs")
    plt.ylabel("True positive rate")
    plt.legend(loc='best')
    plt.grid(axis='both', c='black', alpha=0.5)
    plt.show()
    return tpr


def draw_sasa_res(sasa_file, cutoff):
    """
    Draw residues with SASA >= cutoff in TCL script.
    :param sasa_file:   GetArea Server SASA file 
    :param cutoff: SASA ratio threshold
    :return: TCL script that draws spheres at residue location
    Bugs:
        Non identified as of 4/14/2020
    """
    import numpy as np
    import os
    output_dir = sasa_file.split(os.path.basename(sasa_file))[0]
    res_num, ratio = np.loadtxt(sasa_file, usecols=(1, 6), unpack=True, skiprows=1)
    atomselect = 0
    output_filename = 'draw_sasa_' + str(cutoff) + '.tcl'
    target = open(output_dir + output_filename, 'w')
    target.write("color Display Background gray\n")
    target.write("display projection Orthographic\n")
    target.write("axes location Off\n")
    res_count = 0
    chain = 'C'
    for i in res_num:
        i = int(i)  # needed to make res_num elements callable by index
        # Switches to next chain when residue numbering resets in sasa file
        if res_count > 1 and i == 1:
            chain = 'D'
        # draw the residues that satisfy the cutoff condition
        if cutoff >= ratio[res_count] > 0:
            target.write("set sel%d%s [atomselect top \"chain %s and resid %d\"]\n" % (i, chain, chain, i))
            target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
            atomselect += 1
            if chain == 'D':
                target.write("draw color cyan\n")
            else:
                target.write("draw color green3\n")
            target.write("draw sphere $pos1 radius 2\n")
        res_count += 1
    target.write("mol modselect 0 top \"all\"\n")
    target.write("mol modstyle 0 top newcartoon\n")
    target.write("mol modcolor 0 top chain\n")
    target.close()
    print("Wrote tcl file: %s" % output_dir + output_filename)


def draw_sasa_res_thresh(sasa_file):
    import numpy as np
    resnum, ratio = np.loadtxt(sasa_file, usecols=(1, 6), unpack=True, skiprows=1)
    atomselect = 0
    target = open('files\\XDHB_XDHC\\sasa\\draw_sasa_thresh.tcl', 'w')
    # target.write("color Display Background white\n")
    target.write("display projection Orthographic\n")
    target.write("axes location Off\n")
    # target.write("display depthcue off\n")
    for i in range(len(ratio)):
        target.write("set sel%d%s [atomselect top \"resid %d and chain %s\"]\n"
                     % (i + 1, 'D', i + 1, 'D'))

        target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
        atomselect = atomselect + 1

        si = ratio[i]
        if (si == 0):
            target.write("draw color red\n")
        elif (si > 0 and si <= 5):
            target.write("draw color white\n")
        elif (si > 5 and si <= 20):
            target.write("draw color yellow\n")
        elif (si > 20):
            target.write("draw color green\n")
        target.write("draw sphere $pos1 radius 2\n")

    target.write("mol modselect 0 top \"all\"\n")
    target.write("mol modstyle 0 top newcartoon\n")
    target.write("mol modcolor 0 top chain\n")
    target.close()


def draw_dca_sasa(dca_file, number_of_contacts, length_protein_a, chain, sasa_file, sasa_cutoff):
    """
    :param dca_file: Three-column file of (i,j,score) fashion
    :param number_of_contacts: Number of top DCA pairs to include
    :param length_protein_a: Integer value of the length of first monomer in dimer
    :param chain: List of chain values. (e.g. ['chain_id_a', 'chain_id_b'])
    :param sasa_file: A GetArea server file without footer
    :param sasa_cutoff: Float value
    :return: tcl file

    :bugs: No bugs found as of 4/15/2020
    """
    import numpy as np
    import os
    dca_array = sasa_filter(dca_file, number_of_contacts, sasa_file, sasa_cutoff)
    atomselect = 0
    output_dir = dca_file.split(os.path.basename(dca_file))[0]
    output_file = 'draw_sasa_gt' + str(sasa_cutoff) + '_dca_contacts.tcl'
    target = open(output_dir + output_file, 'w')
    target.write("color Display Background gray\n")
    target.write("axes location Off\n")
    target.write("display rendermode GLSL\n")
    target.write("display ambientocclusion on\n")
    target.write("display depthcue off\n")
    target.write("mol representation VDW\n")
    target.write("mol color restype\n")
    target.write("mol material HardPlastic\n")
    for pair in dca_array[:number_of_contacts]:
        i = int(pair[0])
        j = int(pair[1]) - length_protein_a
        target.write("mol selection {resid %d and chain %s}\n" % (i, chain[0]))
        target.write("mol addrep top\n")
        target.write("mol selection {resid %d and chain %s}\n" % (j, chain[1]))
        target.write("mol addrep top\n")

        # used for drawing bonds between CA atoms of residue pairs
        # P3: Find a way to change CA to heavy atoms
        target.write("set sel%d%s [atomselect top \"name CA and resid %d and chain %s\"]\n"
                     % (i, chain[0], i, chain[0]))
        target.write("set sel%d%s [atomselect top \"name CA and resid %d and chain %s\"]\n"
                     % (j, chain[1], j, chain[1]))
        target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
        atomselect = atomselect + 1
        target.write("lassign [atomselect%d get {x y z}] pos2\n" % atomselect)
        atomselect = atomselect + 1

        target.write("draw color red2\n")
        target.write("draw line $pos1 $pos2 style dashed width 4\n")
    target.write("mol modselect 0 top \"all\"\n")
    target.write("mol modmaterial 0 top Glass3\n")
    target.write("mol modstyle 0 top newcartoon\n")
    target.write("mol modcolor 0 top chain\n")
    target.close()
