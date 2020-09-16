def calc_min_dist(res_a, res_b):
    import numpy as np
    dist = []
    for atom_i in res_a:
        for atom_j in res_b:
            dist.append(np.linalg.norm(atom_i.get_coord() - atom_j.get_coord()))
    dist_array = np.array(dist)
    mindist = min(dist_array[np.nonzero(dist_array)])
    return mindist


def calc_ca_distance(res_a, res_b):
    """
    Calculates the distance between a pair of CA atoms
    :param res_a: Biopython residue object - residue a
    :param res_b: Biopython residue object - residue b
    :return: Distance between CA atoms
    """
    import numpy as np
    a = res_a["CA"].get_coord()
    b = res_b["CA"].get_coord()
    dist = np.linalg.norm(a - b)
    return dist


def distance_matrix(msa, cutoff, all_atom=False):
    """
    Calculates distance matrix.
    :param cutoff:
    :param msa:
    :param all_atom:
    :return: Three Dataframe objects.
    """
    # msa = "3RPF_A_3RPF_D.fas"
    # cutoff = 12
    from pandas import read_csv
    from itertools import combinations_with_replacement
    from get_residues import get_residues
    import time
    fname = "(pdb_map)"
    pdb_path = "PDB_benchmark_structures\\"
    msa_name = msa.strip(".fas")
    chain_ids = [msa_name.split("_")[1], msa_name.split("_")[3]]
    all_atom = False
    if all_atom:
        filename = "{}heavy_atom_distance_matrix_{}_{}A.txt".format(pdb_path, msa_name, cutoff)
    else:
        filename = "{}ca_distance_matrix_{}_{}A.txt".format(pdb_path, msa_name, cutoff)
    fileout = open(filename, 'w')
    fileout.write("i\tj\td\tchain_1\tchain_2\n")
    # output list of residues from pdb
    residues, chain_lengths = get_residues(msa)
    # make each possible pairs of residues
    pair_list = combinations_with_replacement(range(len(residues)), 2)
    start_time = time.time()
    for i, j in pair_list:
        res_a = residues[i]
        res_b = residues[j]
        # get chain id
        if all_atom:
            # if res_a.get_id()[1] - res_b.get_id()[1] > 4:
            chain_a = res_a.get_parent().id
            chain_b = res_b.get_parent().id
            mindist = calc_min_dist(res_a, res_b)
            if mindist <= cutoff:
                fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i, j, mindist, chain_a, chain_b))
        else:
            if res_a.has_id("CA") and res_b.has_id("CA"):
                chain_a = res_a.get_parent().id
                chain_b = res_b.get_parent().id
                dist = calc_ca_distance(res_a, res_b)
                if cutoff >= dist > 0.0:
                    fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i, j, dist, chain_a, chain_b))
            else:
                print("{} NOTE! Res {} \n\tor {} not calculated! (missing CA)\n".format(fname, res_a.get_full_id(),
                                                                                        res_b.get_full_id()))
    fileout.close()
    print("{}\t -- MAIN LOOP TIME -- {}".format(fname, time.time() - start_time))
    # makes a pandas dataframe
    df_pdb = read_csv(filename, delim_whitespace=True)
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
    return df_pdb, df_mon, df_inter, chain_lengths


def plot_pdb(pdb_df_list, cutoff, length_a, length, title=None, atom="ca", msa_name=None, count=0):
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"
    """
    :param pdb_df_list:
    :param msa_name:
    :param cutoff:
    :param length_a:
    :param length:
    :param title:
    :param atom:
    :param count:
    :return:
    """
    import matplotlib.pylab as plt

    df_mon = pdb_df_list[0]
    df_inter = pdb_df_list[1]

    # Plotting
    fig = plt.figure(count, figsize=(10, 10), dpi=100)
    ax = fig.add_subplot(1, 1, 1)

    # monomer
    ax.scatter('i', 'j', data=df_mon, label='PDB monomer', c='xkcd:navy', marker='s')
    # interface
    ax.scatter('i', 'j', data=df_inter, label='PDB dimer', c='olive', marker='s')

    # Plot dimer separator line
    ax.hlines(length_a, 0, length, linestyles='dashed', alpha=0.8, colors='xkcd:off white')
    ax.vlines(length_a, 0, length, linestyles='dashed', alpha=0.8, colors='xkcd:off white')

    # plot design
    ax.legend(loc='lower right')
    ax.set_facecolor('tan')
    plt.title("PDB: {}".format(msa_name))
    plt.minorticks_on()
    plt.grid(alpha=0.3)
    plt.xlabel("residue i"), plt.ylabel("residue j")
    ax.grid(which='major', alpha=0.4, c='black')
    ax.grid(which='minor', linestyle=':', alpha=0.5, c='gray')
    if img_dir:
        imgname = "{}_{}_len_{}_{}{}.png".format(title, msa_name, length, atom, cutoff)
        # plt.savefig(img_dir + imgname, dpi=1000, bbox_inches='tight')
        # plt.close()
