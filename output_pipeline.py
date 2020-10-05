def pipeline_distance_matrix(msa_name, cutoff_type, cutoff, read=False, plot=False):
    """
    Used in calculating distance matrix for given msa input.
    :param msa_name: string - Name for MSA file should follow the following format 'PDBID_CHAIN1_PDBID_CHAIN2'
    :param cutoff_type: string - 'aa' for all atom min distance calculation. 'ca' for c-alpha distance calculation.
    :param cutoff: float - Cutoff defines contact between two inter-chain residues in PDB.
                   Note: 8A is hard-coded for monomer.
    :param read: boolean - If True, reads distance matrix file into pandas Dataframe. (Default: False)
    :param plot: boolean - Plot distance matrix at chosen cutoff using chosen cutoff type. (Default: False)
    :return: List of two pandas Dataframes. Monomer pairs and interface pairs.
    """
    import numpy as np
    import pandas as pd
    from distance_pdb import distance_matrix, read_distance_matrix_file, plot_cm
    if cutoff_type == 'aa':
        all_atom = True
    else:
        all_atom = False

    sifts_table_file = "databases/sifts/pdb_chain_uniprot_plus.csv"
    s = pd.read_csv(sifts_table_file, comment="#")
    pdbid = msa_name[:4].lower()
    chain_1 = msa_name.split("_")[1]
    chain_2 = msa_name.split("_")[3]

    # it's possible to have multiple coord values per chain entry
    pdb_start_chain_1 = s.query("pdb_id == @pdbid and pdb_chain == @chain_1").coord_start.values
    pdb_start_chain_1 = np.array([int(i) for i in pdb_start_chain_1])
    pdb_end_chain_1 = s.query("pdb_id == @pdbid and pdb_chain == @chain_1").coord_end.values
    pdb_end_chain_1 = np.array([int(i) for i in pdb_end_chain_1])

    pdb_start_chain_2 = s.query("pdb_id == @pdbid and pdb_chain == @chain_2").coord_start.values
    pdb_start_chain_2 = np.array([int(i) for i in pdb_start_chain_2])
    pdb_end_chain_2 = s.query("pdb_id == @pdbid and pdb_chain == @chain_2").coord_end.values
    pdb_end_chain_2 = np.array([int(i) for i in pdb_end_chain_2])

    # -- PDB --
    print("\n\t-- |{}| DISTANCE MATRIX CALCULATION at |{}| inter-cutoff: |{}| --".format(pdbid, cutoff_type, cutoff))
    if read:
        df_mon, df_inter, chain_lengths = read_distance_matrix_file(msa_name, all_atom=all_atom)
    else:
        df_pdb, df_mon, df_inter, chain_lengths = distance_matrix(msa_name, all_atom=all_atom)

    df_mon = df_mon[df_mon["d"] <= 8.0]    # hard-coded monomer cutoff
    df_inter = df_inter[df_inter["d"] <= cutoff]

    total_length = sum(chain_lengths)
    print("\t||Chain {}\t||Chain {}\nlengths: {}\t||{}\t||Total length: {}".format(chain_1, chain_2, chain_lengths[0],
                                                                                   chain_lengths[1], total_length))
    # s = SIFTS("databases\\sifts\\pdb_chain_uniprot_plus.csv", "databases\\sifts\\pdb_chain_uniprot_plus.fa")
    # monomer_1 = s.by_pdb_id(pdb_id=pdbid, pdb_chain=chain_1)
    # monomer_2 = s.by_pdb_id(pdb_id=pdbid, pdb_chain=chain_2)
    # mon_map_1 = monomer_1.mapping[0]
    # mon_map_2 = monomer_2.mapping[0]
    pdb_df_list = [df_mon, df_inter]
    if plot:
        df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm
        plot_cm(pdb_df_list, cutoff=cutoff, length_a=chain_lengths[0], length=total_length, atom=cutoff_type,
                df_dca=df_empty, msa_name=msa_name)
    return pdb_df_list, chain_lengths


def pipeline_mapping(msa_name, df_dca, uniprot_lengths, read=False):
    """
    Map DCA indices to PDB-distance-matrix indices
    :param read:
    :param uniprot_lengths:
    :param df_dca:
    :param msa_name:
    :return:
    """
    import pandas as pd
    from parse_fasta import read_msa
    from get_residues import get_residues
    from mapping_functions import align_dca2pdb, apply_map
    from get_region import get_dca_indices
    print("(pipeline mapping)")
    if read:
        import numpy as np
        infile = "results\\reference_maps\\ref_map_{}.txt".format(msa_name.strip(".fas"))
        map_pdb_dca = pd.read_csv(infile, delimiter="\t", header=0)
        # P1: after running all systems remove line 83/4 (was only added to correct col names in align_dca2pdb line 25)
        map_pdb_dca = map_pdb_dca.rename(columns={"dca_res": "dca_i", "dca_i": "dca_res"})
        np.savetxt(infile, map_pdb_dca, header="pdb_i\tpdb_res\tdca_i\tdca_res", fmt="%s\t%s\t%s\t%s", comments='')
        map_pdb_dca = map_pdb_dca.dropna()
        map_to_pdb = dict(zip(map_pdb_dca["dca_i"], map_pdb_dca["pdb_i"]))

    else:
        _, dca_lengths, _ = get_dca_indices(msa_name, uniprot_lengths[0])

        # -- GET MAP FROM MSA TO PDB --
        pdbseq_1, pdbseq_2 = get_residues(msa_name, seq=True)
        pdbseq = [pdbseq_1, pdbseq_2]
        # splits msa sequence based on modified uniprot lengths (removed lowercase)
        msaseq = read_msa(msa_name, split=True, len_a=dca_lengths[0])

        map_to_pdb = align_dca2pdb(msa_name, pdbseq, msaseq)

    print("(map dictionary) {}".format(map_to_pdb))
    mapped_dca_array = apply_map(df_dca.to_numpy(), map_to_pdb)
    df_dca_mapped = pd.DataFrame(mapped_dca_array, columns=['i', 'j', 'score'])

    return df_dca_mapped


def pipeline_interface(df_dca_mapped, msa_name, pdb_chain1_length):
    """

    :param df_dca_mapped: 
    :param msa_name:
    :param pdb_chain1_length:
    :return:
    """
    import numpy as np
    results_dir = "results\\"
    outfile = "{}FN_{}_inter_mapped.txt".format(results_dir, msa_name)
    df_dca_mapped_inter = df_dca_mapped.loc[(df_dca_mapped["i"] < pdb_chain1_length) &
                                            (df_dca_mapped["j"] > pdb_chain1_length)]
    np.savetxt(outfile, df_dca_mapped_inter, header="i\tj\tscore\t", fmt='%d\t%d\t%f', comments='')
    # P2: add interface column to dca dataframe object i.e. "interface": y or n
    return df_dca_mapped_inter


def pipeline_dca_distance(msa_name, dca_df, atom, read=False, n_pairs=None):
    """

    :param n_pairs:
    :param msa_name:
    :param dca_df:
    :param atom:
    :param read:
    :return:
    """
    from distance_dca import distance_dca, read_dca_distance_matrix
    if read:
        df_dca_mapped_inter_dist = read_dca_distance_matrix(msa_name, n_pairs, atom=atom)
    else:
        n_pairs = len(dca_df)
        df_dca_mapped_inter_dist = distance_dca(dca_df, msa_name, atom=atom)

    return df_dca_mapped_inter_dist


def pipeline_fni_scores(msa_name, read=False):
    print("(pipeline_fni_scores)")
    from dca_functions import average_jmatrix, read_matlab_matrix, calculate_fni_score
    if read:
        import pandas as pd
        results_dir = "scrambled_results\\"
        outfile = "{}FNi_apc_{}.txt".format(results_dir, msa_name)
        return pd.read_csv(outfile, delimiter="\t", names=["i", "j", "score"])
    else:
        vanilla_dca_matrix = "results\\matrix_files\\matrix_ising_{}.fas.mat".format(msa_name)
        avg_J = average_jmatrix(msa_name, 5)
        h_paired, J_paired = read_matlab_matrix(vanilla_dca_matrix)
        df_fni, df_fni_apc = calculate_fni_score(msa_name, J_paired, avg_J)
        return df_fni_apc


def pipeline_confusion_matrix(msa_name, dca_df, pdb_df_list,
                              pdb_total_length, atom, read=False):
    """
    Computes confusion matrix for len(dca_df) number of pairs, adds these values to the input df,
    and saves this df to a new file.
    :param msa_name:
    :param dca_df:
    :param pdb_df_list:
    :param pdb_total_length:
    :param atom:
    :param read:
    :return: Input DataFrame with added confusion matrix values.
    """
    import numpy as np
    import pandas as pd
    from analysis_functions import confusion_matrix_list
    from distance_pdb import vectorize_pdb_contacts
    print("Calculating confusion matrix ...")

    n_pairs = len(dca_df)
    filename_df_cfm = "results\\FN_{}_inter_mapped_aa_dist_top{}_tp.txt".format(msa_name, n_pairs)

    pdb_interface_df = pdb_df_list[1]
    pdb_flat_array = vectorize_pdb_contacts(pdb_interface_df, pdb_total_length)

    confusion_list = confusion_matrix_list(pdb_flat_array, dca_df, pdb_total_length)
    df_dca_mapped_inter_dist_tp = dca_df.assign(tp=confusion_list[0], fp=confusion_list[1],
                                                fn=confusion_list[2], tn=confusion_list[3])
    if read:
        df_dca_mapped_inter_dist_tp = pd.read_csv(filename_df_cfm, delimiter="\t")
    if atom == "aa":
        np.savetxt(filename_df_cfm, df_dca_mapped_inter_dist_tp,
                   header="i\tj\tscore\tdist\tchain_1\tchain_2\tresnames\tatom_id\ttp\tfp\tfn\ttn",
                   fmt='%d\t%d\t%f\t%f\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f', comments='')
    else:
        np.savetxt(filename_df_cfm, df_dca_mapped_inter_dist_tp,
                   header="i\tj\tscore\tdist\tchain_1\tchain_2\tresnames\ttp\tfp\tfn\ttn",
                   fmt='%d\t%d\t%f\t%f\t%s\t%s\t%s\t%f\t%f\t%f\t%f', comments='')

    return df_dca_mapped_inter_dist_tp
