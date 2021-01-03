import numpy as np
import pandas as pd


def pipeline_pdb_distance_matrix(msa_name, cutoff_type, cutoff, read=False, plot=False):
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
    from distance_pdb import distance_matrix, import_pdb_distance_matrix, plot_cm
    if cutoff_type == 'aa':
        all_atom = True
    else:
        all_atom = False
    chain1 = msa_name.split("_")[1]
    chain2 = msa_name.split("_")[-1]
    # -- PDB --
    print("\n\t-- |{}| DISTANCE MATRIX CALCULATION at |{}| inter-cutoff: |{}| --".format(msa_name, cutoff_type, cutoff))
    if read:
        df_pdb = import_pdb_distance_matrix(msa_name, all_atom=all_atom)
    else:
        df_pdb = distance_matrix(msa_name, all_atom=all_atom)

    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
    total_length = max(df_mon[df_mon["chain_2"] == chain2].j)
    chain1_length = max(df_mon[df_mon["chain_1"] == chain1].j)
    chain2_length = total_length - chain1_length

    df_mon = df_mon[df_mon["d"] <= 8.0]    # hard-coded monomer cutoff
    df_inter = df_inter[df_inter["d"] <= cutoff]

    # total_length = sum(chain_lengths)
    print("\t||Chain {}\t||Chain {}\nlengths: {}\t||{}\t||Total length: {}".format(chain1, chain2, chain1_length,
                                                                                   chain2_length, total_length))
    pdb_df_list = [df_pdb, df_mon, df_inter]
    if plot:
        df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm
        plot_cm(pdb_df_list, cutoff=cutoff, length_a=chain1_length, length=total_length, atom=cutoff_type,
                df_dca=df_empty, msa_name=msa_name, other_dca=df_empty)

    return pdb_df_list, [chain1_length, chain2_length, total_length]


def pipeline_mapping(msa_name, df_dca, read=False):
    """
    Map DCA indices to PDB-distance-matrix indices
    :param read:
    :param df_dca:
    :param msa_name:
    :return:
    """
    from msa_functions import read_first_sequence_in_msa
    from read_db import get_lengths
    from get_residues import get_residues
    from get_region import get_dca_indices
    from mapping_functions import align_dca2pdb, apply_map
    print("(pipeline mapping)")
    if read:
        infile = "reference_maps\\ref_map_{}.txt".format(msa_name.strip(".fas"))
        map_pdb_dca = pd.read_csv(infile, delimiter="\t", header=0, dtype=str)
        # map_pdb_dca = map_pdb_dca.replace("?", np.nan).dropna()    # some pdbs have unknown seq res UNK
        # map_pdb_dca = pd.read_csv(infile, delimiter="\t", header=0, dtype=str, usecols=(0,1)) # used for HK-RR
        # map_pdb_dca["#HMM"] = map_pdb_dca["#HMM"].astype(int) + 1 # used for HK-RR
        # map_to_pdb = dict(zip(map_pdb_dca["#HMM"], map_pdb_dca["col"])) # used for HK-RR
        map_pdb_dca = map_pdb_dca.replace("X", np.nan).dropna()    # some pdbs have unknown seq res UNK
        map_to_pdb = dict(zip(map_pdb_dca["dca_i"], map_pdb_dca["pdb_i"]))

    else:
        uniprot_lengths = get_lengths(msa_name)
        _, dca_lengths, _ = get_dca_indices(msa_name, uniprot_lengths[0])

        # -- GET MAP FROM MSA TO PDB --
        pdbseq_1, pdbseq_2 = get_residues(msa_name, seq=True)
        pdbseq = [pdbseq_1, pdbseq_2]
        # splits msa sequence based on modified uniprot lengths (removed lowercase)
        msaseq = read_first_sequence_in_msa(msa_name, split=True, len_a=dca_lengths[0])

        map_to_pdb = align_dca2pdb(msa_name, pdbseq, msaseq)

    # print("(map dictionary) {}".format(map_to_pdb))
    mapped_dca_array = apply_map(df_dca.to_numpy(), map_to_pdb)
    df_dca_mapped = pd.DataFrame(mapped_dca_array, columns=['i', 'j', 'fn_apc', 'fn', 'ui', 'uj'])
    df_dca_mapped['i'] = df_dca_mapped['i'].astype(int)
    df_dca_mapped['j'] = df_dca_mapped['j'].astype(int)
    df_dca_mapped['ui'] = df_dca_mapped['ui'].astype(int)
    df_dca_mapped['uj'] = df_dca_mapped['uj'].astype(int)

    return df_dca_mapped


def pipeline_interface(df_dca_mapped, msa_name, pdb_chain1_length):
    """

    :param df_dca_mapped: 
    :param msa_name:
    :param pdb_chain1_length:
    :return:
    """
    df_dca_mapped_inter = df_dca_mapped.loc[(df_dca_mapped["i"] < pdb_chain1_length) &
                                            (df_dca_mapped["j"] > pdb_chain1_length)]
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
        df_dca_dist = read_dca_distance_matrix(msa_name, n_pairs, atom=atom)
    else:
        df_dca_dist = distance_dca(dca_df, msa_name, atom=atom)

    return df_dca_dist


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
