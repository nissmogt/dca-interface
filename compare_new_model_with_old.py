import pandas as pd
from distance_dca import distance_dca, read_dca_distance_matrix
from distance_pdb import plot_cm
from read_db import get_lengths
from output_pipeline import (
    pipeline_distance_matrix, pipeline_mapping, pipeline_interface,
    pipeline_dca_distance, pipeline_fni_scores
)


def compare_new_scrambled_with_paired(msa_name, n_pairs, cutoff):
    """
    PDB protocols:
                    1. Calculate or read distance matrix, plot it if needed (set plot=True).
                    2. Output list of monomer and interface distance matrix Dataframes
                       and chain lengths for each monomer in pdb.
    DCA protocols:
                    1. Get uniprot lengths for each chain in msa (NOTE: chains may not match msa
                       if columns were removed from msa as a result of post-processing.)
                    2. Load Matlab matrix file, calculate FNi scores, output to a Dataframe, rank by top score,
                       and remove i,j pairs with sequence distance < 5.
                    3. Map DCA indices to PDB-distance-matrix indices
                    4. Load vanilla plmDCA predictions and distances
    Plots:
                    1. PDB distance matrix (optional - see PDB protocol)
                    2. PDB distance matrix, DCAi and vanilla DCA top [n_pairs]

    :param msa_name:
    :param n_pairs:
    :param cutoff:
    :return:
    """
    # PDB PROTOCOL ====================================================================================================
    atom = 'aa'  # cutoff type aa: all-atom, ca: c-alpha
    pdb_df_list, chain_lengths = pipeline_distance_matrix(msa_name, cutoff_type=atom, cutoff=cutoff,
                                                          read=True, plot=False)
    pdb_chain1_length = chain_lengths[0]
    pdb_total_length = sum(chain_lengths)

    # DCA PROTOCOL ====================================================================================================
    uniprot_chain_lengths = get_lengths(msa_name)
    df_fni_apc = pipeline_fni_scores(msa_name, read=False)
    df_fni_apc_mapped = pipeline_mapping(msa_name, df_fni_apc, uniprot_chain_lengths, read=True)
    df_fni_apc_mapped = df_fni_apc_mapped[:n_pairs]
    df_fni_apc_mapped_dist = distance_dca(df_fni_apc_mapped, msa_name, atom=atom, other_name=True)
    # df_fni_apc_mapped_dist = read_dca_distance_matrix(msa_name, n_pairs, atom=atom, other_name=True)

    # top_pairs = int(np.ceil(3 * pdb_chain1_length / 2.0))
    # print("\n3L_a/2 Pairs: |{}|\tN_pairs: {}".format(top_pairs, n_pairs))
    # df_dca_mapped_inter_dist = pipeline_dca_distance(msa_name, None, atom, read=True, n_pairs=50)

    # PLOTS PROTOCOL ==================================================================================================
    # df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm if you dont want to plot dca

    # plot_cm(pdb_df_list, cutoff, pdb_chain1_length, pdb_total_length, atom=atom,
    #         df_dca=df_dca_mapped_inter_dist[:n_pairs], msa_name=msa_name, other_dca=df_fni_apc_mapped_dist)
    return df_fni_apc_mapped_dist
