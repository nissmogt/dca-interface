import numpy as np
import pandas as pd
from distance_dca import distance_dca, read_dca_distance_matrix
from distance_pdb import plot_cm
from read_db import get_lengths
from load_dca import load_dca_to_df
from analysis_functions import plot_true_positives
from output_pipeline import (
    pipeline_distance_matrix, pipeline_mapping, pipeline_interface,
    pipeline_dca_distance
)


def run_analysis(msa_name, n_pairs, cutoff, results_dir, count=0):
    """
    PDB protocols:
                    1. Calculate or read distance matrix, plot it if needed (set plot=True).
                    2. Output list of monomer and interface distance matrix Dataframes
                       and chain lengths for each monomer in pdb.
    DCA protocols:
                    1. Get uniprot lengths for each chain in msa (NOTE: chains may not match msa
                       if columns were removed from msa as a result of post-processing.)
                    2. Load results from Matlab plmDCA script into a Dataframe, rank by top score,
                       and remove i,j pairs with sequence distance < 5.
                    3. Map DCA indices to PDB-distance-matrix indices
                    4. Filter DCA pairs that lie on the interface
    Plots:
                    1. PDB distance matrix (optional - see PDB protocol)
                    2. PDB distance matrix and DCA top [n_pairs]
                    3. True positives vs. top [n_pairs]
                    4. DCA score vs. PDB distance

    :param msa_name:
    :param n_pairs:
    :param cutoff:
    :param results_dir:
    :param count:
    :return:
    """

    # PDB PROTOCOL ====================================================================================================
    atom = 'aa'  # cutoff type aa: all-atom, ca: c-alpha
    pdb_df_list, chain_lengths = pipeline_distance_matrix(msa_name, cutoff_type=atom, cutoff=cutoff,
                                                          read=True, plot=False)
    pdb_chain1_length = chain_lengths[0]
    pdb_total_length = sum(chain_lengths)

    # DCA PROTOCOL ====================================================================================================
    top_pairs = int(np.ceil(3 * pdb_chain1_length / 2.0))
    print("\n3L_a/2 Pairs: |{}|\tN_pairs: {}".format(top_pairs, n_pairs))
    uniprot_chain_lengths = get_lengths(msa_name)

    df_dca_unmapped = load_dca_to_df(msa_name, results_dir)
    df_dca_mapped = pipeline_mapping(msa_name, df_dca_unmapped, uniprot_chain_lengths)
    df_dca_mapped_inter = pipeline_interface(df_dca_mapped, msa_name, pdb_chain1_length)
    df_dca_mapped_inter = df_dca_mapped_inter[:n_pairs]
    df_dca_mapped_inter_dist = pipeline_dca_distance(msa_name, df_dca_mapped_inter, atom, read=True, n_pairs=n_pairs)

    # PLOTS PROTOCOL ==================================================================================================
    df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm if you dont want to plot dca
    plot_cm(pdb_df_list, cutoff, pdb_chain1_length, pdb_total_length, atom=atom, df_dca=df_dca_mapped_inter,
            msa_name=msa_name, other_dca=df_empty)

    return df_dca_mapped_inter_dist

