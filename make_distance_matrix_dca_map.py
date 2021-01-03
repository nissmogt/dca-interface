import os
from dca_functions import cm_make
import numpy as np
import pandas as pd
from distance_pdb import plot_cm
from load_dca import load_dca_to_df
from process_coupling_matrix import process_coupling_matrix_output_scores, rank_hamming
from output_pipeline import (
    pipeline_pdb_distance_matrix, pipeline_mapping, pipeline_interface,
    pipeline_dca_distance
)


def run_analysis(msa_name, n_pairs, cutoff, results_dir, matrix, fni=False, plot=False):
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
                    4. Filter DCA pairs that lie on the pre-defined interface
    Plots:
                    1. PDB distance matrix (optional - see PDB protocol)
                    2. PDB distance matrix and DCA top [n_pairs]

    :param matrix:
    :param msa_name:
    :param n_pairs:
    :param cutoff:
    :param results_dir:
    :param fni:
    :param plot:
    :return:
    """

    # PDB PROTOCOL ====================================================================================================
    atom = 'aa'  # cutoff type aa: all-atom, ca: c-alpha
    pdb_df_list, chain_lengths = pipeline_pdb_distance_matrix(msa_name, cutoff_type=atom, cutoff=cutoff,
                                                              read=True, plot=False)
    df_pdb = pdb_df_list[0]
    pdb_chain1_length = chain_lengths[0]
    pdb_total_length = chain_lengths[2]

    # DCA PROTOCOL ====================================================================================================

    if fni:
        fni_matrix_dir = "scrambled_results\\fni_matrices\\"
        npy_matrix = np.load("{}matrix_FNi_{}.npy".format(fni_matrix_dir, msa_name))
        df_dca_mapped = cm_make(npy_matrix)

    df_dca = process_coupling_matrix_output_scores(matrixFile=matrix, freqFile=False)
    df_dca_ranked = rank_hamming(df_dca)
    df_dca_mapped = pipeline_mapping(msa_name, df_dca_ranked, read=True)
    df_dca_mapped_dist = df_dca_mapped.merge(df_pdb, how='inner', on=['i', 'j'])

    # PLOTS PROTOCOL ==================================================================================================
    if plot:
        df_plot = df_dca_mapped[:n_pairs]
        df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm if you dont want to plot dca

        folder = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\"
        img_dir = "{}{}".format(folder, results_dir)
        if not os.path.exists(img_dir):
            os.makedirs(img_dir)
        imgName = '{}APC_{}_top{}'.format(img_dir, msa_name, n_pairs)

        plot_cm(pdb_df_list, cutoff, pdb_chain1_length, pdb_total_length,
                atom=atom, df_dca=df_plot, msa_name=msa_name,
                other_dca=df_empty, img_name=imgName)

    return df_dca_mapped_dist
