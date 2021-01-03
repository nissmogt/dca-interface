import pandas as pd
import numpy as np
import os
from distance_dca import distance_dca, read_dca_distance_matrix
from dca_functions import cm_make
from distance_pdb import plot_cm
from read_db import get_lengths
from output_pipeline import (
    pipeline_pdb_distance_matrix, pipeline_mapping,
    pipeline_dca_distance
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
    pdb_df_list, chain_lengths = pipeline_pdb_distance_matrix(msa_name, cutoff_type=atom, cutoff=cutoff,
                                                              read=True, plot=False)
    pdb_chain1_length = chain_lengths[0]
    pdb_total_length = sum(chain_lengths)

    pdb_file = "distance_matrix\\heavy_atom_distance_matrix_{}.txt".format(msa_name)
    df_pdb = pd.read_csv(pdb_file, delimiter='\t')

    # DCA PROTOCOL ====================================================================================================
    uniprot_chain_lengths = get_lengths(msa_name)
    outDir = "scrambled_results\\fni_matrices\\"
    npy_matrix = np.load("{}matrix_FNi_{}.npy".format(outDir, msa_name))
    df_fni_apc = cm_make(npy_matrix)
    df_fni_apc_mapped = pipeline_mapping(msa_name, df_fni_apc, uniprot_chain_lengths, read=True)
    # np.savetxt('{}{}_FNi_mapped.txt'.format(outDir, msa_name), df_fni_apc_mapped, delimiter='\t')
    df_fni_apc_mapped = df_fni_apc_mapped

    # and this name - dont forget
    fni_dist_file = "{}{}_FNi_mapped_top{}.txt".format(outDir, msa_name, n_pairs)
    if os.path.exists(fni_dist_file):
        df_fni_apc_mapped_dist = pd.read_csv(fni_dist_file, delimiter="\t")
        df_fni_apc_mapped_dist[:n_pairs]
    else:
        df_fni_apc_mapped_dist = distance_dca(df_fni_apc_mapped, msa_name, atom=atom, other_name=True)
        np.savetxt(fni_dist_file, df_fni_apc_mapped_dist[:n_pairs],
                   header="i\tj\tscore\tdist_aa\tchain_1\tchain_2\tresnames\tatom_id",
                   fmt='%d\t%d\t%f\t%f\t%s\t%s\t%s\t%s', comments='')

    df_dca_mapped_inter_dist = pipeline_dca_distance(msa_name, None, atom, read=True, n_pairs=50)
    df_dca_mapped_inter_dist = df_dca_mapped_inter_dist

    # PLOTS PROTOCOL ==================================================================================================
    # df_empty = pd.DataFrame({'A': []})  # an empty Dataframe to use in plot_cm if you dont want to plot dca

    # DONT FORGET TO CHANGE IMG NAME
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\SEPT_2020\\plmDCA_scrambled\\"
    imgName = img_dir + os.path.basename(fni_dist_file).strip(".txt")
    plot_cm(pdb_df_list, cutoff, pdb_chain1_length, pdb_total_length, atom=atom, df_dca=df_dca_mapped_inter_dist[:10],
            msa_name=msa_name, other_dca=df_fni_apc_mapped_dist, img_name=imgName)
    return df_fni_apc_mapped_dist
