from plot_cm import plot_cm, draw_dca_tcl
from pdb import pdb_map, make_list, get_lengths_seq, map_msa_to_pdb, cat_seq, cm_make, apply_map
from dca_performance import plot_performance
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import glob
import csv
import time
import logging

logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.basicConfig(filename='log_benchmark_results.log', level=logging.DEBUG)
# relevant directories
msa_directory = 'PDB_benchmark_alignments\\'
pdb_directory = 'PDB_benchmark_structures\\'
msa_name = '3IP4_A_3IP4_C'
logging.debug("Loading {}...".format(msa_name))

# dca files
# gfile = 'gremlin\\GRM_2Y69_B_2Y69_C.txt'
try:
    vanilla_dca = 'results\\fn_' + msa_name + '_plmdca_rt2.txt'
    dca_score_matrix = 'results\\FNi_' + msa_name + '.txt'
    contact_map = 'results\\cm_' + os.path.basename(dca_score_matrix)

    # Cat PDB chain seqs together and find msa to pdb alignment
    print("\tcatting seqs...")
    full_seq, msa_seq, chains = cat_seq(msa_name, msa_directory, pdb_directory, get_length=None)
    print("\tcreating map...")
    dca_indices, pdb_indices, map_to_dca, map_to_pdb = map_msa_to_pdb(full_seq, msa_seq)
    total_length = len(pdb_indices)

    # parameters
    n_pairs = 50
    step = 10
    cutoff = 15
    pdbfile = pdb_directory + msa_name[:4] + ".cif"
    df_pdb, df_mon, df_inter = pdb_map(pdbfile, chains, cutoff)

    # make contact map and apply backmapping
    print("\tapplying map...")
    df_umap = cm_make(dca_score_matrix)
    df_map = cm_make(dca_score_matrix, map_to_pdb, dca_indices[0])


    # Get length of first chain
    print("\tgetting lengths...")
    with open('lengths.csv', 'r', newline='\n') as csvfile:
        r = csv.reader(csvfile)
        for row in r:
            if msa_name == row[0]:
                length_a = int(row[2]) + 1

    # Plot for a range of number of DCA pairs
    for i in range(10, n_pairs + step, step):
        logging.info("Plotting top {}...".format(i))
        contact_map_file = "results\\cm_FNi_" + msa_name + ".txt"
        draw_dca_tcl(contact_map_file, n_pairs, length_a, chains)
        plot_cm(pdbfile, df_umap, df_map, chains, i, cutoff, length_a, total_length,
            gremlinfile=None, vanillafile=None, title="map_vs_umap")

    # Plots TPR for top DCA predictions (compared to PDB interface pairs)
    plot_performance(df_inter, df_map, n_pairs, total_length, msa_name, cutoff)

except OSError:
    logging.debug("{}.fas has not been analysed.".format(msa_name))
    print("Not yet analyzed.")


def multiple_msa_cm(msa_directory):
    msa_list = glob.glob(msa_directory + '*.fas')
    start_time = time.time()
    for msa_file in msa_list:
        msa_name = os.path.basename(msa_file).strip(".fas")
        logging.debug("Loading {}...".format(msa_name))

        # dca files
        # gfile = 'gremlin\\GRM_2Y69_B_2Y69_C.txt'
        try:
            logging.captureWarnings(True)
            vanilla_dca = 'results\\fn_' + msa_name + '_plmdca_rt2.txt'
            dca_score_matrix = 'results\\FNi_' + msa_name + '.txt'
            contact_map = 'results\\cm_' + os.path.basename(dca_score_matrix)

            # Cat PDB chain seqs together and find msa to pdb alignment
            full_seq, msa_seq, chains = cat_seq(msa_name, msa_directory, pdb_directory)
            dca_indices, pdb_indices, map_to_dca, map_to_pdb = map_msa_to_pdb(full_seq, msa_seq)
            length = len(pdb_indices)

            # make contact map and apply backmapping
            df_umap = cm_make(dca_score_matrix)
            df_map = cm_make(dca_score_matrix, map_to_pdb, dca_indices[0])

            # parameters
            n_pairs = 40
            step = 10
            cutoff = 12
            pdbfile = pdb_directory + msa_name[:4] + ".cif"

            # Get length of first chain
            with open('lengths.csv', 'r', newline='\n') as csvfile:
                r = csv.reader(csvfile)
                for row in r:
                    if msa_name == row[0]:
                        length_a = int(row[2]) + 1

            # Plot for a range of number of DCA pairs
            for i in range(20, n_pairs + step, step):
                logging.info("Plotting top {}...".format(i))
                plot_cm(pdbfile, df_umap, df_map, chains, i, cutoff, length_a, length,
                        gremlinfile=None, vanillafile=None, title="map_vs_umap")
        except OSError:
            logging.debug("{}.fas has not been analysed.".format(msa_name))
            continue

    logging.info("Total time: {}".format(time.time() - start_time))
    # draw_dca_tcl(contact_map, n_pairs, length, chain)
