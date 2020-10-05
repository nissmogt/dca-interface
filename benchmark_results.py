import logging
import pandas as pd
import os
import glob
import time
from make_distance_matrix_dca_map import run_analysis
from compare_new_model_with_old import compare_new_scrambled_with_paired
from scramble_sequence import scramble_sequence


def main(n_pairs, cutoff, batch=False, msa_name=None):
    """
    Analyzes DCA results for every msa. Saves contact map and TPR plots to img_dir.
    :param n_pairs: int - Number of DCA predictions
    :param cutoff: int - Contact map distance cutoff
    :param batch: boolean - if True, runs analysis for all MSAs in msa_directory
    :param msa_name: optional - str - Name of single MSA file. Used iff batch=False
    :return: List of MSAs for DCA to run.
    """

    # list of all MSAs in directory
    msa_directory = 'PDB_benchmark_alignments\\'
    msa_directory = 'PDB_benchmark_alignments\\'
    pdb_directory = 'PDB_benchmark_structures\\'
    results_directory = 'results\\'
    msa_list = glob.glob('{}*.fas'.format(msa_directory))
    length_file = "uniprot_lengths_pdbid_chains.csv"

    # -- Batch run --
    if batch:
        count = 0
        systemNames = []
        start_time = time.time()
        for msa_file in msa_list[:100]:
            msa_name = os.path.basename(msa_file).strip(".fas")
            if os.path.exists("scrambled_results\\matrix_ising_{}_rep0_scrambled.fas.mat".format(
                    os.path.basename(msa_name))):
            # if os.path.exists("results\\FN_{}.txt".format(os.path.basename(msa_file))):
                if os.path.exists("PDB_benchmark_alignments\\a2m\\{}.a2m".format(msa_name)):
                    count += 1
                    print("CURRENTLY ON MSA: {} ({}/{})".format(msa_name, count, len(msa_list[:100])))
                    # dca_df = run_analysis(msa_name, n_pairs, cutoff, results_directory)
                    compare_new_scrambled_with_paired(msa_name, n_pairs, cutoff)
                    # k = scramble_sequence(msa_name, n_replicates=5)
                    # systemNames.append("{}.fas\n".format(msa_name))
        # --END LOOP--
        print("-- Total time to run: {} --".format(time.time() - start_time))
        print("=" * 72)
        # f = open("benchmark_systems_count{}.txt".format(count), 'w')
        # f.writelines(systemNames)
        # f.close()
        return systemNames

    # Single MSA run
    else:
        msa_name = os.path.basename(msa_name).strip(".fas")
        # dca_df = run_analysis(msa_name, n_pairs, cutoff, results_directory)
        dca_df = compare_new_scrambled_with_paired(msa_name, n_pairs, cutoff)
        return dca_df


if __name__ == '__main__':
    n = 50
    cut = 12
    pdbid = '1GL2'
    c1 = 'A'
    c2 = 'D'
    msa = '{}_{}_{}_{}.fas'.format(pdbid, c1, pdbid, c2)
    # Single run
    df_dca = main(n, cut, msa_name=msa)
    # Run batch
    # x = main(n, cut, batch=True)
