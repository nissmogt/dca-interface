import logging
import os
import glob
import time
from dca_performance import run_analysis, map_cm
from eigen_decomp import eigen_decomp
from score_compare import score_compare


def main(n_pairs, cutoff, sasa=False, batch=False, msa_name=None, calc_eff=False):
    """
    Analyzes DCA results for every msa. Saves contact map and TPR plots to img_dir.
    :param calc_eff:
    :param n_pairs: int - Number of DCA predictions
    :param cutoff: int - Contact map distance cutoff
    :param sasa:
    :param batch: boolean - if True, runs analysis for all MSAs in msa_directory
    :param msa_name: optional - str - Name of single MSA file. Used iff batch=False
    :return: List of MSAs for DCA to run.
    """

    # list of all MSAs in directory
    from score_vs_distance import score_vs_distance
    msa_list = glob.glob('{}*.fas'.format(msa_directory))

    if calc_eff:
        import PreProcess as pp
        print("Calculating Beff...")
        for msa_file in msa_list:
            print("{}".format(msa_file))
            pdbid = os.path.basename(msa_file).strip(".fas")[:4].lower()
            p = pp.Preprocess("PDB_benchmark_structures\\{}.cif".format(pdbid), msa_file, cutoff)
            head, seq = p.parse_fasta()
            k = p.mk_msa(seq)
            Beff = k["neff"]
            print("Beff: {}".format(Beff))
    # -- Batch run --
    else:
        # --MAIN LOOP--
        if batch:
            from dca_distance_plot import dca_distance_plot
            et = []  # list of times for each msa file
            start_time = time.time()
            count = 0
            for msa_file in msa_list:
                # --MAIN LOOP--
                print("\n-- Running analysis for {} --\n".format(msa_file))
                df_map, df_map_g = run_analysis(msa_file, n_pairs, cutoff, results_directory, msa_directory,
                                                pdb_directory, sasa_calc=sasa, count=count)
                score_vs_distance(msa, n, cut)
                # id = os.path.basename(msa_file).strip('.fas')
                # eigen_decomp(id)
                # score_compare(df_map, df_map_g, n_pairs=n, msa_name=id, title='plmDCA_GRM_c{}'.format(cutoff))
                # map_cm(msa_file, cutoff, results_directory, pdb_directory)

                count += 1
            # --END LOOP--
            print("-- Total time to run: {} --".format(time.time() - start_time))
            print("=" * 72)

        # -- Single MSA run --
        else:
            df_map = run_analysis(msa_name, n_pairs, cutoff, results_directory, msa_directory,
                                  pdb_directory, sasa_calc=sasa)
            score_vs_distance(msa, n, cut)
            # map_cm(msa_name, cutoff, results_directory, pdb_directory)
    return df_map


if __name__ == '__main__':
    import PreProcess as pp
    import shannon

    # relevant directories
    # msa_directory = 'PDB_benchmark_alignments\\'
    msa_directory = 'PDB_benchmark_alignments\\'
    pdb_directory = 'PDB_benchmark_structures\\'
    results_directory = 'results\\'
    # msa_list = ['1EFP_A_1EFP_B.fas', '1I1Q_A_1I1Q_B.fas', '1QOP_A_1QOP_B.fas', '1RM6_A_1RM6_C.fas', '1W85_A_1W85_B.fas',
    #             '1TYG_B_1TYG_A.fas', '2D1P_B_2D1P_C.fas', '2NU9_A_2NU9_B.fas', '2Y69_B_2Y69_C.fas', '2Y69_A_2Y69_B.fas',
    #             '3A0R_A_3A0R_B.fas', '3IP4_A_3IP4_B.fas', '3IP4_A_3IP4_C.fas', '3IP4_B_3IP4_C.fas', '3PNL_A_3PNL_B.fas',
    #             '3G5O_A_3G5O_B.fas', '3OAA_H_3OAA_G.fas', '3PNL_A_3PNL_B.fas']
    # msa_list = ["1EFP_A_1EFP_B.fas", "1QOP_A_1QOP_B.fas", "1I1Q_A_1I1Q_B.fas", "2ONK_A_2ONK_C.fas", "2NU9_A_2NU9_B.fas",
    #             "3G5O_A_3G5O_B.fas", "3A0R_A_3A0R_B.fas", "3RPF_A_3RPF_D.fas"]

    n = 300
    cut = 12  # 8, 10, 12
    id = '1EM8'
    c1 = 'D'
    c2 = 'C'
    msa = msa_directory + '{}_{}_{}_{}.fas'.format(id, c1, id, c2)
    df_map = main(n, cut, msa_name=msa)

    # Calc Shannon Entropy
    # aln, l_seq, ind = shannon.parseMSA(msa, alnformat='fasta', verbose=1)
    # la, gl = shannon.shannon_entropy_list_msa(aln)
    # shannon.plot(ind, la, verbose=1, msa_name='{}_{}_{}_{}'.format(id, c1, id, c2))

    # i=6
    # df_map, df_map_g = main(n, cut, msa_name=msa_directory+msa_list[i])
    # compare_index_plot(df_map, df_map_g, n_pairs=n, msa_name=msa_list[i].strip('.fas'), title='IDX_GRM_GRMmatlabprior'.format(n))

    from dca_distance_plot import dca_distance_plot
    # for msa in msa_list:
    #     score_vs_distance(msa, n, cut)
    # d = dca_distance_plot(msa, n, cut)
    #     id = msa.strip('.fas')
    #     df_map, df_map_g = main(n, cut, msa_name=msa_directory+msa)
    #     compare_index_plot(df_map, df_map_g, n_pairs=n, msa_name=id, title='IDX_GRM_GRMmatlabprior')

    # Run batch
    # dca, grm = main(n, cut, batch=True)
