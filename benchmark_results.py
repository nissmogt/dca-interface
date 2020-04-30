import logging
import os
import glob
import time
from dca_performance import run_analysis


# msa_name = '3IP4_A_3IP4_C'


def main_program(n_pairs, cutoff, run_again=False):
    """
    Analyzes DCA results for every msa. Saves contact map and TPR plots to img_dir.
    :param n_pairs: int - Number of DCA predictions
    :param cutoff: int - Contact map distance cutoff
    :param run_again: Optional - Boolean (Input True to ignore image check and rerun analysis)
    :return: List of MSAs for DCA to run.
    """
    # relevant directories
    msa_directory = 'PDB_benchmark_alignments\\'
    pdb_directory = 'PDB_benchmark_structures\\'
    results_directory = 'results\\'
    img_directory = "C:\\Users\\kmehr\\Google Drive\\work\\images\\2020\\APR_2020\\Week_4\\"

    # list of all MSAs in directory
    msa_list = glob.glob('{}*.fas'.format(msa_directory))
    # the name is based on how plmDCA outputs are named
    raw_img_list = '\t'.join(glob.glob('{}FNi_*'.format(img_directory)))
    raw_results_list = '\t'.join(glob.glob("{}FNi_*".format(results_directory)))

    # Log stuff
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.basicConfig(filename='log_benchmark_analysis.log', level=logging.DEBUG)

    to_analyze_list = []  # initialize list of MSAs to-be-analyzed
    analyzed_list = []    # initialize list for analyzed MSAs

    start_time = time.time()
    # --MAIN LOOP--
    for msa_file in msa_list:
        msa_name = os.path.basename(msa_file).strip(".fas")
        print("\tLoading {}...".format(msa_name))

        # check if DCA results for current MSA exists, if exists go to next check
        if msa_name in raw_results_list:
            print("\n(result check)\t{}{} EXISTS...".format(results_directory, msa_name))

            # check if an image file exists, if not exists continue with analysis
            if msa_name not in raw_img_list and not run_again:
                print("(new_run)\tRunning analysis for {}...".format(msa_name))
                # -- ANALYSIS FUNCTION CALL --
                run_analysis(msa_name, n_pairs, cutoff, results_directory, msa_directory,
                             pdb_directory, img_directory)

            elif msa_name in raw_img_list and run_again:
                print("(rerun)\tRun again - analysis for {}...".format(msa_name))
                # -- ANALYSIS FUNCTION CALL --
                run_analysis(msa_name, n_pairs, cutoff, results_directory, msa_directory,
                             pdb_directory, img_directory)

            else:
                # make a list of analyzed MSAs
                analyzed_list.append(msa_name)
        else:
            # make a list of MSAs to-be-analyzed
            to_analyze_list.append(msa_name)

    # --END LOOP--
    print("-- Total time: {} --".format(time.time() - start_time))
    print("=" * 72)
    if len(analyzed_list) > 0:
        print("Already analyzed and in\n{}\n{}\n".format(img_directory, analyzed_list))
    # print out list of un-analyzed MSAs if list exists
    if len(to_analyze_list) > 0:
        print("NOTE!\tMSAs NOT analyzed yet:\n{}".format('\n'.join(to_analyze_list)))

    return to_analyze_list


n = 50
c = 16
analysis_list = main_program(n, c, run_again=True)
