import os
import pandas as pd
from make_distance_matrix_dca_map import run_analysis
from calculate_sasa import batch_calculate
from draw_contacts import batch_draw
from compare_new_model_with_old import compare_new_scrambled_with_paired


def main(n_pairs, cutoff, resDir, matFile, msa_name=None):
    """
    Analyzes DCA results for every msa. Saves contact map and TPR plots to img_dir.
    :param matFile:
    :param n_pairs: int - Number of DCA predictions
    :param cutoff: int - Contact map distance cutoff
    :param resDir:
    :param msa_name: optional - str - Name of single MSA file. Used iff batch=False
    """
    dca_df = run_analysis(msa_name, n_pairs, cutoff, resDir, matFile, fni=False, plot=False)
    return dca_df


results_directory = ["nonbonded_restraints_results\\20A\\", "nonbonded_restraints_results\\15A\\",
                     "nonbonded_restraints_results\\12A\\", "nonbonded_restraints_results\\8A\\",
                     "sasa_restraints_results\\sasa_5\\", "vanilla_results\\"]

matrix_directory = ["coupling_matrices\\nonbonded\\20A\\", "coupling_matrices\\nonbonded\\15A\\",
                    "coupling_matrices\\nonbonded\\12A\\", "coupling_matrices\\nonbonded\\8A\\",
                    "coupling_matrices\\sasa_restraints\\", "coupling_matrices\\vanilla\\"]

dimers = ["1EM8_D_1EM8_C", "1FM0_E_1FM0_D", "1KA9_H_1KA9_F", "1ZT2_A_1ZT2_B", "2NQ2_C_2NQ2_A", "2OXG_Z_2OXG_Y",
          "4NQW_A_4NQW_B", '5WY5_B_5WY5_A', '5L8H_B_5L8H_A', '5UNI_B_5UNI_A', '5F5S_A_5F5S_B',
          '5MU7_B_5MU7_A', '5M72_A_5M72_B', 'HKRR_C_HKRR_A']
# for rd in results_directory:
#     print(results_directory)
rd = results_directory[-1]
m = matrix_directory[-1]
if not os.path.exists(rd):
    os.makedirs(rd)
n = 100
cut = 12
# batch_calculate(calc_type='dca', result_dir=rd, inputList=dimers)
# for i in range(11):
i = -3
msa = '{}'.format(dimers[i])
matrix_file = "{}matrix_ising_{}.fas.mat".format(m, msa)
df_dca = main(n, cut, resDir=rd, matFile=matrix_file, msa_name=msa)
# output
fn_dist_file_out = "{}FN_all_{}_mapped_aa_dist.txt".format(rd, msa.strip('.fas'))
header = "i\tj\tfn_apc\tfn\tdist_aa\tui\tuj\tsi\tsj\tchain_1\tchain_2\tresnames\tatom_id"
df_dca.to_csv(fn_dist_file_out, sep='\t', index=False, header=header, float_format='%.5f')

# batch_draw(dimers, resDir=rd, nPairs=n)
