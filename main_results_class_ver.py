import PreProcess as pp
import os
import glob
n = 30
cut = [10, 15, 20, 25, 30]

msa_directory = 'PDB_benchmark_alignments\\'
pdb_directory = 'PDB_benchmark_structures\\'
results_directory = 'results\\'
msa_list = glob.glob("{}*.fas".format(msa_directory))

msa = "PDB_benchmark_alignments\\3G5O_A_3G5O_B.fas"
msa_name = os.path.basename(msa)
pdb = "PDB_benchmark_structures\\{}.cif".format(msa_name[:4].lower())
# score_matrix = "results\\FNi_{}.txt".format(msa_name.strip(".fas"))
print("(pdb filename)\t{}".format(os.path.basename(pdb)))
p = pp.Preprocess(pdb, msa, cut[2])
seqa, seqb = p.get_residues(seq=True)
s = seqa + seqb
ma, mb = p.msa_template(split=True, len_a=p.chain_lengths[0])
m = ma+mb
# dina, pina, mapa, mappa = p.get_msa_to_pdb_dictionary(seqa, ma)
# dinb, pinb, mapb, mappb = p.get_msa_to_pdb_dictionary(seqb, mb)
# dinb[0] + len(dina)
di, pi, md, mp = p.get_msa_to_pdb_dictionary(s, m)
# df_mon, df_inter = p.read_distance_matrix_file()
# print(p.read_length_file())
# for c in cut:
#     for msa_file in msa_list:
#         msa_name = os.path.basename(msa_file)
#         pdbid = msa_name[:4]
#         pdbfile = "{}{}.cif".format(pdb_directory, pdbid.lower())
#         print("(pdb filename)\t{}".format(os.path.basename(pdbfile)))
#         p = pp.Preprocess(pdbfile, msa_file, c)
#         df_pdb, df_mon, df_inter = p.distance_matrix()

# df_gremlin = p.read_dca_matrix(score_matrix, gremlin_file=True)
# pdbseq = p.combine_chain_sequence()
# msaseq = p.msa_template()
# dca_index, pdb_index, map_to_dca, map_to_pdb = p.get_msa_to_pdb_dictionary(pdbseq, msaseq)
# df_map_dca = p.backmap(map_to_pdb, dca_index[0], score_matrix)
