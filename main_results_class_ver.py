import PreProcess as pp
n = 30
c = 16

msa = '1B70_A_1B70_B.fas'
pdb = '1b70.cif'
score_matrix = "results\\FNi_1B70_A_1B70_B.txt"
p = pp.Preprocess(pdb, msa, c)
pdbseq = p.combine_pdb_seq()
msaseq = p.read_msa()
dca_index, pdb_index, map_to_dca, map_to_pdb = p.map_msa_to_pdb(pdbseq, msaseq)
df_map_dca = p.backmap(map_to_pdb, dca_index[0], score_matrix)
