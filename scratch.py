def map_to_uniprot(msa_name, msaseq, dca_indices):
    import pandas as pd
    # msa_name = "2OXG_Z_2OXG_Y"
    sifts_table_file = "databases/sifts/pdb_chain_uniprot_plus.csv"
    s = pd.read_csv(sifts_table_file, comment="#")
    pdbid = msa_name[:4].lower()
    chain_1 = msa_name.split("_")[1]
    chain_2 = msa_name.split("_")[3]
    pdbstart_1 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_1").coord_start.values[0])
    pdbend_1 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_1").coord_end.values[0])
    pdbstart_2 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_2").coord_start.values[0])
    pdbend_2 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_2").coord_end.values[0])

    p1 = range(pdbstart_1, pdbend_1)
    p2 = range(pdbstart_2, pdbend_2)
    p_length = len(p1) + len(p2)

    uniprot_start_1 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_1").uniprot_start.values[0])
    uniprot_end_1 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_1").uniprot_end.values[0])
    uniprot_start_2 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_2").uniprot_start.values[0])
    uniprot_end_2 = int(s.query("pdb_id == @pdbid and pdb_chain == @chain_2").uniprot_end.values[0])

    m1 = range(uniprot_start_1, uniprot_end_1)
    m2 = range(uniprot_start_2, uniprot_end_2)

    m_length = len(m1) + len(m2)

    # n = min(m_length, p_length)
    dca_start = dca_indices[0]

    n = min(len(msaseq), p_length)
    pdb_index = range(pdbstart_1, n + pdbstart_1)
    uniprot_index = range(uniprot_start_1, n + uniprot_start_1)
    dca_index = range(dca_start, n + dca_start)

    return dict(zip(dca_index, pdb_index))
    # d = dict(zip(pdb_index, uniprot_index))
