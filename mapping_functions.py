def align_dca2pdb(msa_name, pdbseq, msaseq):
    """
    Some code taken from
    https://github.com/bsir/dca-frustratometer/blob/master/dca_frustratometer.py
    :param msa_name:
    :param pdbseq: PDB seq string
    :param msaseq: MSA seq string
    :return:
    """
    import numpy as np
    from Bio import pairwise2
    from evcouplings.compare import mapping as evmp
    print("(map_dca2pdb)\taligning dca sequence to pdb coord sequence...")

    # need to penalize for opening and adding gaps otherwise mapping is off (s param {-.5,-.1})
    alignments_1 = pairwise2.align.globalxs(pdbseq[0], msaseq[0], -.5, -.1)
    alignments_2 = pairwise2.align.globalxs(pdbseq[1], msaseq[1], -.5, -.1)
    print(pairwise2.format_alignment(*alignments_1[0], full_sequences=True))
    print(pairwise2.format_alignment(*alignments_2[0], full_sequences=True))

    map_1 = evmp.map_indices(alignments_1[0][0], 1, 0, alignments_1[0][1], 1, 0)
    map_2 = evmp.map_indices(alignments_2[0][0], 1 + len(pdbseq[0]), 0,
                             alignments_2[0][1], 1 + len(msaseq[0]), 0)
    map_pdb_dca = map_1.append(map_2)
    map_pdb_dca = map_pdb_dca.rename(columns={"i": "pdb_i", "A_i": "pdb_res", "j": "dca_i", "A_j": "dca_res"})

    outfile = "results\\reference_maps\\ref_map_{}.txt".format(msa_name.strip(".fas"))
    np.savetxt(outfile, map_pdb_dca, header="pdb_i\tpdb_res\tdca_i\tdca_res", fmt="%s\t%s\t%s\t%s", comments='')
    print("(map_dca2pdb)\tWrote {}".format(outfile))

    map_pdb_dca = map_pdb_dca.dropna()
    map_dca2pdb_dict = dict(zip(map_pdb_dca["dca_i"], map_pdb_dca["pdb_i"]))
    return map_dca2pdb_dict


def apply_map(array, map_pdb_dca):
    """

    :param array:
    :param map_pdb_dca:
    :return:
    """
    import numpy as np
    print("(apply_map)")

    map_dca_list = []
    for i, j, fn_apc, fn in array:
        i = int(i)
        j = int(j)
        # if i in map_pdb_dca.keys() and j in map_pdb_dca.keys():
        if str(i) in map_pdb_dca.keys() and str(j) in map_pdb_dca.keys():
            map_index_i = int(map_pdb_dca[str(i)])
            map_index_j = int(map_pdb_dca[str(j)])
            # map_index_i = int(map_pdb_dca[i])
            # map_index_j = int(map_pdb_dca[j])
            if map_index_i and map_index_j >= 0:
                map_dca_list.append([map_index_i, map_index_j, fn_apc, fn, i, j])

    return np.array(map_dca_list)


def map_dict(msa_name):
    import pandas as pd
    import numpy as np
    from get_region import get_dca_indices
    sifts_table_file = "databases/sifts/pdb_chain_uniprot_plus.csv"
    s = pd.read_csv(sifts_table_file, comment="#")
    pdbid = msa_name[:4].lower()
    chain_1 = msa_name.split("_")[1]
    chain_2 = msa_name.split("_")[3]

    pdb_start_chain_1 = s.query("pdb_id == @pdbid and pdb_chain == @chain_1").coord_start.values
    pdb_start_chain_2 = s.query("pdb_id == @pdbid and pdb_chain == @chain_2").coord_start.values
    pdb_end_chain_1 = s.query("pdb_id == @pdbid and pdb_chain == @chain_1").coord_end.values
    pdb_end_chain_2 = s.query("pdb_id == @pdbid and pdb_chain == @chain_2").coord_end.values

    uniprot_start_chain_1 = s.query("pdb_id == @pdbid and pdb_chain == @chain_1").uniprot_start.values
    uniprot_end_chain_1 = s.query("pdb_id == @pdbid and pdb_chain == @chain_1").uniprot_end.values
    uniprot_start_chain_2 = s.query("pdb_id == @pdbid and pdb_chain == @chain_2").uniprot_start.values
    uniprot_end_chain_2 = s.query("pdb_id == @pdbid and pdb_chain == @chain_2").uniprot_end.values

    # pdb
    pdb_start_chain_1 = ([int(i) for i in pdb_start_chain_1])
    pdb_end_chain_1 = ([int(i) for i in pdb_end_chain_1])
    # add last index of end index + 1 for chain 2
    pdb_start_chain_2 = ([(int(i) + pdb_end_chain_1[-1] + 1) for i in pdb_start_chain_2])
    pdb_end_chain_2 = ([(int(i) + pdb_end_chain_1[-1] + 1) for i in pdb_end_chain_2])

    # uniprot
    uniprot_start_chain_1 = ([int(i) for i in uniprot_start_chain_1])
    uniprot_end_chain_1 = ([int(i) for i in uniprot_end_chain_1])
    # add last index of end index + 1 for chain 2
    uniprot_start_chain_2 = ([(int(i) + uniprot_end_chain_1[-1] + 1) for i in uniprot_start_chain_2])
    uniprot_end_chain_2 = ([(int(i) + uniprot_end_chain_1[-1] + 1) for i in uniprot_end_chain_2])

    pdb_start_indices = pdb_start_chain_1 + pdb_start_chain_2
    pdb_end_indices = pdb_end_chain_1 + pdb_end_chain_2
    uniprot_start_indices = uniprot_start_chain_1 + uniprot_start_chain_2
    uniprot_end_indices = uniprot_end_chain_1 + uniprot_end_chain_2

    pdb_indices = make_indices(pdb_start_indices, pdb_end_indices)
    # pdb_indices = range(1, len(pdb_indices))
    uniprot_indices = make_indices(uniprot_start_indices, uniprot_end_indices)

    dca_indices = get_dca_indices(msa_name)
    uni2pdb = dict(zip(uniprot_indices, pdb_indices))
    dca2uni = dict(zip(dca_indices, uniprot_indices))
    dca2pdb = dict(zip(dca_indices, pdb_indices))
    pdb2uni = dict(zip(pdb_indices, uniprot_indices))
    # print(dca2pdb)
    return uni2pdb, dca2uni, dca2pdb, pdb2uni


def make_indices(indices_1, indices_2):
    import numpy as np
    x = []
    for i in range(len(indices_1)):
        x.append(np.array(range(indices_1[i], indices_2[i])))
    return np.concatenate(x)


