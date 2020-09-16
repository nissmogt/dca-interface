def get_msa_to_pdb_dictionary(msa_name, pdbseq, msaseq):
    """
    Taken from
    https://github.com/bsir/dca-frustratometer/blob/master/dca_frustratometer.py
    :param msa_name:
    :param pdbseq: PDB seq string
    :param msaseq: MSA seq string
    :return:
    """
    import numpy as np
    import re
    from Bio import pairwise2
    fname = "(map_msa_to_pdb)"
    print("Using pairwise2...")
    alignments = pairwise2.align.localxx(msaseq, pdbseq)

    fastastart = re.search("[A-Z]", alignments[0][1]).start()
    pdbstart = re.search("[A-Z]", alignments[0][0]).start()
    print("fastastart: {}\tpdbstart: {}".format(fastastart, pdbstart))
    print(pairwise2.format_alignment(*alignments[0], full_sequences=True))

    n = min(len(msaseq), len(pdbseq))
    pdb_indices = range(pdbstart, n + pdbstart)
    dca_indices = range(fastastart, n + fastastart)
    outfile = "ref_map_{}.txt".format(msa_name.strip(".fas"))
    np.savetxt(outfile, np.transpose([dca_indices, pdb_indices]), fmt="%s", header="MSA  PDB")
    print("{}\tWrote {}".format(fname, outfile))
    map_to_dca = dict(zip(pdb_indices, dca_indices))
    map_to_pdb = dict(zip(dca_indices, pdb_indices))

    print(fname + "\t\tn: %s pdb indices: %s dca indices: %s" % (n, pdb_indices, dca_indices))
    print(fname + "\tFinished aligning msa sequences to pdb...")
    return alignments, dca_indices, pdb_indices, map_to_dca, map_to_pdb
    # return dca_indices, pdb_indices, map_to_dca, map_to_pdb
