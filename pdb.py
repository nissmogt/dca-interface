"""
Functions for calculating a PDB distance matrix.
    get_residues(pdbfilename, ['Chain A', 'Chain B'])
        Returns: list of all residues in PDB chains
    calc_distances(residue a, residue b)
        Returns: distance between residue a and residue b
    pdb_map(pdbfilename, ['Chain A', 'Chain B'], distance cutoff)
        Returns:
"""
import logging


def get_pdb(pdb_list):
    from Bio.PDB import PDBList

    out_dir = "PDB_benchmark_structures\\"
    pdb = pdb_list
    number_ids = len(pdb)

    print("Downloading in %s:\n" % out_dir)
    for ids in pdb:
        print('%s' % ids)
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(ids, file_format='pdb', pdir=out_dir)
    return pdbl


def get_lengths_seq(pdbid_list, id_chain_dict):
    """
    Get lengths of individual chains in mmcif file.
    :param pdbid_list: List of pdb ids
    :param id_chain_dict: Dictionary of msa name and relevant chains
    :return: Text file of pdb id, first chain, and first chain length
    """
    from Bio.PDB import MMCIFParser, MMCIF2Dict
    import time
    import logging
    fname = "(get_lengths_seq)"
    pdb_dir = 'PDB_benchmark_structures\\'
    out = open('lengths.csv', 'w')
    out.write("PDBid,Chain,Length\n")
    start_time = time.time()
    logging.info('PDBid\tChain\tLength\n')
    # loop through each pdb in list, extract seq and seq length for relevant chain
    for ids in pdbid_list:
        pdb_id = ids[:4]
        seq_out = open(pdb_dir + str(pdb_id) + '.fasta', 'w')  # seq file name
        m = MMCIF2Dict.MMCIF2Dict(pdb_dir + pdb_id + '.cif')  # create a dict of mmcif codes
        chains = m['_entity_poly.pdbx_strand_id']  # string list of chain ids
        if '_entity_poly.pdbx_seq_one_letter_code_can' in m.keys():
            full_sequence = m['_entity_poly.pdbx_seq_one_letter_code_can']  # raw sequence for each chain
            first_chain = id_chain_dict[ids][0]  # gets first chain id from msa file
            # loop thru each chain and write seq in fasta format
            for c in chains:
                # write every seq for every chain
                for c_id in c.split(','):
                    seq_out.write('> Chain_%s\n' % c_id)
                    seq = full_sequence[chains.index(c)].replace('\n', '')
                    seq_out.write('%s\n' % seq)
                # extract only first chain length P3: maybe extract both?
                if first_chain in c.split(','):
                    length = len(seq)
                    out.write("%s,%s,%s\n" % (ids, first_chain, length))
                    logging.debug('%s\t%s\t%s' % (ids, first_chain, length))
    out.close()
    seq_out.close()
    print(fname + "\t\tloop time: ", time.time() - start_time)


def make_list(msa_dir):
    """
    Makes a list of pdbids and a dictionary of ids and chains

    :param msa_dir: Directory of filename: ID1_chain1_ID2_chain2.fas
    :return: List of pdbids and a dictionary file storing ids and chains
    """
    fname = "(make_list)"
    import glob
    import os
    import csv
    import logging
    logging.info(fname + "\twriting dictionary of pdb id and chain id...")
    # list of msa files in msa dir
    list_fasta_files = glob.glob(msa_dir + '*.fas')
    pdb_id_list = []
    chain_list = []
    # loop through each msa file and append to a master list of pdb id and relevant chains
    for ids in list_fasta_files:
        fasta_names = os.path.basename(ids).split('.fas')[0]
        pdb_id_list.append(fasta_names)
        chain_list.append([fasta_names.split('_')[1], fasta_names.split('_')[3]])
    # write pdb id and relevant chains to csv file
    with open("id_chains.csv", "w", newline='\n') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(zip(pdb_id_list, chain_list))
    # make pdb id and chain dictionary
    id_chain_dict = dict(zip(pdb_id_list, chain_list))
    logging.info(fname + "\tFinished writing dictionary of pdb id and chain id.")
    return pdb_id_list, id_chain_dict


def cat_seq(msa_name, msa_directory, pdb_directory, get_length=None):
    """
    Concatenates PDB sequences
    :param msa_name:
    :param msa_directory:
    :param pdb_directory:
    :param get_length:
    :return:
    """
    fname = "(cat_seq)"
    logging.info(fname + "\tcat-ting sequences...")
    # Get lengths of first chain for all msa files
    list_ids, pdbid_chain_dict = make_list(msa_directory)  # outputs pdbid and chains of each msa
    # writes file of pdbid and length of first chain and pdb fasta file
    if get_length:
        get_lengths_seq(list_ids, pdbid_chain_dict)
    # length_file = 'lengths.csv'

    # Create MSA-seq-template file object
    msa_file_obj = open(msa_directory + msa_name + '.fas', 'r')
    msa_header = msa_file_obj.readline().rstrip()
    msa_seq = msa_file_obj.readline().rstrip()
    msa_file_obj.close()

    # Create PDB fasta file object
    logging.info(fname + "\t\tFasta file\tNumber of chains")
    fasta_file = pdb_directory + msa_name[:4] + '.fasta'
    fasta_file_obj = open(fasta_file, 'r')
    num_of_chains = int(fasta_file_obj.read().count('\n') / 2)
    fasta_file_obj.seek(0)  # return to top of file
    # concatenate relevant chains
    fasta_header = []
    fasta_seq = []
    # Loop thru each chain and append header and seq to a list, then make into lookup table
    for line in range(num_of_chains):
        fasta_header.append(fasta_file_obj.readline().rstrip()[-1])
        fasta_seq.append(fasta_file_obj.readline().rstrip())
    logging.debug(fname + "\t\t%s\t%s" % (fasta_file, num_of_chains))
    fasta_file_obj.close()
    msa_file_obj.close()

    # Error check
    if len(fasta_seq) != len(fasta_header):
        logging.error("(cat_seq_error) Error! Number of seqs not equal to number of chains\n")

    header_seq_dict = dict(zip(fasta_header, fasta_seq))

    # For each chain in msa find relevant pdb seq and cat the two seqs
    assert len(pdbid_chain_dict[msa_name]) == 2  # length is always 2 unless code changes
    first_chain_id = str(pdbid_chain_dict[msa_name][0])
    second_chain_id = str(pdbid_chain_dict[msa_name][1])
    logging.debug(fname + "\t\tChain 1: %s  Chain 2: %s" % (first_chain_id, second_chain_id))
    first_seq = header_seq_dict[first_chain_id]
    second_seq = header_seq_dict[second_chain_id]
    full_seq = first_seq + second_seq
    logging.debug(fname + "\t\tPDB seq length: %s  MSA seq length: %s" % (len(full_seq), len(msa_seq)))
    logging.info(fname + "\tFinished cat-ting sequences.")
    return full_seq, msa_seq, [first_chain_id, second_chain_id]


def map_msa_to_pdb(pdbseq, msaseq):
    """
    Taken from
    https://github.com/bsir/dca-frustratometer/blob/master/dca_frustratometer.py
    :param pdbseq: PDB seq string
    :param msaseq: MSA seq string
    :return:
    """
    fname = "(map_msa_to_pdb)"
    logging.debug(fname + "\taligning msa sequences to pdb...")
    if msaseq in pdbseq:
        pdbstart = pdbseq.find(msaseq)
        fastastart = 0
    elif pdbseq in msaseq:
        fastastart = msaseq.find(pdbseq)
        pdbstart = 0

    else:
        import re
        from Bio import pairwise2
        alignments = pairwise2.align.globalxx(msaseq, pdbseq)

        fastastart = re.search("[A-Z]", alignments[0][1]).start()
        pdbstart = re.search("[A-Z]", alignments[0][0]).start()
        logging.debug("fastastart: {}\tpdbstart: {}".format(fastastart, pdbstart))

    n = min(len(msaseq), len(pdbseq))

    pdb_indices = range(pdbstart, n + pdbstart)
    dca_indices = range(fastastart, n + fastastart)
    map_to_dca = dict(zip(pdb_indices, dca_indices))
    map_to_pdb = dict(zip(dca_indices, pdb_indices))

    logging.debug(fname + "\t\tn: %s pdb indices: %s dca indices: %s" % (n, pdb_indices, dca_indices))
    logging.info(fname + "\tFinished aligning msa sequences to pdb...")
    return dca_indices, pdb_indices, map_to_dca, map_to_pdb


def cm_make(score_matrix, map_dictionary=None, dca_start=None):
    """
    Makes a contact map from a DCA score matrix file. Also maps dca indices to pdb.
    :param score_matrix: Frobenius norm matrix
    :param map_dictionary: Optional: Dictionary of mapping
    :return: Three-column dataframe composed of pair i, j, and fn score
    """
    fname = "(cm_make)"
    import numpy as np
    import pandas as pd
    import os
    dca = np.loadtxt(score_matrix)
    filename = score_matrix.strip(".txt")
    basename = os.path.basename(filename)
    logging.info(fname + "\t\tScore matrix filename: %s" % filename)
    x_output = []
    n = dca.shape[0]
    logging.debug(fname + "\tshape of matrix: %s" % n)
    for i in range(n - 1):
        for j in range(i + 1, n):
            x_output.append([i + 1, j + 1, dca[i, j]])
    dca_array = np.array(x_output)
    if map_dictionary and dca_start:
        logging.debug(fname + "\t\tMap dictionary given - PROCEED with mapped contact map.")
        map_dca_array = apply_map(dca_array, map_dictionary, dca_start)
        df_map_dca = pd.DataFrame(map_dca_array, columns=['i', 'j', 'score'])
        df_map_dca = df_map_dca.sort_values(ascending=False, by=['score'])
        logging.debug("sorted df_map_dca head {}".format(df_map_dca.head()))
        np.savetxt('results\\mapped_cm_' + basename + '.txt', df_map_dca, fmt='%d\t%d\t%f')
        return df_map_dca
    else:
        logging.debug(fname + "\t\tNo map dictionary given - PROCEED with unmapped contact map.")
        df_dca = pd.DataFrame(dca_array, columns=['i', 'j', 'score'])
        df_dca = df_dca.sort_values(ascending=False, by=['score'])
        logging.debug("sorted df_dca head {}".format(df_dca.head()))
        np.savetxt('results\\cm_' + basename + '.txt', df_dca, fmt='%d\t%d\t%f')
        return df_dca


def apply_map(dca_array, map_dictionary, dca_start):
    import numpy as np
    map_dca_list = []
    for i, j, score in dca_array:
        if i-1 >= dca_start and j-1 >= dca_start:
            map_index_i = map_dictionary[int(i) - 1]
            map_index_j = map_dictionary[int(j) - 1]
            map_dca_list.append([map_index_i, map_index_j, score])
    map_dca_array = np.array(map_dca_list)
    return map_dca_array


# dev branch test
def get_residues(pdbfile, chain_ids=None):
    import Bio.PDB
    import os
    import sys
    """Build a simple list of residues from a single chain of a PDB file.
    Args:
        pdb_fn: The path to a PDB file.
        chain_ids: A list of single-character chain identifiers.
    Returns:
        A list of Bio.PDB.Residue objects.
    """
    fname = "(get_residues)"
    pdb_id = os.path.splitext(os.path.basename(pdbfile))[0]
    if os.path.basename(pdbfile).split('.')[-1] == 'pdb':
        print(fname + '\tprocessing .pdb file...')
        parser = Bio.PDB.PDBParser()
        struct = parser.get_structure(pdb_id, pdbfile)
        model = struct[0]
    else:
        print(fname + '\tprocessing .cif file...')
        parser = Bio.PDB.MMCIFParser()
        struct = parser.get_structure(pdb_id, pdbfile)
        model = struct[0]

    if chain_ids is None:
        # get residues from every chain.
        chains = model.get_list()
    else:
        chains = [model[ch_id] for ch_id in chain_ids]

    residues = []

    for ch in chains:
        # make sure res are standard AA
        for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
            if Bio.PDB.is_aa(res, standard=True):
                residues.append(res)
            else:
                sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                 (res.get_id(), os.linesep))

    return residues


def calc_ha_distance(res_a, res_b, min_dist=8.0):
    """Calculates the distance between a pair of heavy atoms
        :param res_a: Biopython residue object - residue a
        :param res_b: Biopython residue object - residue b
        :param min_dist: float - minimum distance cutoff
        :return: Minimum distance between heavy atoms
    """
    import numpy as np
    for a in res_a.get_iterator():
        for b in res_b.get_iterator():
            coord_a = a.get_coord()
            coord_b = b.get_coord()
            dist = np.linalg.norm(coord_a - coord_b)
            if dist < min_dist:
                min_dist = dist

    return min_dist


def calc_ca_distance(res_a, res_b):
    """Calculates the distance between a pair of CA atoms
    :param res_a: Biopython residue object - residue a
    :param res_b: Biopython residue object - residue b
    :return: Distance between CA atoms
    """
    import numpy as np
    a = res_a["CA"].get_coord()
    b = res_b["CA"].get_coord()
    dist = np.linalg.norm(a - b)
    return dist


def pdb_map(pdbfile, chain_ids, cutoff):
    import os
    from pandas import read_csv
    from itertools import combinations_with_replacement
    import time
    fname = "(pdb_map)"
    pdbname = os.path.splitext(pdbfile)[0]
    # output filename
    filename = pdbname + "_" + str(cutoff) + "A.txt"
    fileout = open(filename, 'w')
    # write header for output file
    fileout.write("i\tj\td\tchain_1\tchain_2\n")

    # create list of residues from pdb
    residues = get_residues(pdbfile, chain_ids=chain_ids)

    # make each possible pairs of residues
    pair_indices = combinations_with_replacement(range(len(residues)), 2)
    start_time = time.time()
    for i, j in pair_indices:
        logging.captureWarnings(True)
        res_a = residues[i]
        res_b = residues[j]
        # get chain id
        chain_a = res_a.get_parent().id
        chain_b = res_b.get_parent().id
        dist = calc_ca_distance(res_a, res_b)
        if cutoff >= dist > 0.0:
            fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i + 1, j + 1, dist, chain_a, chain_b))

    print(fname + "\tloop time: %s" % (time.time() - start_time))
    # makes a pandas dataframe
    df_pdb = read_csv(filename, delim_whitespace=True)
    df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
    df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]

    return df_pdb, df_mon, df_inter
