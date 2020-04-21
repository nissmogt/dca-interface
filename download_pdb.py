def get_pdb(pdb_list):
    from Bio.PDB import PDBList

    out_dir = "PDB_benchmark_structures\\"
    pdb = pdb_list
    number_ids = len(pdb)

    print("Downloading in %s:\n" % out_dir)
    for ids in pdb:
        print('%s' % ids)
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(ids, pdir=out_dir)
    return pdbl


def get_lengths(pdbid_list, id_chain_dict):
    from Bio.PDB import MMCIFParser, MMCIF2Dict
    import glob
    import csv
    import time
    directory = 'PDB_benchmark_structures\\'

    for id in pdbid_list:
        start_time = time.time()
        parser = MMCIFParser()
        m = MMCIF2Dict.MMCIF2Dict(directory + id + '.cif')
        chains = m['_entity_poly.pdbx_strand_id']
        if '_entity_poly.pdbx_seq_one_letter_code' in m.keys():
            full_structure = m['_entity_poly.pdbx_seq_one_letter_code']
            for c in chains:
                print('Chain %s' % (c))
                print('Sequence: %s' % (full_structure[chains.index(c)]))

    print("Total time: ", time.time() - start_time)


def make_list(directory):
    """
    Makes a list of pdbids and a dictionary of ids and chains

    :param directory: Directory of filename: ID1_chain1_ID2_chain2.fas
    :return: List of pdbids and a dictionary file storing ids and chains
    """
    import glob
    list_fasta_files = glob.glob(directory + '*.fas')
    list_pdb_ids = []
    chain_list = []
    for ids in list_fasta_files:
        fasta_names = ids.split('\\')[-1]
        name_chain = fasta_names.split('.fas')[0]
        list_pdb_ids.append(name_chain.split('_')[0])
        chain_list.append([name_chain.split('_')[1], name_chain.split('_')[3]])
    # Create a dictionary of PDBids and corresponding chains and write to csv file
    dic_id_chain = dict(zip(list_pdb_ids, chain_list))
    with open("PDB_benchmark_structures\\id_chains.csv", 'w') as f:
        for key in dic_id_chain.keys():
            f.write("%s,%s\n" % (key, dic_id_chain[key]))
    return list_pdb_ids, dic_id_chain


pdbid_list, id_chain_dict = make_list('PDB_benchmark_alignments\\')
# get_pdb(list_id)
# get_lengths(pdbid_list, id_chain_dict)
from Bio.PDB import MMCIFParser, MMCIF2Dict
import time

directory = 'PDB_benchmark_structures\\'
out = open('PDB_benchmark_alignments\\' + 'lengths.csv', 'w')
out.write("PDBid,Chain,Length\n")
start_time = time.time()
for ids in pdbid_list:
    parser = MMCIFParser()
    m = MMCIF2Dict.MMCIF2Dict(directory + ids + '.cif')
    chains = m['_entity_poly.pdbx_strand_id']
    if '_entity_poly.pdbx_seq_one_letter_code' in m.keys():
        full_structure = m['_entity_poly.pdbx_seq_one_letter_code']
        chain_index = id_chain_dict[ids][0]
        for c in chains:
            if chain_index in c.split(','):
                length = len(full_structure[chains.index(c)].replace('\n', ''))
                out.write("%s,%s,%s\n" % (ids, chain_index, length))
                print('ID: %s Chain %s' % (ids, chain_index))
                print('Chain %s length: %s' % (chain_index, length))
out.close()
print("Total time: ", time.time() - start_time)
