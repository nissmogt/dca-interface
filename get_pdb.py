def get_pdb(pdb_list):
    import os
    from Bio.PDB import PDBList

    out_dir = "PDB_benchmark_structures\\"
    pdb = pdb_list
    number_ids = len(pdb)

    filename = []
    not_found = []
    print("Downloading in %s:\n" % out_dir)
    for i, pdbid in enumerate(pdb):
        print('%s' % pdbid[:4])
        pdbl = PDBList()
        try:
            if not os.path.exists("{}{}.pdb".format(out_dir, pdbid)):
                x = pdbl.retrieve_pdb_file(pdbid[:4], file_format='pdb', pdir=out_dir)
                filename.append(x)

        except FileNotFoundError:
            not_found.append(pdbid)
            print("(NOTE) {} not found.".format(pdbid))
    return filename


import pandas as pd
import os

db_file = "media-3.xlsx"
df_db = pd.read_excel(db_file)

pdb_ls = []
for ids in df_db["PDB"]:
    if ids not in pdb_ls:
        pdb_ls.append(ids)

# f = get_pdb(pdb_ls)
# nf = []
# for i, record in enumerate(f):
#     pdbid = os.path.basename(record).strip(".ent")[-4:]
#     new_filename = "{}\\{}.pdb".format(os.path.dirname(record), pdbid.upper())
#     try:
#         print("Renaming to {}".format(new_filename))
#         os.replace(record, new_filename)
#     except FileNotFoundError:
#         nf.append(pdbid)
#         print("(NOTE) {} not found.".format(pdbid))


