import os
import shutil
import pandas as pd
import urllib.request
import tarfile

db_file = "media-3.xlsx"
df_db = pd.read_excel(db_file)
outdir = "Sequences/"

filenames = []
N = 1
for i, row in df_db.iterrows():
    if i < N:
        prefix = row["prefix"]
        pdbid = row["PDB"]
        chain_1 = row["Chain 1"]
        chain_2 = row["Chain 2"]
        filenames.append("{}_{}_{}_{}".format(pdbid, chain_1, pdbid, chain_2))

        filedir = "{}{}/".format(outdir, pdbid)
        if not os.path.exists(filedir):
            os.mkdir(filedir)

        print("Downloading {} to {}".format(filenames[i], filedir))
        url = ("https://marks.hms.harvard.edu/ecolicomplex/data/calibration/{}/{}.tar.gz".format(prefix, prefix))
        urllib.request.urlretrieve(url, "{}{}.tar.gz".format(filedir, filenames[i]))

        print("Extracting archive and saving...")
        t = tarfile.open("{}{}.tar.gz".format(filedir, filenames[i]), "r")
        msa_name = "output/{}/concatenate/{}.a2m".format(prefix, prefix)
        t.extractall(filedir, members=[t.getmember(msa_name)])
        t.close()
        os.replace(filedir+msa_name, "{}{}.a2m".format(outdir, filenames[i]))


with open("msa_names.txt", "w") as fileout:
    fileout.write("\n".join(filenames))
