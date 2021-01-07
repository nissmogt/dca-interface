df_file = "media-3.xlsx"


def get_lengths(msa_name):
    import pandas as pd
    df_db = pd.read_excel(df_file)

    # used in query() pandas function
    pdbid = msa_name[:4]
    chain_1 = msa_name.split("_")[1]
    chain_2 = msa_name.split("_")[3]

    try:
        x = df_db.query("PDB == @pdbid and Chain_1 == @chain_1 and Chain_2 == @chain_2")
        chain_1_length = x.Sequence_Length_1.values[0]
        chain_2_length = x.Sequence_Length_2.values[0]
        return [chain_1_length, chain_2_length]
    except ValueError:
        print("{} not found in query database.".format(msa_name))
        return [1, 1]
