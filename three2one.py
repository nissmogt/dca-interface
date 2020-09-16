def three2one(sequence):
    """ Lookup table - translate a protein sequence from 3 to 1 letter code
    """

    code = {"GLY": "G", "ALA": "A", "LEU": "L", "ILE": "I",
            "ARG": "R", "LYS": "K", "MET": "M", "CYS": "C",
            "TYR": "Y", "THR": "T", "PRO": "P", "SER": "S",
            "TRP": "W", "ASP": "D", "GLU": "E", "ASN": "N",
            "GLN": "Q", "PHE": "F", "HIS": "H", "VAL": "V",
            "M3L": "K", "MSE": "M", "CAS": "C"}

    newprot = ""
    for a in sequence:
        newprot += code.get(a, "?")

    return newprot
