def draw_sasa_res(msaName):
    """
    Draw residues with SASA >= cutoff in TCL script.
    :param sasa_file:   GetArea Server SASA file
    :param cutoff: SASA ratio threshold

    """
    import os
    import pandas as pd
    sasa_file = "sasa\\total_{}_freesasa.txt".format(msaName)
    # msaName = os.path.basename(sasa_file).strip(".txt")
    pdbid = msaName[:4]
    chains = [msaName.split("_")[1], msaName.split("_")[3]]
    output_dir = os.path.dirname(sasa_file)
    threshold = 50
    df = pd.read_csv(sasa_file, delimiter='\t', usecols=(1, 2, 7))
    df = df.dropna().reset_index(drop=True)
    output_filename = "tcl_scripts\\draw_sasa\\draw_totalSASA_{}_thresh_20_50.tcl".format(msaName)
    target = open("{}".format(output_filename), 'w')
    target.write("color Display Background gray\n")
    target.write("display projection Orthographic\n")
    target.write("axes location Off\n")
    target.write("display rendermode GLSL\n")
    target.write("display ambientocclusion on\n")
    target.write("display depthcue off\n")
    for ch in chains:
        print("Chain {}".format(ch))
        target.write("mol representation Surf\n")
        target.write("mol material Diffuse\n")
        df_chain = df[df["Chain"] == ch].reset_index(drop=True)
        nRes = len(df_chain)
        atomselect = 0
        nSasa_res = 0
        outside = ''
        middle = ''
        inside = ''
        for i in range(nRes):
            resnum = int(df_chain["Residue_num"][i])
            ratio = df_chain["Ratio"][i]
            if ratio >= threshold:
                nSasa_res += 1
                outside += str(resnum) + ' '
            elif threshold > ratio >= 20:
                middle += str(resnum) + ' '
            elif ratio < 20:
                inside += str(resnum) + ' '
        if ch == chains[0]:
            target.write("mol color ColorID 22\n")
        else:
            target.write("mol color ColorID 23\n")
        target.write("mol selection {resid %s and chain %s}\n" % (outside, ch))
        target.write("mol addrep top\n")
        # middle sasa values
        if ch == chains[0]:
            target.write("mol color ColorID 9\n")
        else:
            target.write("mol color ColorID 8\n")
        target.write("mol selection {resid %s and chain %s}\n" % (middle, ch))
        target.write("mol addrep top\n")
        # buried sasa values
        if len(inside) > 0:     # sometimes this list is empty so gotta check
            if ch == chains[0]:
                target.write("mol material EdgyGlass\n")
                target.write("mol color ColorID 12\n")
            else:
                target.write("mol material EdgyGlass\n")
                target.write("mol color ColorID 17\n")
            target.write("mol selection {resid %s and chain %s}\n" % (inside, ch))
            target.write("mol addrep top\n")
        fractionSasa = 100*(nSasa_res/nRes)
        print("Percent of residues outside: {}".format(fractionSasa))
        print("Wrote tcl file: %s" % output_dir + output_filename)
    target.write("mol modselect 0 top \"chain {} {}\"\n".format(chains[0], chains[1]))
    target.write("mol modmaterial 0 top AOEdgy\n")
    target.write("mol modcolor 0 top Secondary Structure\n")
    target.write("mol modstyle 0 top NewCartoon\n")
    target.write("display resetview\n")
    target.close()

# draw_sasa_res("1EM8_D_1EM8_C")