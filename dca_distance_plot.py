def dca_distance_plot(msa, n, cutoff):
    import matplotlib.pyplot as plt
    from distance_dca import distance_dca

    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\AUG_2020\\"
    imgname = "{}dca_distance_distribution_{}_top{}.png".format(img_dir, msa, n)
    # n = 15
    # msa = "1EFP_A_1EFP_B.fas"
    distances = distance_dca(msa, n_pairs=n, cutoff=cutoff)
    plt.figure(figsize=(8, 8))
    plt.hist(distances['d'], edgecolor='black', color='gold')
    plt.xlabel("distance (A)")
    plt.ylabel("counts")
    plt.title("{} Top {} DCA".format(msa, n))
    plt.savefig(imgname, dpi=500, bbox_inches='tight')
    plt.show()
    plt.close()
