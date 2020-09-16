def score_vs_distance(msa, n, cutoff):
    import matplotlib.pyplot as plt
    from distance_dca import distance_dca

    # img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\AUG_2020\\"
    # imgname = "{}dca_distance_distribution_{}_top{}.png".format(img_dir, msa, n)
    distances = distance_dca(msa, n_pairs=n, cutoff=cutoff)
    # plt.figure(figsize=(8, 8))
    plt.scatter(x="d", y="score", data=distances, label=msa.strip(".fas"), s=90)
    plt.legend(loc='best')
    plt.xlabel("distance (A)")
    plt.ylabel("FN score")
    plt.show()
