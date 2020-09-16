def score_compare(df_dca, df_gremlin, n_pairs, msa_name, title):
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\JUL_2020\\"
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    match = pd.merge(df_gremlin[:n_pairs], df_dca[:n_pairs], on=['i', 'j'], how='inner')
    # dca_array = df_dca.to_numpy()
    # g_array = df_gremlin.to_numpy()
    # index_g = []
    # for i in range(n_pairs):
    #     for j in range(400):
    #         if dca_array[j][0] == g_array[i][0] and dca_array[j][1] == g_array[i][1]:
    #             index_g.append(j)
    fig, ax = plt.subplots()
    plt.scatter(x='score_x', y='score_y', data=match, s=30, edgecolors='black', c=range(len(match)), cmap='tab20c')
    # plt.scatter(range(len(index_g)), index_g, edgecolors='black')
    # plt.plot(range(len(index_g)), range(len(index_g)), c='black', linestyle='dashed')
    plt.xlabel('Gremlin (webserver) score')
    # plt.ylabel('Gremlin (local-svm prior) score')
    plt.ylabel('DCAi score')
    plt.title("{} {}".format(msa_name, n_pairs))
    # plt.minorticks_on()
    plt.grid(axis='both', alpha=0.3)
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.show()
    imgname = "{}_{}_top{}_tab20c.png".format(msa_name, title, len(match))
    plt.savefig(img_dir + imgname, dpi=1000, bbox_inches='tight')
    plt.close()
