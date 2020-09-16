import numpy as np
import matplotlib.pyplot as plt
import os

def eigen_decomp(msa_name):
    img_dir = "C:\\Users\\kmehr\\OneDrive\\Documents\\phd_research\\images\\2020\\JUL_2020\\"
    matrix_file_1 = "results\\FNi_apc_{}.txt".format(msa_name)
    matrix_file_2 = "results\\FNi_{}.txt".format(msa_name)
    title = os.path.basename(matrix_file_2).strip('FNi_')
    A = np.loadtxt(matrix_file_1)
    B = np.loadtxt(matrix_file_2)
    eval1, evec1 = np.linalg.eig(A)
    eval2, evec2 = np.linalg.eig(B)

    eval_diag = np.diag(eval1)
    x = evec1.dot(eval_diag).dot(evec1.T)

    n = 20
    plt.scatter(range(n), eval1[:n], label="APC", edgecolors='black', marker='s')
    plt.scatter(range(n), eval2[:n], label="no-APC", marker='^')
    plt.title("First {} eigenvalues {}".format(n, title))
    plt.legend(loc='best')
    plt.grid(which='both', alpha=0.5)
    plt.show()
    imgname = "eigvals{}_{}.png".format(n, title.strip('.txt'))
    plt.savefig(img_dir + imgname, dpi=1000, bbox_inches='tight')
    plt.close()
    return x, A
