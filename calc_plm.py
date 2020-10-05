import numpy as np
import pandas as pd
import os
import glob
import argparse

# import matlab.engine

### DIRECTORIES ###
projectdir = '/home/kmm5/dca_interface/'
msadir = '/home/kmm5/dca_interface/PDB_benchmark_1/'
resultsdir = '/scratch/kmm5/benchmark_results_1/'
plmdir = '/home/kmm5/plm-code/'
# mfdir = '/home/kmm5/mfdca/'

### FILES ###
gene_list = '/home/kmm5/dca_interface/gene_list.csv'
complex_list = '/home/kmm5/dca_interface/complex_list.csv'


def process_id_chain_list(length_list):
    df = pd.read_csv(length_list, delimiter=',', header=0)
    return df


def process_list(gene_list, complex_list):
    df_glist = pd.read_csv(gene_list, delimiter=',')
    df_clist = pd.read_csv(complex_list, delimiter=',')
    df_geneA = df_clist.loc[:, ['gene_a']]
    df_geneB = df_clist.loc[:, ['gene_b', 'Length']]
    df_geneA = df_geneA.rename(columns={'gene_a': 'Gene'})
    df_geneB = df_geneB.rename(columns={'gene_b': 'Gene'})
    df_genes = df_glist.loc[:, ['Gene', 'Length']]

    df_geneA = df_geneA.merge(df_genes, on='Gene')
    df_geneB = df_geneB.merge(df_genes, on='Gene')
    df_geneB = df_geneB.drop(columns=['Length_x'])
    df_geneB = df_geneB.rename(columns={'Length_y': 'Length'})
    df_sgenes = df_geneA.append(df_geneB, ignore_index=True)
    df_sgenes.drop_duplicates(subset='Gene', inplace=True)
    return df_sgenes


def split_headerandseq(msa_file, len_a, limit=-1):
    """Function to split fasta headers and sequences"""
    header_a = []
    header_b = []
    seq_a = []
    seq_b = []
    lines = open(msa_file, "r")
    next(lines)  # skips null template
    next(lines)
    for idx, line in enumerate(lines):
        line = line.rstrip()  # removes whitespace from the right
        if line[0] == '>':
            if len(header_a) == -1:
                break
            header_entry = line[1:].split('_')
            if len(header_entry) < 3:
                header_a.append(header_entry[0])
                header_b.append(header_entry[1])
            else:
                header_a.append(header_entry[1])
                header_b.append(header_entry[2])
        else:
            seq_a.append(line[:len_a])  # sequence A
            seq_b.append(line[len_a:])  # sequence B

    lines.close()
    return np.array(header_a), np.array(header_b), np.array(seq_a), np.array(seq_b)


def permute_index(n_seqs, n_replicates):
    # creates 2 lists of random indices for seq A and B
    #     rindex1 = []
    #     rindex2 = []
    for seed in range(n_replicates):
        R1 = np.random.RandomState(seed)
        R2 = np.random.RandomState(seed + 2)
        yield R1.permutation(n_seqs), R2.permutation(n_seqs)


#         rindex1.append(R1.permutation(n_seqs))
#         rindex2.append(R2.permutation(n_seqs))


def scramble_sequence(msa_file, len_a, n_replicates):
    """Randomly pair sequences"""
    msa_name = msa_file.split('/')[-1]
    header_a, header_b, seq_a, seq_b = split_headerandseq(msa_file, len_a)
    n_seqs = len(seq_b)
    # creates 2 lists of random indices for seq A and B
    index = list(permute_index(n_seqs, n_replicates))
    outfile = []
    for rep in range(n_replicates):
        scramble_seq = []
        scramble_header = []
        for i in range(n_seqs):
            scramble_header.append(header_a[index[rep][0][i]] + '_' + \
                                   header_b[index[rep][1][i]])
            scramble_seq.append(seq_a[index[rep][0][i]] + seq_b[index[rep][1][i]])
        scramble_msa = dict(zip(scramble_header, scramble_seq))
        ### Write MSA replicates to file ###
        outfile.append(('MSA_rep%d_scrambled_' + msa_name) % rep)
        with open(outfile[rep], 'w') as f:
            for key in scramble_msa.keys():
                f.write(">%s\n%s\n" % (key, scramble_msa[key]))
    return outfile


def read_dca(filein):
    ###Function that reads in the couplings and fields Matlab matrix file.
    import h5py
    mat = {}
    f = h5py.File(filein, 'r')
    for k, v in f.items():
        mat[k] = np.array(v)
    h = mat['h']
    J = mat['J']
    q = mat['h'].shape[1]
    N = mat['h'].shape[0]
    fields = np.zeros((N, q))
    # couplings = np.swapaxes(J,1,2)
    couplings = J
    fields = h
    f.close()
    return fields, couplings


## AVERAGE MATRICES
def avg_scrmbl(n_replicates, out_scrmbl_matrix):
    ###Function that averages over coupling matrices.
    # initialize coupling sum to zero
    j_sum = 0.0
    for i in range(n_replicates):
        # filein = ("scrmbl_rep%d_matrix_NUOJ_NUOK_plmdca_rt2.mat" % i)
        h_scrmbl, j_scrmbl = read_dca(out_scrmbl_matrix[i])
        j_sum += j_scrmbl
    j_scrmbl_avg = j_sum / n_replicates
    return j_scrmbl_avg


## Calculate Scores
def calc_score(msa_name, j_paired, j_scrmbl_avg):
    ###Function that calculates Frobenius Norm score of coupling matrices.'''
    q = j_paired.shape[0]
    N = j_scrmbl_avg.shape[-1]
    norms_a = np.zeros((N, N))
    norms_b = np.zeros((N, N))
    x = np.zeros((N, N))
    print("Size %d" % N)

    for i in range(N - 1):
        for j in range(i + 1, N):
            #         print("Row:%d Col:%d" % (i,j))
            norms_a[i, j] = np.linalg.norm(
                (j_paired[:q, :q, i, j]), 'fro')
            norms_b[i, j] = np.linalg.norm(
                (j_scrmbl_avg[:q, :q, i, j]), 'fro')
            # Take the difference between paired and scrambled FN
            x[i, j] = (norms_a[i, j] - norms_b[i, j])
            x[j, i] = x[i, j]
    #         print(norms_a[i,j], norms_b[i,j])
    #         print("Xij=%f" % x[i,j])
    x_i = np.mean(x, axis=0) * N / (N - 1)
    x_j = np.mean(x, axis=1) * N / (N - 1)
    x_mean = np.mean(x) * N / (N - 1)
    ###APC correction
    corrnorms = x - np.outer(x_i, x_j) / x_mean
    xsquared = x * x - np.outer(x_i, x_j) / x_mean
    np.savetxt('FNi_' + msa_name + '.txt', x, delimiter='\t')
    np.savetxt('FNi_apc_' + msa_name + '.txt', corrnorms, delimiter='\t')
    np.savetxt('FNi_xsq_' + msa_name + '.txt', xsquared, delimiter='\t')
    return x_mean, corrnorms, xsquared


### Calculate KLD-interface ###
# def kld_interface(pij_paired, pij_scrambled, kld_out, kld_apc):
#    len_complex = 

# def filter_msa(msa_files, len_limit, n_limit):
#    from Bio import AlignIO
#    msa_list = []
#    for i, msa in enumerate(msa_files):
#        filein = open(msa)
#        alignment = AlignIO.read(filein, 'fasta')
#        align_len = alignment.get_alignment_length()
#        num_seqs = len([1 for line in open(msa) if line.startswith(">")])
#        if num_seqs > n_limit:
#            if align_len < len_limit:
#                msa_list.append(msa)
#        filein.close()
#    return msa_list


### MAIN ###
def main():
    ### HARD-CODED VALUES ###
    PC = 0.2  # pseudocount
    NR = 5  # number of replicates
    ### HEADER ###
    TITLE = "|| WELCOME TO THE DCA-INTERFACE PROGRAM. ||"
    width = 43
    DASH = "-" * width
    print(DASH)
    print(TITLE)
    print(DASH)
    print("\nCurrent working directory is %s" % os.getcwd())

    # try:
    #     os.mkdir(resultsdir)
    # except OSError:
    #     print("\tCreation of the directory %s failed" % resultsdir)
    # else:
    #     print("\tSuccessfully created the directory %s " % resultsdir)
    #
    # os.chdir(resultsdir)
    # print("\tCurrent working directory is:", os.getcwd())
    ### START MATLAB ENGINE VROOM! ###
    # eng = matlab.engine.start_matlab()
    # eng.addpath('%s' % resultsdir, nargout=0)
    # eng.addpath('%s' % mfdir, nargout=0)
    # eng.addpath('%s' % plmdir, nargout=0)
    # eng.addpath('%s' % msadir, nargout=0)
    # eng.addpath('%s/functions' % plmdir, nargout=0)
    # eng.addpath('%s/3rd_party_code/minFunc' % plmdir, nargout=0)
    ### PARSE NUMBER OF CORES FROM SLURM ###
    # parser = argparse.ArgumentParser(description='Calcs plm')
    # parser.add_argument('nc', help='Number of cores to run plm code')
    # args = parser.parse_args()
    # nr_of_cores = args.nc

    # Process gene list and complex list to get len A
    # len_limit = 400
    # n_limit = 800
    # msa_files = filter_msa(msa_list, len_limit, n_limit)

    # make a dataframe of every gene's length and drop duplicates
    # df_sgenes = process_list(gene_list, complex_list)
    msa_list = glob.glob("PDB_benchmark_alignments\\*.fas")
    msa_files = msa_list
    df_id_list = process_id_chain_list("lengths.csv")

    # header_items = ["Name", "Length"]
    # column_width = width // len(header_items)
    # column_width_items = [column_width] * len(header_items)
    # header_fmt_content = [None] * (len(header_items) + len(column_width_items))
    # header_fmt_content[::2] = header_items
    # header_fmt_content[1::2] = column_width_items
    # header = "{:^{}s}{:^{}s}".format(*header_fmt_content)
    # print(DASH)
    # print(header)
    # print(DASH)

    # MAIN LOOP ###
    # for msa_in in msa_files:
    #     msa_ext = '.fas'
    #     msa_name = os.path.basename(msa_in).strip(msa_ext)
    #     prot_len = int(df_id_list[df_id_list['PDBid'] == msa_name].iloc[0]['Length'])

        # gene_a = (msa_name.split('_')[0])
        # prot_len = find_length().get(str(msa_name.split('_')[0]))
        # prot_len = int(df_sgenes[df_sgenes['Gene']==gene_a].iloc[0]['Length'])

        # paired output filenames 
        # out_paired_fn = (resultsdir+'fn_%s_plmdca_rt%d.txt' \
        #         % (msa_name,PC*10))
        # out_paired_dist = (resultsdir+'pdist_%s_plmdca_rt%d.mat' \
        #         % (msa_name,PC*10))
        # out_paired_mat = (resultsdir+'matrix_%s_plmdca_rt%d.mat' \
        #         % (msa_name,PC*10))
        #
        # print('\t%s\t%s\t' % (msa_proc, prot_len))


       # scrmbl_msa = scramble_sequence(msa_in, prot_len, NR)

# RUN DCA MATLAB CODE FOR PAIRED MSA ###
# PLM-DCA
# eng.plmDCA_symmetric_mod7(msa_proc,out_paired_fn,0.01,0.01,PC,
#         nr_of_cores,out_paired_dist,out_paired_mat, nargout=0)
# MF-DCA
# eng.dca_pairdist(msa_proc, out_paired_fn, out_paired_mat, PC, nargout=0)

#        print('\nStarting Replicate loop...')

# RUN DCA MATLAB CODE FOR SCRAMBLED MSA ###
#        out_scrmbl_mat = []
#        for i in range(NR):
#            print('Wrote %s to %s dir...' % (scrmbl_msa[i], resultsdir))
# PLM-DCA
#            out_scrmbl_fn = (resultsdir+'scrmbl_rep%d_fn_%s_plmdca_rt%d.txt' \
#                    % (i, msa_name,PC*10))
#            out_scrmbl_dist = (resultsdir+'scrmbl_rep%d_pdist_%s_plmdca_rt%d.mat' \
#                    % (i, msa_name,PC*10))
#            out_scrmbl_mat.append(resultsdir+'scrmbl_rep%d_matrix_%s_plmdca_rt%d.mat' \
#                    % (i, msa_name,PC*10))
#            eng.plmDCA_symmetric_mod7(scrmbl_msa[i],out_scrmbl_fn,0.01,0.01,PC,
#                    nr_of_cores,out_paired_dist,out_scrmbl_mat[i], nargout=0)
# MF-DCA
# out_scrmbl_fn = (resultsdir+'scrmbl_rep%d_fn_%s_plmdca_rt%d.txt' \
#        % (i, msa_name,PC*10))
# out_scrmbl_dist = (resultsdir+'scrmbl_rep%d_pdist_%s_plmdca_rt%d.mat' \
#        % (i, msa_name,PC*10))
# out_scrmbl_mat.append(resultsdir+'scrmbl_rep%d_matrix_%s_plmdca_rt%d.mat' \
#        % (i, msa_name,PC*10))
# eng.dca_pairdist(scrmbl_msa[i], out_scrmbl_fn, out_scrmbl_mat[i], PC, nargout=0)

#        j_scrmbl_avg = avg_scrmbl(NR, out_scrmbl_mat)
#        h_paired, j_paired = read_dca(out_paired_mat)
#        x_mean, corrnorms, xsquared = calc_score(msa_name, j_paired, j_scrmbl_avg)

# eng.quit()
if __name__ == '__main__':
    main()
