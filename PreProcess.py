import os
import sys
import glob
import logging


class Preprocess:

    def __init__(self, msa_name, cutoff):

        self.msa_name = msa_name
        self.cutoff = cutoff

        self.chain_ids = [(self.msa_name.strip(".fas")).split("_")[1], (self.msa_name.strip(".fas")).split('_')[3]]
        self.pdbid = self.msa_name[:4].lower()
        self.pdb_path = "PDB_benchmark_structures\\"
        self.pdbfile = "{}.pdb".format(self.pdbid)
        # self.msa_path = os.path.dirname(msa_file) + "\\"
        self.msa_path = "PDB_benchmark_alignments\\"

        self.chain_lengths = []
        self.length_a = 0

        self.result_path = ""

        alphabet = "ARNDCQEGHILKMFPSTWYV-"
        self.states = len(alphabet)
        self.a2n = {}
        for a, n in zip(alphabet, range(self.states)):
            self.a2n[a] = n

        ################

    def aa2num(self, aa):
        """convert aa into num"""
        if aa in self.a2n:
            return self.a2n[aa]
        else:
            return self.a2n['-']

    def read_pdb(self):
        """
        Reads PDB file and extracts relevant chain sequence.
        :return: Two lists - Fasta header list and PDB sequence list. len == 2
        """
        fname = "(read_pdb)"
        pdb_fasta_file = "{}{}.fasta".format(self.pdb_path, self.msa_name[:4])
        pdb_fasta_obj = open(pdb_fasta_file, 'r')
        num_of_chains = int(pdb_fasta_obj.read().count('>'))
        pdb_fasta_obj.seek(0)  # return to top of file

        pdb_fasta_header = []
        pdb_fasta_seq = []
        # Loop thru each chain and append header and seq to a list, then make into lookup table
        for line in range(num_of_chains):
            pdb_fasta_header.append(pdb_fasta_obj.readline().rstrip()[-1])
            pdb_fasta_seq.append(pdb_fasta_obj.readline().rstrip())
        logging.debug("{}\t\t{}\tnum of chains: {}".format(fname, pdb_fasta_file, num_of_chains))
        pdb_fasta_obj.close()

        # Error check
        if len(pdb_fasta_seq) != len(pdb_fasta_header):
            print("Error! Number of seqs not equal to number of chains\n".format(fname))
        return pdb_fasta_header, pdb_fasta_seq

    def read_length_file(self):
        import csv
        print("\treading lengths...")
        with open('lengths.csv', 'r', newline='\n') as csvfile:
            r = csv.reader(csvfile)
            for row in r:
                if self.msa_name.strip(".fas") == row[0]:
                    self.length_a = int(row[-1])
        print("\tlength of first chain: {}".format(self.length_a))
        return self.length_a

    def msa_template(self, split=False, len_a=None):
        """
        Reads first sequence in MSA which should be a template
        :return: List - MSA sequence
        """
        # Create MSA-seq-template file object
        msa_file_obj = open("{}{}.fas".format(self.msa_path, self.msa_name), 'r')
        msa_header = msa_file_obj.readline().rstrip()
        msa_seq = msa_file_obj.readline().rstrip()
        msa_file_obj.close()
        if split:
            return msa_seq[:len_a], msa_seq[len_a:]
        else:
            return msa_seq

    def parse_fasta(self, null=True, limit=-1):
        """
        function to parse fasta
        :return header and sequence of fasta
        """
        import numpy as np
        header = []
        sequence = []
        lines = open("{}{}".format(self.msa_path, self.msa_name), 'r')
        if null:
            # used to skip first fasta sequence and header
            skip_null = [next(lines) for x in range(2)]
        for line in lines:
            line = line.rstrip()
            if line[0] == ">":
                if len(header) == limit:
                    break
                header.append(line[1:])
                sequence.append([])
            else:
                sequence[-1].append(line)
        lines.close()
        sequence = [''.join(seq) for seq in sequence]
        return np.array(header), np.array(sequence)

    def filter_gaps(self, msa, gap_cutoff=0.5):
        """filters alignment to remove gappy positions"""
        import numpy as np
        tmp = (msa == self.states - 1).astype(np.float)
        non_gaps = np.where(np.sum(tmp.T, -1).T / msa.shape[0] < gap_cutoff)[0]
        return msa[:, non_gaps], non_gaps

    def get_eff(self, msa, eff_cutoff=0.8):
        '''compute effective weight for each sequence'''
        from scipy.spatial.distance import pdist, squareform
        import numpy as np
        ncol = msa.shape[1]

        # pairwise identity
        msa_sm = 1.0 - squareform(pdist(msa, "hamming"))

        # weight for each sequence
        msa_w = (msa_sm >= eff_cutoff).astype(np.float)
        msa_w = 1 / np.sum(msa_w, -1)

        return msa_w

    def mk_msa(self, seqs):
        """converts list of sequences to msa"""
        import numpy as np
        msa_ori = []
        for seq in seqs:
            msa_ori.append([self.aa2num(aa) for aa in seq])
        msa_ori = np.array(msa_ori)

        # remove positions with more than > 50% gaps
        # msa, v_idx = self.filt_gaps(msa_ori, 0.5)
        msa = msa_ori

        # compute effective weight for each sequence
        msa_weights = self.get_eff(msa, 0.8)

        # compute effective number of sequences
        ncol = msa.shape[1]  # length of sequence
        # w_idx = v_idx[np.stack(np.triu_indices(ncol, 1), -1)]

        return {"msa_ori": msa_ori,
                "msa": msa,
                "weights": msa_weights,
                "neff": np.sum(msa_weights),
                "nrow": msa.shape[0],
                "ncol": ncol,
                "ncol_ori": msa_ori.shape[1]}
        # "v_idx": v_idx,
        # "w_idx": w_idx
        # }

    def combine_chain_sequence(self, split=False):
        """
        Concatenates PDB sequences
        :return: str - combined chain 1 and chain 2 protein sequence from pdb
        """
        fname = "(combine_pdb_seq)"
        print("{}\tcat-ting sequences...".format(fname))
        # writes file of pdbid and length of first chain and pdb fasta file
        # Create PDB fasta file object
        logging.info(fname + "\t\tFasta file\tNumber of chains")
        pdb_fasta_header, pdb_fasta_seq = self.read_pdb()
        header_seq_dict = dict(zip(pdb_fasta_header, pdb_fasta_seq))

        # For each chain in msa find relevant pdb seq and cat the two seqs
        first_seq = header_seq_dict[self.chain_ids[0]]
        second_seq = header_seq_dict[self.chain_ids[1]]
        # self.chain_lengths = [len(first_seq), len(second_seq)]
        full_seq = first_seq + second_seq
        logging.debug("{}\t\tPDB seq length: {}".format(fname, len(full_seq)))
        logging.info("\tFinished cat-ting sequences.".format(fname))
        if split:
            return first_seq, second_seq
        else:
            return full_seq


    def distance_matrix(self, all_atom=False):
        """
        Calculates distance matrix.
        :param all_atom:
        :return: Three Dataframe objects.
        """
        from pandas import read_csv
        from itertools import combinations_with_replacement
        import time
        fname = "(pdb_map)"

        if all_atom:
            filename = "{}heavy_atom_distance_matrix_{}_{}A.txt".format(self.pdb_path, self.pdbfile.strip(".cif"),
                                                                        self.cutoff)
        else:
            filename = "{}ca_distance_matrix_{}_{}A.txt".format(self.pdb_path, self.pdbfile.strip(".cif"), self.cutoff)

        fileout = open(filename, 'w')
        fileout.write("i\tj\td\tchain_1\tchain_2\n")
        # if ".cif" in self.pdbfile:
            # output filename
            # filename = "{}distance_matrix_{}_{}A.txt".format(self.pdb_path, self.pdbfile.strip(".cif"), self.cutoff)
        # else:
        #     output filename
            # filename = "{}distance_matrix_{}_{}A.txt".format(self.pdb_path, self.pdbfile.strip(".pdb"), self.cutoff)

        # output list of residues from pdb
        residues = self.get_residues()
        # make each possible pairs of residues
        pair_indices = combinations_with_replacement(range(len(residues)), 2)
        start_time = time.time()
        for i, j in pair_indices:
            res_a = residues[i]
            res_b = residues[j]
            # get chain id
            if all_atom:
                # if res_a.get_id()[1] - res_b.get_id()[1] > 4:
                chain_a = res_a.get_parent().id
                chain_b = res_b.get_parent().id
                mindist = self.calc_min_dist(res_a, res_b)
                if mindist <= self.cutoff:
                    fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i, j, mindist, chain_a, chain_b))
                else:
                    if res_a.has_id("CA") and res_b.has_id("CA"):
                        chain_a = res_a.get_parent().id
                        chain_b = res_b.get_parent().id
                        dist = self.calc_ca_distance(res_a, res_b)
                        if self.cutoff >= dist > 0.0:
                            fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i, j, dist, chain_a, chain_b))
                    else:
                        print("{} NOTE! Res {} \n\tor {} not calculated! (missing CA)\n".format(fname, res_a.get_full_id(),
                                                                                                res_b.get_full_id()))
        fileout.close()
        print("{}\t -- MAIN LOOP TIME -- {}".format(fname, time.time() - start_time))
        # makes a pandas dataframe

        df_pdb = read_csv(filename, delim_whitespace=True)
        df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
        df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
        return df_pdb, df_mon, df_inter

    def read_distance_matrix_file(self, all_atom=False):
        from pandas import read_csv
        pdbid = self.pdbfile.strip(".cif")
        if all_atom:
            filename = "{}heavy_atom_distance_matrix_{}_{}A.txt".format(self.pdb_path, pdbid, self.cutoff)
        else:
            filename = "{}ca_distance_matrix_{}_{}A.txt".format(self.pdb_path, pdbid, self.cutoff)
        df_pdb = read_csv(filename, delim_whitespace=True)
        df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
        df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]
        return df_mon, df_inter

    def calc_min_dist(self, res_a, res_b):
        import numpy as np
        dist = []
        for atom_i in res_a:
            for atom_j in res_b:
                dist.append(np.linalg.norm(atom_i.get_coord() - atom_j.get_coord()))
        dist_array = np.array(dist)
        mindist = min(dist_array[np.nonzero(dist_array)])
        return mindist

    def calc_ca_distance(self, res_a, res_b):
        """
        Calculates the distance between a pair of CA atoms
        :param res_a: Biopython residue object - residue a
        :param res_b: Biopython residue object - residue b
        :return: Distance between CA atoms
        """
        import numpy as np
        a = res_a["CA"].get_coord()
        b = res_b["CA"].get_coord()
        dist = np.linalg.norm(a - b)
        return dist

    def get_residues(self, seq=False):
        """
        Build a simple list of residues from a single chain of a PDB file.
        :param seq:
        :return: A list of Bio.PDB.Residue objects.
        """
        import Bio.PDB
        fname = "(get_residues)"
        pdb_id = self.pdbfile.strip('.pdb')
        print('\t{}\tprocessing {} file...'.format(fname, self.pdbfile))
        parser = Bio.PDB.PDBParser()

        struct = parser.get_structure(pdb_id, self.pdb_path + self.pdbfile)
        model = struct[0]
        # if len(self.chain_ids) == 0:
        # get residues from every chain.
        #    chains = model.get_list()
        # else:
        chains = [model[ch_id] for ch_id in self.chain_ids]
        print("\t{} CHAIN IDs:\t{}".format(fname, self.chain_ids))

        residues = []
        sequence = []
        for ch in chains:
            # make sure res are standard AA
            num_residues = 0
            for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
                # if Bio.PDB.is_aa(res, standard=True):
                is_regular_res = res.has_id('CA') and res.has_id('O')
                res_id = res.get_id()[0]
                if (res_id == ' ' or res_id == 'H_MSE' or res_id == 'H_M3L' or res_id == 'H_CAS') and is_regular_res:
                    residues.append(res)
                    sequence.append(res.get_resname())
                    num_residues += 1
                else:
                    sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                     (res.get_id(), os.linesep))
            self.chain_lengths.append(num_residues)  # P1: BUG if runs twice, chain_length list will be twice as long.

        if seq:
            sequence = self.three2one(sequence)
            seq_a = sequence[:self.chain_lengths[0]]
            seq_b = sequence[self.chain_lengths[0]:]
            return seq_a, seq_b
        else:
            return residues

    def three2one(self, prot):
        """ Lookup table - translate a protein sequence from 3 to 1 letter code
        """

        code = {"GLY": "G", "ALA": "A", "LEU": "L", "ILE": "I",
                "ARG": "R", "LYS": "K", "MET": "M", "CYS": "C",
                "TYR": "Y", "THR": "T", "PRO": "P", "SER": "S",
                "TRP": "W", "ASP": "D", "GLU": "E", "ASN": "N",
                "GLN": "Q", "PHE": "F", "HIS": "H", "VAL": "V",
                "M3L": "K", "MSE": "M", "CAS": "C"}

        newprot = ""
        for a in prot:
            newprot += code.get(a, "?")

        return newprot

    def plot_pdb_map(self, read=None, calc=None, heatmap=None, all_atom=False):
        """
        Used for plotting a pdb distance map. For testing and debugging purposes only.
        :param all_atom:
        :param calc:
        :param read:
        :param heatmap:
        :return:
        """
        import matplotlib.pylab as plt

        print("starting pdb_map...")
        if read:
            df_mon, df_inter = self.read_distance_matrix_file()
        if calc:
            if all_atom:
                df_pdb, df_mon, df_inter = self.distance_matrix(all_atom=True)
            else:
                df_pdb, df_mon, df_inter = self.distance_matrix()

        print("plotting...")
        # Plotting
        fig = plt.figure(figsize=(10, 10), dpi=100)
        ax = fig.add_subplot(1, 1, 1)

        # monomer
        ax.scatter('i', 'j', data=df_mon, label='PDB monomer', c='xkcd:navy', cmap='coolwarm', marker='s')
        # interface
        ax.scatter('i', 'j', data=df_inter, label='PDB dimer', c='olive', cmap='Viridis', marker='s')

        if heatmap:
            # monomer
            ax.scatter('i', 'j', data=df_mon, label='PDB monomer', c=df_mon['distance'], marker='s')
            # interface
            ax.scatter('i', 'j', data=df_inter, label='PDB dimer', c=df_inter['distance'], marker='s')

        # Plot dimer separator line
        length = self.chain_lengths[0] + self.chain_lengths[1]
        plt.hlines(self.chain_lengths[0], 0, length, linestyles='dashed', alpha=0.6)
        plt.vlines(self.chain_lengths[0], 0, length, linestyles='dashed', alpha=0.6)

        # plot design
        ax.legend(loc='lower right')
        plt.minorticks_on()
        plt.grid(alpha=0.3)
        plt.xlabel("residue i"), plt.ylabel("residue j")
        ax.grid(which='major', alpha=0.4, c='black')
        ax.grid(which='minor', linestyle=':', alpha=0.5, c='gray')
        plt.show()

    def read_dca_matrix(self, score_matrix):
        """
        Makes a contact map from a DCA score matrix file.
        :param score_matrix: Frobenius norm matrix
        :return: Three-column dataframe composed of pair i, j, and fn score
        """
        fname = "(read_dca_matrix)"
        import numpy as np
        dca_matrix = np.loadtxt(score_matrix)
        filename = score_matrix.strip(".txt")
        logging.info("{}\t\tScore matrix filename: {}".format(fname, filename))

        x_output = []
        n = dca_matrix.shape[0]
        logging.debug("{}\tshape of matrix: {}".format(fname, n))
        # index starts from zero
        for i in range(n - 1):
            for j in range(i + 1, n):
                x_output.append([i, j, dca_matrix[i, j]])
        dca_array = np.array(x_output)
        return dca_array



    def apply_map_g(self, dca_array, map_dictionary):
        import numpy as np
        print("(apply_map)")
        map_dca_list = []
        for i, j, score, i_score in dca_array:
            if int(i) in map_dictionary.keys() and int(j) in map_dictionary.keys():
                map_index_i = map_dictionary[int(i)]
                map_index_j = map_dictionary[int(j)]
                # map_dca_list.append([map_index_i, map_index_j, score, i_score])
                map_dca_list.append([map_index_i, map_index_j, score, i_score])
        return np.array(map_dca_list)

    def split_header_and_seq(self, limit=-1):
        """Function to split fasta headers and sequences"""
        import numpy as np
        header_a = []
        header_b = []
        seq_a = []
        seq_b = []
        lines = open(self.msa_file, "r")
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
                seq_a.append(line[:self.chain_lengths[0]])  # sequence A
                seq_b.append(line[self.chain_lengths[0]:])  # sequence B
        lines.close()
        return np.array(header_a), np.array(header_b), np.array(seq_a), np.array(seq_b)

    def permute_index(self, n_seqs, n_replicates):
        # creates 2 lists of random indices for seq A and B
        import numpy as np
        for seed in range(n_replicates):
            R1 = np.random.RandomState(seed)
            R2 = np.random.RandomState(seed + 2)
            yield R1.permutation(n_seqs), R2.permutation(n_seqs)

    def scramble_sequence(self, n_replicates):
        """Randomly pair sequences"""
        msa_name = self.msa_name
        header_a, header_b, seq_a, seq_b = self.split_header_and_seq()
        n_seqs = len(seq_b)
        # creates 2 lists of random indices for seq A and B
        index = list(self.permute_index(n_seqs, n_replicates))
        outfile = []
        for rep in range(n_replicates):
            scramble_seq = []
            scramble_header = []
            for i in range(n_seqs):
                scramble_header.append(header_a[index[rep][0][i]] + '_' + \
                                       header_b[index[rep][1][i]])
                scramble_seq.append(seq_a[index[rep][0][i]] + seq_b[index[rep][1][i]])
            scramble_msa = dict(zip(scramble_header, scramble_seq))
            # write MSA replicates to file
            outfile.append("MSA_rep{}_scrambled_{}".format(rep, self.msa_name))
            with open(outfile[rep], 'w') as f:
                for key in scramble_msa.keys():
                    f.write(">%s\n%s\n" % (key, scramble_msa[key]))
        return outfile
