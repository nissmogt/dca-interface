import os
import sys
import glob
import logging


class Preprocess:

    def __init__(self, pdbfile, msa_name, cutoff):
        self.pdbfile = pdbfile
        self.msa_name = msa_name  # ID_CHAIN1_ID_CHAIN2.fas
        self.cutoff = cutoff

        self.chain_ids = [(msa_name.strip(".fas")).split("_")[1], (msa_name.strip(".fas")).split('_')[3]]
        self.pdb_path = "PDB_benchmark_structures\\"
        self.msa_path = "PDB_benchmark_alignments\\"
        self.chain_lengths = []

        self.result_path = ""
        self.dca_array = 0
        self.map_dca_array = 0

    def read_pdb(self):
        """
        Reads PDB file and extracts relevant chain sequence.
        :return: Two lists - Fasta header list and PDB sequence list. len == 2
        """
        fname = "(read_pdb)"
        pdb_fasta_file = "{}{}.fasta".format(self.pdb_path, self.msa_name[:4])
        pdb_fasta_obj = open(pdb_fasta_file, 'r')
        num_of_chains = int(pdb_fasta_obj.read().count('\n') / 2)
        pdb_fasta_obj.seek(0)  # return to top of file

        assert num_of_chains == len(self.chain_ids)  # asserts the chain count is correct
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

    def read_msa(self):
        # Create MSA-seq-template file object
        msa_file_obj = open("{}{}".format(self.msa_path, self.msa_name), 'r')
        msa_header = msa_file_obj.readline().rstrip()
        msa_seq = msa_file_obj.readline().rstrip()
        msa_file_obj.close()
        return msa_seq

    def combine_pdb_seq(self):
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
        full_seq = first_seq + second_seq
        logging.debug("{}\t\tPDB seq length: {}".format(fname, len(full_seq)))
        logging.info("\tFinished cat-ting sequences.".format(fname))
        return full_seq

    def distance_matrix(self):
        from pandas import read_csv
        from itertools import combinations_with_replacement
        import time
        fname = "(pdb_map)"

        if ".cif" in self.pdbfile:
            # output filename
            filename = "{}distance_matrix_{}_{}A.txt".format(self.pdb_path, self.pdbfile.strip(".cif"), self.cutoff)
        else:
            # output filename
            filename = "{}distance_matrix_{}_{}A.txt".format(self.pdb_path, self.pdbfile.strip(".pdb"), self.cutoff)

        fileout = open(filename, 'w')

        # write header for output file
        fileout.write("i\tj\td\tchain_1\tchain_2\n")

        # output list of residues from pdb
        residues = self.get_residues()

        # make each possible pairs of residues
        pair_indices = combinations_with_replacement(range(len(residues)), 2)

        start_time = time.time()
        for i, j in pair_indices:
            logging.captureWarnings(True)
            res_a = residues[i]
            res_b = residues[j]
            # get chain id
            if res_a.has_id("CA") and res_b.has_id("CA"):
                chain_a = res_a.get_parent().id
                chain_b = res_b.get_parent().id
                dist = self.calc_ca_distance(res_a, res_b)
                if self.cutoff >= dist > 0.0:
                    fileout.write("%d\t%d\t%f\t%s\t%s\n" % (i + 1, j + 1, dist, chain_a, chain_b))
            else:
                print("{} NOTE! Res {} \n\tor {} not calculated! (missing CA)\n".format(fname, res_a.get_full_id(),
                                                                                        res_b.get_full_id()))

        print("{}\t -- MAIN LOOP TIME -- %s".format(fname, time.time() - start_time))
        # makes a pandas dataframe
        df_pdb = read_csv(filename, delim_whitespace=True)
        df_mon = df_pdb[df_pdb['chain_1'] == df_pdb['chain_2']]
        df_inter = df_pdb[df_pdb['chain_1'] != df_pdb['chain_2']]

        return df_pdb, df_mon, df_inter

    @staticmethod
    def calc_ca_distance(res_a, res_b):
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

    def get_residues(self):
        """
        Build a simple list of residues from a single chain of a PDB file.
        :return: A list of Bio.PDB.Residue objects.
        """
        import Bio.PDB
        fname = "(get_residues)"
        pdb_id = self.pdbfile.strip('.cif')
        print('{}\tprocessing {} file...'.format(fname, self.pdbfile))
        if pdb_id == 'pdb':
            parser = Bio.PDB.PDBParser()
        else:
            parser = Bio.PDB.MMCIFParser()

        struct = parser.get_structure(pdb_id, "{}{}".format(self.pdb_path, self.pdbfile))
        model = struct[0]
        if len(self.chain_ids) == 0:
            # get residues from every chain.
            chains = model.get_list()
        else:
            chains = [model[ch_id] for ch_id in self.chain_ids]

        residues = []
        for ch in chains:
            # make sure res are standard AA
            num_residues = 0
            for res in filter(lambda r: Bio.PDB.is_aa(r), ch.get_residues()):
                if Bio.PDB.is_aa(res, standard=True):
                    residues.append(res)
                    num_residues += 1
                else:
                    sys.stderr.write("WARNING: non-standard AA at %r%s" %
                                     (res.get_id(), os.linesep))
            self.chain_lengths.append(num_residues)

        return residues

    def plot_pdb_map(self, heatmap=None):
        """

        :param heatmap:
        :return:
        """
        import matplotlib.pylab as plt

        print("starting pdb_map...")
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
        fname = "(cm_make)"
        import numpy as np
        dca_matrix = np.loadtxt(score_matrix)
        filename = score_matrix.strip(".txt")
        logging.info("{}\t\tScore matrix filename: {}".format(fname, filename))

        x_output = []
        n = dca_matrix.shape[0]
        logging.debug("{}\tshape of matrix: {}".format(fname, n))
        for i in range(n - 1):
            for j in range(i + 1, n):
                x_output.append([i + 1, j + 1, dca_matrix[i, j]])
        self.dca_array = np.array(x_output)
        return self.dca_array

    @staticmethod
    def map_msa_to_pdb(pdbseq, msaseq):
        """
        Taken from
        https://github.com/bsir/dca-frustratometer/blob/master/dca_frustratometer.py
        :param pdbseq: PDB seq string
        :param msaseq: MSA seq string
        :return:
        """
        fname = "(map_msa_to_pdb)"
        logging.debug(fname + "\taligning msa sequences to pdb...")
        if msaseq in pdbseq:
            pdbstart = pdbseq.find(msaseq)
            fastastart = 0
        elif pdbseq in msaseq:
            fastastart = msaseq.find(pdbseq)
            pdbstart = 0

        else:
            import re
            from Bio import pairwise2
            alignments = pairwise2.align.globalxx(msaseq, pdbseq)

            fastastart = re.search("[A-Z]", alignments[0][1]).start()
            pdbstart = re.search("[A-Z]", alignments[0][0]).start()
            logging.debug("fastastart: {}\tpdbstart: {}".format(fastastart, pdbstart))

        n = min(len(msaseq), len(pdbseq))

        pdb_indices = range(pdbstart, n + pdbstart)
        dca_indices = range(fastastart, n + fastastart)
        map_to_dca = dict(zip(pdb_indices, dca_indices))
        map_to_pdb = dict(zip(dca_indices, pdb_indices))

        logging.debug(fname + "\t\tn: %s pdb indices: %s dca indices: %s" % (n, pdb_indices, dca_indices))
        logging.info(fname + "\tFinished aligning msa sequences to pdb...")
        return dca_indices, pdb_indices, map_to_dca, map_to_pdb

    def backmap(self, map_dictionary, dca_start, score_matrix):
        """

        :param map_dictionary:
        :param dca_start:
        :param score_matrix:
        :return: Dataframe
        """
        import numpy as np
        import pandas as pd
        fname = "(backmap)"
        self.dca_array = self.read_dca_matrix(score_matrix)
        basename = os.path.basename(score_matrix)
        self.result_path = "{}\\".format(os.path.dirname(score_matrix))
        if map_dictionary:
            logging.debug("{}\t\tMap dictionary given - PROCEED with mapped contact map.".format(fname))
            map_dca_array = self.apply_map(map_dictionary, dca_start)
            df_map_dca = pd.DataFrame(map_dca_array, columns=['i', 'j', 'score'])
            df_map_dca = df_map_dca.sort_values(ascending=False, by=['score'])
            logging.debug("sorted df_map_dca head {}".format(df_map_dca.head()))
            np.savetxt("{}mapped_cm_{}".format(self.result_path, basename), df_map_dca, fmt='%d\t%d\t%f')
            return df_map_dca
        else:
            logging.debug(fname + "\t\tNo map dictionary given - PROCEED with unmapped contact map.")
            df_dca = pd.DataFrame(self.dca_array, columns=['i', 'j', 'score'])
            df_dca = df_dca.sort_values(ascending=False, by=['score'])
            logging.debug("sorted df_dca head {}".format(df_dca.head()))
            np.savetxt("{}cm_{}".format(self.result_path, basename), df_dca, fmt='%d\t%d\t%f')
            return df_dca

    def apply_map(self, map_dictionary, dca_start):
        import numpy as np
        map_dca_list = []
        for i, j, score in self.dca_array:
            if i - 1 >= dca_start and j - 1 >= dca_start:
                map_index_i = map_dictionary[int(i) - 1]
                map_index_j = map_dictionary[int(j) - 1]
                map_dca_list.append([map_index_i, map_index_j, score])
        self.map_dca_array = np.array(map_dca_list)
        return self.map_dca_array

