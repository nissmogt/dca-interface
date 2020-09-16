def straighten_msa(msa_file):
    from Bio import SeqIO
    print(msa_file)
    records = list(SeqIO.parse(msa_file, 'fasta'))
    with open("{}s".format(msa_file), "w") as output:
        for i in range(1, len(records)):
            header_a = records[i].id.split("|")[1]
            header_b = records[i].id.split("|")[3]
            # header = records[i].id.split("/")[0]
            # output.write(">{}\n{}\n".format(header, records[i].seq))
            output.write(">{}_{}\n{}\n".format(header_a, header_b, records[i].seq))
    output.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Rewrites MSA file so that sequence is in one row. '
                                                 'NOTE: Removes first sequence.')
    parser.add_argument('msa', help='msa file')
    args = parser.parse_args()
    msa = args.msa
    straighten_msa(msa)


if __name__ == '__main__':
    main()
