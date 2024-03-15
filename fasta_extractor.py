import argparse
from Bio import SeqIO

def fasta_extractor(pathfile, outfile, extract_list):
    """
    This function let me chose the sequences I want from the fasta by their number.

    Packages needed:
       - ``from Bio import SeqIO``
       - ``from pathlib import Path``

    :param pathfile: Path to the fasta file we want to subset sequences
    :type pathfile: string

    :param outfile: Path to the output file. Include name and ".fasta" suffix.

    :param extract_list: list with our index numbers to extract the sequence. For example [2, 4, 6] will extract sequences 2, 4 and 6 from ``pathfile``
    :type extract_list: list

    :return: A fasta file with the sequence we want
    :rtype: fasta file
    """
    with open(outfile, "w") as out_file:
        # Remember "enumerate" starts in "1"
        for count, fasta in enumerate(SeqIO.parse(open(pathfile), "fasta"), 1):  # from Bio import SeqIO
            # name, sequence = fasta.id, str(fasta.seq)
            if count in extract_list:
                SeqIO.write(fasta, out_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract sequences from a fasta file.')
    parser.add_argument('pathfile', type=str, help='Path to the fasta file we want to subset sequences')
    parser.add_argument('outfile', type=str, help='Path to the output file. Include name and ".fasta" suffix.')
    parser.add_argument('extract_list', type=lambda s: list(map(int, s.split(','))), help='List with our index numbers to extract the sequence. For example "2,4,6" will extract sequences 2, 4 and 6 from pathfile')

    args = parser.parse_args()

    fasta_extractor(args.pathfile, args.outfile, args.extract_list)
