import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Define fasta_creator function
def fasta_creator(data, fasta_output_path):
    matrix = []
    for index, sequence in data.iterrows():
        rec = SeqRecord(Seq(sequence['sseq']),
                        id=f"Seq_{index}_{sequence['sseqid']}",
                        description="Leishmania infantum"
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")

# Main function to parse arguments and call fasta_creator
def main():
    parser = argparse.ArgumentParser(description="Create a FASTA file from a CSV file containing sequence data.")
    
    # Define arguments
    parser.add_argument('input_csv', type=str, help="Path to the input CSV file")
    parser.add_argument('output_fasta', type=str, help="Path to the output FASTA file")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Read the input CSV file
    data = pd.read_csv(args.input_csv)
    
    # Call the fasta_creator function
    fasta_creator(data, args.output_fasta)

if __name__ == "__main__":
    main()
