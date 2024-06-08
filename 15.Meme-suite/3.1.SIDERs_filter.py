import argparse
import numpy as np
import pandas as pd
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter elements based on chromosome count and e-value.")
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('-d', '--dict_path', type=str, required=True, help='Path to the genome fasta file.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save the output CSV files.')
    return parser.parse_args()

def blastn_blaster(query_path, dict_path, evalue):
    cmd = "blastn -word_size 11 -query " \
        + query_path + " -db " \
        + dict_path \
        + " -evalue " + str(evalue) \
        + " -outfmt 10"
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    return data

def fasta_creator(sequence, index, fasta_output_path):
    rec = SeqRecord(Seq(sequence),
                    id="Seq_" + str(index),
                    description="Leishmania infantum"
                    )
    SeqIO.write(rec, fasta_output_path, "fasta")

if __name__ == "__main__":
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)

    data = pd.read_csv(args.file, sep=",", header=0)
    matches = pd.Series([False] * data.shape[0])
    not_matches = pd.Series([False] * data.shape[0])
    accepted = 0
    rejected = 0

    for index, row in data.iterrows():
        print("="*50)
        print(f"Analyzing row {index + 1} of {data.shape[0]}")

        fasta_path = os.path.join(args.output_dir, "mySequence.fasta")
        print(row)
        fasta_creator(row["sseq"], index, fasta_path)
        blastn_data = blastn_blaster(fasta_path, args.dict_path, 1.0E-09)

        if blastn_data.count("\n") <= 1:
            not_matches[index] = True
            rejected += 1
            print("\t\tREJECTED")
        else:
            blastn_data = blastn_data.strip().split("\n")
            blast_data_df = pd.DataFrame([x.split(",") for x in blastn_data if x])
            if blast_data_df[1].nunique() >= 5:
                matches[index] = True
                accepted += 1
                print("\t\tACCEPTED")
            else:
                not_matches[index] = True
                rejected += 1
                print("\t\tREJECTED")
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")

    print(f"The total number of matches is: {matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / data.shape[0] * 100, 2)}%")
    print("~"*50)
    print(f"The total number of not matches is: {not_matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of not matches is: {round(not_matches.sum() / data.shape[0] * 100, 2)}%")

    yes_data = data[matches]
    yes_data.to_csv(os.path.join(args.output_dir, "positives_testing_elements.csv"), index=False, header=True)
    no_data = data[~matches]
    no_data.to_csv(os.path.join(args.output_dir, "negatives_testing_elements.csv"), index=False, header=True)
