# Script to filter the elements by:
## They should be in >= 5 chromosomes
## The evalue should be <= 1.0E-09

# Import needed modules
import numpy as np
import pandas as pd
import subprocess
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Let's setup the working directory, just in case
os.chdir("/home/rfpacheco/Desktop/Projects/Testing_Leishmania_project/BEDOPS_join_strands")

# =============================================================================
# Import data
# =============================================================================
file = "./sequences.csv"
dict_path = "../Data/genome/TriTrypDB-67_LinfantumJPCM5_Genome.fasta"  # Path to the genome
data = pd.read_csv(file, sep=",", header=None)  # Read the file

# =============================================================================
# Prepare functions
# =============================================================================
def blastn_blaster(query_path, dict_path, evalue):
    cmd = "blastn -word_size 11 -query " \
        + query_path + " -db " \
        + dict_path \
        + " -evalue " + str(evalue) \
        + " -outfmt 10"
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)  # Important the E value
    return data

def fasta_creator(sequence, index, fasta_output_path):
    rec = SeqRecord(Seq(sequence),
                    id="Seq_" + str(index),
                    description="Leishmania infantum"
                    )
    SeqIO.write(rec, fasta_output_path, "fasta")

# =============================================================================
# Main loop
# =============================================================================

if __name__ == "__main__":
    matches = pd.Series([False] * data.shape[0])  # Create a series with False values
    not_matches = pd.Series([False] * data.shape[0])  # Create a series with False values
    accepted = 0
    rejected = 0

    for index, row in data.iterrows():
        print("="*50)
        print(f"Analyzing row {index + 1} of {data.shape[0]}")

        fasta_creator(row.iloc[5], index, "./filtered_data/mySequence.fasta")  # Create the tmp fasta file
        blastn_data = blastn_blaster("./filtered_data/mySequence.fasta", dict_path, 1.0E-09)  # Run the blastn tool and get the data

        if blastn_data.count("\n") <= 1:
            not_matches[index] = True
            rejected += 1
            print("\t\tREJECTED")
        else:
            blastn_data = blastn_data.strip().split("\n")
            blast_data_df = pd.DataFrame([x.split(",") for x in blastn_data if x])  # Create a DataFrame from the data
            if blast_data_df[1].nunique() >= 5:
                matches[index] = True
                accepted += 1
                print("\t\tACCEPTED")
            else:
                not_matches[index] = True
                rejected += 1
                print("\t\tREJECTED")
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")  # For a counting method
    
    print(f"The total number of mathches is: {matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / data.shape[0] * 100, 2)}%")
    print("~"*50)
    print(f"The total number of not mathches is: {not_matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of not matches is: {round(not_matches.sum() / data.shape[0] * 100, 2)}%")

    
    yes_data = data[matches]  # Filter the data
    yes_data.to_csv("./filtered_data/positives_testing_elements.csv", index=False, header=False)  # Save the data
    no_data = data[~matches]  # Filter the data
    no_data.to_csv("./filtered_data/negatives_testing_elements.csv", index=False, header=False)  # Save the data     
