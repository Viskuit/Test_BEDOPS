import pandas as pd
import numpy as np
import subprocess
import os


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# print(os.getcwd())
os.chdir("/home/rfpacheco/Desktop/Projects/Testing_Leishmania_project/Filter_Test_1")

# =============================================================================
# Import data
# =============================================================================
path = "../Data/To_analyze/Results_software/SIDER2_chr32/RUNS/run_20.csv"  # Path to the file
data = pd.read_csv(path, sep=",", header=None)  # Read the file
data = data.dropna(axis=1, how='all')  # Drop columns with all NaN values only if the whole column is NaN
data.columns = range(data.columns.size)  # Reset the columns index

# =============================================================================
# Prepare functions
# =============================================================================
dict_path = "../Data/genome/TriTrypDB-67_LinfantumJPCM5_Genome.fasta"
query_path = "."

# Main BLASTn tool
def blastn_blaster(query_path, dict_path, perc_identity, evalue):
    data = subprocess.check_output("blastn -word_size 11 -query " 
                + query_path + " -db " 
                + dict_path 
                + " -perc_identity " + str(perc_identity)
                + " -evalue " + str(evalue) 
                + " -outfmt 10", shell=True, universal_newlines=True)  # Important the E value
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
# matches_start = pd.Series([False] * group_data_1.shape[0])
    
if __name__ == "__main__":

    matches = pd.Series([False] * data.shape[0])  # Create a series with False values
    accepted = 0
    rejected = 0
    for index, row in data.iterrows():
        
        print("="*50)
        print(f"Analyzing row {index + 1} of {data.shape[0]}")

        fasta_creator(row.iloc[5], index, "./mySequence_2.fasta")  # Create the tmp fasta file
        blastn_data = blastn_blaster("./mySequence_2.fasta", dict_path, 0, 1.0E-20)  # Run the blastn tool
        
        if blastn_data.count("\n") > 2:  # If there is a match
            print("\t\tACCEPTED")
            accepted += 1
            matches[index] = True
        else: 
            print("\t\tREJECTED")
            rejected += 1
        print(f"\t\t\t\t\tAccepted: {accepted} - Rejected: {rejected}")

    print(f"The toal number of mathches is: {matches.sum()} out of {data.shape[0]}")
    print(f"The percentage of matches is: {round(matches.sum() / data.shape[0] * 100, 2)}%")
    data = data[matches]  # Filter the data
    data.to_csv("./filtered_data_2.csv", index=False, header=False)  # Save the data

