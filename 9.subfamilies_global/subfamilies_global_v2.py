# R. Pacheco #
# =============================================================================
# Importing needed modules
# =============================================================================
# import numpy as np
import pandas as pd
import subprocess
import os 
import sys

# from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =============================================================================
# Importing from subfamilies_global_functions.py
# =============================================================================
from subfamilies_global_functions import blastn_dic, blastn_blaster, fasta_creator, save_sequences_to_csv_pandas, find_maximal_sets, count_sequences, join_conflicted_sequences, subfamily_naming

# =============================================================================
# Main action
# =============================================================================
def grouping_subfamily(input_data, output_path, naming):
    repeated_chr = []
    for index, group in input_data:
        folder_path = os.path.join(output_path, index)  # normalizing the path

        cmd = ["mkdir", "-p", folder_path]
        subprocess.run(cmd)  # creating a folder for each chromosome
        print(f"Chromosome {index} folder created")

        fasta_path = os.path.join(folder_path, f"positives_{index}.fasta")
        fasta_creator(group, fasta_path)  # creating a for each chromosome
        blastn_dic(fasta_path)  # creating a BLASTN dict for each chromosome

        data = blastn_blaster(fasta_path, fasta_path, 85)  # BLASTN for each chromosome
        data = data[data["length"].astype(int) > 100]  # filtering by length
        # data = data.query("qseqid != sseqid")  # filtering by qseqid != sseqid

        dict = {}
        for seq in data["qseqid"].unique():
            values = data[data["qseqid"] == seq].loc[:, ["sseqid"]].values.flatten().tolist()  # get values where key = seq. And get the values in a list.
            values = list(set(values))  # remove duplicates
            values = sorted(values)  # sort values.

            dict[seq] = values
        dict = {key: sorted(value, key=lambda x: int(x.split('_')[1])) for key, value in dict.items()}   # Sort the values by the number of the sequence

        dataset = []
        for _, value in dict.items():
            if value in dataset:
                continue  # If it's in the dataset, it doesn't do anything.
            else:
                dataset.append(value)  # If it's not in the dataset, it appends it.
        
        dataset = find_maximal_sets(dataset)  # Find the maximal sets. If there's a list [1, 2] and a list [1, 2, 3], it will remove the first one.

        dataset = join_conflicted_sequences(dataset)  # Join the conflicted sequences. If there's a list [1, 2] and a list [2, 3], it will join them in a list [1, 2, 3]

        named_elements = subfamily_naming(index, naming, dataset)  # Naming the subfamilies

        counters = count_sequences(dataset)  # Count the sequences if there is an element in more than one subfamily, i.e.,  >1
        if len(counters) > 0:
            repeated_chr.append(index)
            dataset.append("Repeated values:")
            combined_list = [item for pair in zip(counters.keys(), counters.values()) for item in pair]
            dataset.append(combined_list)  # Append the sequences that appear more than once

        save_sequences_to_csv_pandas(dataset, os.path.join(folder_path, f"families_{index}.csv"))  # saving the dataset
        save_sequences_to_csv_pandas(named_elements, os.path.join(folder_path, f"NAMED_families_{index}.csv"))  # saving the dataset
    
    print(f"{len(repeated_chr)} chromosomes with repeated values: {repeated_chr}")


if __name__ == "__main__":
        # =============================================================================
    # Importing data
    # =============================================================================
    input_path = os.path.normpath(sys.argv[1])  # Import data
    data = pd.read_csv(input_path, sep="\,", header=None)  # Import data
    data_grouped = data.groupby(0)  # Group data by the first column, i.e., chromosomes.

    path_out = os.path.normpath(sys.argv[2])  # Path to save the files. Finish with "/", e.g. "./dict/chromosomes/"

    grouping_subfamily(data_grouped, path_out, sys.argv[3])  # sys.argv[3] is the naming of the elements like "sre_c10.20A"