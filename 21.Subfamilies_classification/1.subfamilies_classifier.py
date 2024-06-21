# Needed modules
import argparse
import pandas as pd
import subprocess
import os
import re

from Bio import Seq
from Bio.SeqRecord import SeqRecord

from subfamilies_global_functions import blastn_dic, blastn_blaster, fasta_creator, save_sequences_to_csv_pandas, find_maximal_sets, count_sequences, join_conflicted_sequences, subfamily_naming, numbering_dict

# =============================================================================
# Defining arparse arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Classify SIDER elements in subfamilies.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the CSV file with the positive SIDER elements.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data.')
    parser.add_argument('-n', '--naming', type=str, required=True, help='Naming of the subfamilies.')
    return parser.parse_args()

# =============================================================================
# Defining needed functions
# =============================================================================
def grouping_subfamily(input_df, path_output, naming):
    # Group <input_df> by chromosome
    repeated_chr = []  # List to store the repeated chromosomes
    input_df_grouped = input_df.groupby('sseqid')  # Data type: pandas groupby object
    for index, group in input_df_grouped:  # group is a pandas dataframe
        # Set working folder
        path_chr_folder = os.path.join(path_output, index)   # Folder for each chr, i.e., index
        os.makedirs(path_chr_folder, exist_ok=True)  # Create the folder
        print(f'Chromosome {index} folder created at {path_chr_folder}')

        # Create the fasta file
        path_fasta_file = os.path.join(path_chr_folder, f'{index}.fasta')
        fasta_creator(data=group, fasta_output_path=path_fasta_file)

        # Create the BLAST database
        path_dict_folder = os.path.join(path_chr_folder, 'blast_dict')
        os.makedirs(path_dict_folder, exist_ok=True)
        path_dict_file = os.path.join(path_dict_folder, os.path.basename(path_fasta_file))
        blastn_dic(path_input=path_fasta_file, path_output=path_dict_file)

        # BLASTN
        data = blastn_blaster(query=path_fasta_file, path_genome=path_dict_file, identity=85)  # in each chr launch the BLASTn sequence against each other.

        # Skip if the dataframe is empty
        if data.empty:
            continue  # If the dataframe is empty, continue with the next iteration
        else:
            pass

        # Filter by length
        data = data[data['length'] > 100].copy()  # Filter by length

        # Create a dictionary with the sequences
        main_dict = {}
        for query in data['qseqid'].unique():
            values = data[data['qseqid'] == query].loc[:, ['sseqid']].values.flatten().tolist()  # For each query, get the values from 'sseqid' that match with that query.
            values = list(set(values))  # Remove duplicates
            values = sorted(values)
            main_dict[query] = values
        main_dict = {key: sorted(value, key=lambda x: int(re.findall(r'_(\d+)_', x)[0])) for key, value in main_dict.items()}  # Inside each key in the dicctionary, sort the values by the number of the sequence

        unique_values = []
        [unique_values.append(value) for _, value in main_dict.items() if value not in unique_values]
        # Get the unique values from the dictionary
        unique_values = find_maximal_sets(sequences=unique_values)  # Find the maximal sets. # Find the maximal sets. If there's a list [1, 2] and a list [1, 2, 3], it will remove the first one.
        unique_values = join_conflicted_sequences(sequences=unique_values)  # Join the conflicted sequences. # Join the conflicted sequences. If there's a list [1, 2] and a list [2, 3], it will join them in a list [1, 2, 3]
        numering_dict = numbering_dict(list_array=unique_values)
        named_elements = subfamily_naming(chromosome=index, naming=naming, sequences=unique_values, naming_dict=numering_dict)
        counters = count_sequences(sequences=unique_values)

        if len(counters) > 0:
            repeated_chr.append(index)
            unique_values.append("Repeated values:")
            combined_list = [item for pair in zip(counters.keys(), counters.values()) for item in pair]
            unique_values.append(combined_list)  # Append the sequences that appear more than once

        save_sequences_to_csv_pandas(data=unique_values, filename=os.path.join(path_chr_folder, f'families_{index}.csv'))
        save_sequences_to_csv_pandas(data=named_elements, filename=os.path.join(path_chr_folder, f'NAMED_families_{index}.csv'))

    print(f'{len(repeated_chr)} chromosomes with repeated values: {repeated_chr}')

        
# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()  # Parse arguments

    # Read data needed
    data = pd.read_csv(args.file, sep=',', header=0)

    # Create the output folder
    path_folder = os.path.join(args.output, 'subfamilies_classification')
    os.makedirs(path_folder, exist_ok=True)
    print(f'Principal folder created in {path_folder}')

    # call the function
    grouping_subfamily(input_df=data, path_output=path_folder, naming=args.naming)


