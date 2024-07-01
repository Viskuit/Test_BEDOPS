# Testing negative data relationship
import argparse
import pandas as pd
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from subfamilies_global_functions import blastn_dic, fasta_creator, blastn_blaster, find_maximal_sets, join_conflicted_sequences, count_sequences, save_sequences_to_csv_pandas

# =============================================================================
# Defining argparser arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Make the relationship between the negative database')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to the negative CSV data SIDER file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save te data')
    return parser.parse_args()

# =============================================================================
# Defining needed functions
# =============================================================================




# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # Get the working folder
    path_working_folder = os.path.join(args.output, 'negative_noCDS_relationship')
    os.makedirs(args.output, exist_ok=True)

    # Get the negative data csv
    data = pd.read_csv(args.csv_file, sep=',', header=0)

    # Make the BLASTn dict
    path_dict_folder = os.path.join(path_working_folder, 'blast_dict')
    os.makedirs(path_dict_folder, exist_ok=True)
    path_dict_file = os.path.join(path_dict_folder, 'negative_noCDS.fasta')

    # Create the fasta file
    fasta_creator(data=data, fasta_output_path=path_dict_file)

    # Create the BLAST database
    blastn_dic(path_input=path_dict_file, path_output=path_dict_file)

    # BLASTN
    blastn_df = blastn_blaster(query=path_dict_file, path_genome=path_dict_file, identity=85)

    # Filter by length
    blastn_df = blastn_df[blastn_df['length'] > 100].copy()

    # Create a dictionary with the sequences an their respective relationships
    main_dict = {}
    for query in blastn_df['qseqid'].unique():
        values = blastn_df[blastn_df['qseqid'] == query].loc[:, ['sseqid']].values.flatten().tolist()  # For each query, get the values from 'sseqid' that match with that query.
        values = list(set(values))  # Remove duplicates
        values = sorted(values)
        main_dict[query] = values
    main_dict = {key: sorted(value, key=lambda x: int(re.findall(r'_(\d+)_', x)[0])) for key, value in main_dict.items()}

    # Remove duplicated values. The ones that are exactly the same
    unique_values = []
    [unique_values.append(value) for _, value in main_dict.items() if value not in unique_values]

    # Now remove the values that are a subset of a bigger set
    unique_values = find_maximal_sets(sequences=unique_values)

    # Now join the conflicted values. When there are two sets that share a value, join them
    unique_values = join_conflicted_sequences(sequences=unique_values)

    # Make the warning for repeated elements inside a set, if it appears more than once.
    repeated_elements = []
    counters = count_sequences(sequences=unique_values)
    if len(counters) > 0:
        repeated_elements.append("Repeated values:")
        combined_list = [item for pair in zip(counters.keys(), counters.values()) for item in pair]
        repeated_elements.append(combined_list)  # Append the sequences that appear more than once
    save_sequences_to_csv_pandas(data=unique_values, filename=os.path.join(path_working_folder, 'negative_noCDS_relationship.csv'))
    print(f'{len(repeated_elements)} repeated values: {repeated_elements}')