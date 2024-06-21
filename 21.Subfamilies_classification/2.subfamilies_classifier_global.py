import argparse
import pandas as pd
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from subfamilies_global_functions import blastn_dic, blastn_blaster, find_maximal_sets, join_conflicted_sequences, numbering_dict, subfamily_naming, count_sequences, save_sequences_to_csv_pandas

# =============================================================================
# Defining argparser arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Classify SIDER elements in subfamilies but global.')
    parser.add_argument('-p', '--main_path', type=str, required=True, help='Path to the main folder.')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to the CSV file with the positive SIDER elements.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data.')
    return parser.parse_args()

# =============================================================================
# Defining needed functions
# =============================================================================
def fasta_creator_dict(main_dict, fasta_output_path):
    matrix = []
    for key, value in main_dict.items():
        rec = SeqRecord(Seq(value[1]),
                        id=f"{key}|{value[0]}",  # Key is the ID like "Seq_1_LinJ.01" amd value[0] is the last subfamily identifier.
                        description="Leishmania infantum"
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")


# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # get the working folder
    os.makedirs(args.output, exist_ok=True)

    #Get the folders name
    folders = [folder for folder in os.listdir(args.main_path) if os.path.isdir(os.path.join(args.main_path, folder))]
    folders = sorted(folders, key=lambda x: int(x.split('.')[1]))

    # Get the files with the ID form the main data base and the NAMED data base
    files_dirs = []
    files_dirs_named = []
    for folder in folders:
        path_chr_folder = os.path.join(args.main_path, folder)
        file = [os.path.join(path_chr_folder, file) for file in os.listdir(path_chr_folder) if file.endswith(".csv") and file.startswith("families")]  # There will be only one file in the folder
        file_named = [os.path.join(path_chr_folder, file) for file in os.listdir(path_chr_folder) if file.endswith(".csv") and file.startswith("NAMED_families")]
        files_dirs.append(file[0])
        files_dirs_named.append(file_named[0])

    # Get the main data frame
    main_df = pd.read_csv(args.csv_file, sep=',', header=0)

    # Lets get the dictionary with the sequences and IDs
    main_dict = {}
    for file, file_named in zip(files_dirs, files_dirs_named):
        file_df = pd.read_csv(file, sep=',', header=None)
        file_named_df = pd.read_csv(file_named, sep=',', header=None)
        for (index, row), (index_named, row_named) in zip(file_df.iterrows(), file_named_df.iterrows()):
            row = [value for value in row if str(value) != 'nan']
            row_named = [value for value in row_named if str(value) != 'nan']
            # Start the process
            if len(row) == 1:  # For the orphans
                slice_id = int(re.search(r'_(\d+)_', row[0]).group(1))
                main_dict[row[0]] = [row_named[0], main_df.loc[slice_id]["sseq"]]  # Get the sequence
            else:  # For the families
                subfamilies_seqs = {}
                for value in row:
                    slice_id = int(re.search(r'_(\d+)_', value).group(1))
                    subfamilies_seqs[value] = [row_named[0], main_df.loc[slice_id]["sseq"]]
                # get the index dict in subfamilies_seqs that has the longer seq
                max_len = 0
                max_index = ""
                for key, value in subfamilies_seqs.items():
                    if len(value[1]) > max_len:
                        max_len = len(value[1])
                        max_index = key
                main_dict[max_index] = subfamilies_seqs[max_index]

    # Now let's create a fasta file with the sequences
    path_fasta_file = os.path.join(args.output, 'sequences.fasta')
    fasta_creator_dict(main_dict=main_dict, fasta_output_path=path_fasta_file)

    # Now create a dictionary
    path_dict_folder = os.path.join(args.output, 'blastn_dict')
    os.makedirs(path_dict_folder, exist_ok=True)
    path_dict_file = os.path.join(path_dict_folder, os.path.basename(path_fasta_file))

    # Create the BLAST database
    blastn_dic(path_input=path_fasta_file, path_output=path_dict_file)

    # BLASTN the sequences
    data = blastn_blaster(query=path_fasta_file, path_genome=path_dict_file, identity=85)

    # Skip if the dataframe is empty
    if data.empty:
        print('Blaster returned an empty dataframe.')
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
    unique_values = find_maximal_sets(sequences=unique_values)  # Find the maximal sets. # Find the maximal sets. If there's a list [1, 2] and a list [1, 2, 3], it will remove the first one.
    unique_values = join_conflicted_sequences(sequences=unique_values)  # Join the conflicted sequences. # Join the conflicted sequences. If there's a list [1, 2] and a list [2, 3], it will join them in a list [1, 2, 3]

    repeated_elemns = []
    counters = count_sequences(sequences=unique_values)
    if len(counters) > 0:
        repeated_elemns.append("Repeated values:")
        combined_list = [item for pair in zip(counters.keys(), counters.values()) for item in pair]
        repeated_elemns.append(combined_list)

    save_sequences_to_csv_pandas(data=unique_values, filename=os.path.join(args.output, 'families.csv'))
    print(f'{len(repeated_elemns)} repeated values: {repeated_elemns}')


    


