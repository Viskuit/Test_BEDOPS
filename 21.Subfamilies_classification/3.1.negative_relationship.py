# Testing negative data relationship
import argparse
import pandas as pd
import os
import re
import string

from subfamilies_global_functions import blastn_dic, fasta_creator, blastn_blaster2, find_maximal_sets, \
    join_conflicted_sequences, count_sequences, save_sequences_to_csv_pandas, numbering_dict


# =============================================================================
# Defining "argparse" arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Make the relationship between the negative database that does not have CDS match')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to the negative CSV data SIDER file that does not have CDS match')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save te data')
    return parser.parse_args()


# =============================================================================
# Defining needed functions
# =============================================================================
def generate_sequence():
    letters = string.ascii_uppercase
    single_letter_sequences = list(letters) # Generate single letter sequences
    double_letter_sequences = [a + b for a in letters for b in letters] # Generate double letter sequences
    triple_letter_sequences = [a + b + c for a in letters for b in letters for c in letters] # Generate triple letter sequences
    all_sequences = single_letter_sequences + double_letter_sequences + triple_letter_sequences # Combine all sequences
    return all_sequences

def subfamily_naming_mod(naming, sequences, naming_dict):
    abc_seq = generate_sequence()  # it will generate from A to ZZZ
    letter_slider = 0
    named_elements = []
    # subfamily_abc_counter_id = "A"
    for _, seqs in enumerate(sequences):  # selecting list
        if len(seqs) > 1:
            subfamily = []
            for _, element in enumerate(seqs):  # selecting element inside a list.
                chr_id = element.split('_')[2].split(".")[1] # get the chromosome id
                number_id = naming_dict[element]  # get the correct number id
                number_id = str(number_id)[::-1].zfill(len(str(number_id)) + 1)[::-1]  # fill the number with 0
                subfamily.append(f"{naming}_c{chr_id}.{number_id}{abc_seq[letter_slider]}")
            letter_slider += 1
            named_elements.append(subfamily)
        else: # if len(seqs) == 1
            chr_id = chr_id = seqs[0].split('_')[2].split(".")[1] # get the chromosome id
            number_id = naming_dict[seqs[0]]  # get the correct number id
            number_id = str(number_id)[::-1].zfill(len(str(number_id)) + 1)[::-1]  # fill the number with 0
            named_elements.append([f"{naming}_c{chr_id}.{number_id}"])

    return named_elements

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
    blastn_df = blastn_blaster2(query=path_dict_file, path_genome=path_dict_file, identity=85)

    # Filter by length
    blastn_df = blastn_df[blastn_df['length'] > 100].copy()
    
    # Save data for future coordinates correction
    blastn_df.to_csv(os.path.join(path_working_folder, 'blastn_neg_noMatchCds_plus100nt.csv'), index=False, header=True, sep=',')

    # Create a dictionary with the sequences In their respective relationships
    main_dict = {}
    for subject in sorted(list(blastn_df['sseqid'].unique()), key=lambda x: int(x.split('_')[1])):
        values = blastn_df.loc[
            blastn_df['sseqid'] == subject,
            ['qseqid']
            ].values.flatten().tolist()  # For each subject, get the values from 'qseqid' that match with that subject.
            
        values = list(set(values))  # Remove duplicates
        values = sorted(values)
        main_dict[subject] = values
    main_dict = {key: sorted(value, key=lambda x: int(re.findall(r'_(\d+)_', x)[0])) for key, value in
                 main_dict.items()}

    # Remove duplicated values. The ones that are exactly the same
    unique_values = []
    [unique_values.append(value) for _, value in main_dict.items() if value not in unique_values]

    # Now remove the values that are a subset of a bigger set
    unique_values = find_maximal_sets(sequences=unique_values)

    # Now join the conflicted values. When there are two sets that share a value, join them
    unique_values = join_conflicted_sequences(sequences=unique_values)
    
    numering_dict = numbering_dict(list_array=unique_values)
    named_elements = subfamily_naming_mod(naming="rejected_noCDS", sequences=unique_values, naming_dict=numering_dict)
    merge_naming =[[f"{x}|{y}" for x, y in zip(sublist1, sublist2)] for sublist1, sublist2 in zip(unique_values, named_elements)]

    # Make the warning for repeated elements inside a set, if it appears more than once.
    repeated_elements = []
    counters = count_sequences(sequences=unique_values)
    if len(counters) > 0:
        repeated_elements.append("Repeated values:")
        combined_list = [item for pair in zip(counters.keys(), counters.values()) for item in pair]
        repeated_elements.append(combined_list)  # Append the sequences that appear more than once
    save_sequences_to_csv_pandas(data=unique_values,
                                 filename=os.path.join(path_working_folder, 'negative_noCDS_relationship.csv'))
    save_sequences_to_csv_pandas(data=merge_naming,
                                 filename=os.path.join(path_working_folder, 'negative_noCDS_relationship_named.csv')
    )
    print(f'{len(repeated_elements)} repeated values: {repeated_elements}')
