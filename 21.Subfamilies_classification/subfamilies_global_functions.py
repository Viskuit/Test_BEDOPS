# R. Pacheco #
# =============================================================================
# Importing needed modules
# =============================================================================
import pandas as pd
import subprocess
import os
import re
import string  # For the subfamily_naming function

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =============================================================================
# Defining needed functions
# =============================================================================
def blastn_dic(path_input, path_output):
    cmd = f'makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}'
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# -----------------------------------------------------------------------------
def blastn_blaster(query, path_genome, identity):
    cmd = (
        f'blastn -word_size 11 '
        f'-query {query} '
        f'-db {path_genome} '
        f'-perc_identity {identity} '
        f'-outfmt 10'
    )
    data = subprocess.run(cmd, shell=True, capture_output=True, text=True, universal_newlines=True, executable='/usr/bin/bash')
    data = data.stdout
    data_df = pd.DataFrame(
        [x.split(',') for x in data.split('\n') if x],
        columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'] 
    )
    if not data_df.empty:
        data_df[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']] = data_df[['pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']].apply(pd.to_numeric)  # Convert to numeric
    else:  # If empty, return an empty dataframe
        return pd.DataFrame()

    return data_df  # I only need to return the data, no need to transform it into a pandas dataframe since I want to check the number of lines with the line separator '\n'


# -----------------------------------------------------------------------------
# Get fasta f# Define fasta_creator file from the first chromosome
def fasta_creator(data, fasta_output_path):
    matrix = []
    for index, sequence in data.iterrows():
        rec = SeqRecord(Seq(sequence['sseq']),
                        id=f"Seq_{index}_{sequence['sseqid']}",
                        description="Leishmania infantum"
                        )
        matrix.append(rec)
    SeqIO.write(matrix, fasta_output_path, "fasta")

# -----------------------------------------------------------------------------
def save_sequences_to_csv_pandas(data, filename):
    # Convert the list of lists to a DataFrame
    df = pd.DataFrame(data)
    # Save the DataFrame to a CSV file
    df.to_csv(filename, header=False, index=False)

# -----------------------------------------------------------------------------
def find_maximal_sets(sequences):
    # Identify maximal sets that are not subsets of any other set in the list
    # If there are repeated lists. And those lists are the maximun sets, it will break.
    maximal_sets = []
    for i, seq in enumerate(sequences):
        is_subset = False
        for j, other_seq in enumerate(sequences):
            # print(f"For index {i} and {j}, sequences are {seq} and {other_seq}")
            if i != j and set(seq).issubset(set(other_seq)):
                # print("\tChecking if it's a subset")
                is_subset = True
                break
        if not is_subset:
            # print(f"\tAdding {seq} to maximal_sets")
            maximal_sets.append(seq)
    return maximal_sets

# -----------------------------------------------------------------------------
def join_conflicted_sequences(sequences):
    filtered_data = []  # Save the data after all the filtering in this code
    joines_rel = {}  # Save the index relations that we have already joined.
    for i, seq1 in enumerate(sequences):
        counter_no_coincidence = 0  # When seq1 loops through every seq2 and there is no coincidence, we will increase the counter by 1. If the counter is equal to the length of the sequences, we will append the seq1 to the filtered_data, since there is no coincidence with any other seq2.

        for j, seq2 in enumerate(sequences):
            if i == j:  # Skip the same list
                counter_no_coincidence += 1
                continue
            
            # Number of coincidence of each element in seq1 with seq2
            count_coincidences = sum(1 for x in seq1 if x in seq2)

            if count_coincidences >= 1:
                # Since the two seqs have relations, we will join them
                new_seq = list(set(seq1 + seq2))  # Join the two lists and remove the repeated elements
                new_seq = sorted(new_seq)  # Sort the new list
                filtered_data.append(new_seq)  # Append the new list to the filtered_data
                
                # Foolproof for the next iterations
                seq1 = new_seq  # Update the values for the comparison with the next seq2 in the loop
                sequences[j] = new_seq  # Update the value for the next iteration when seq2 is called again the the "j" index. It has to be changed in the original sequences list.

                # Save the relations
                ## Save which elements did "j/seq2" matched with "i/seq1"
                if j in joines_rel.keys():
                    joines_rel[j].append(i)
                else:
                    joines_rel[j] = [i]
            else:  # If there is no coincidence, we will increase the counter by 1
                counter_no_coincidence += 1
        if counter_no_coincidence == len(sequences):
            filtered_data.append(seq1)
    
    # Remove duplicates lists inside the lists of lists
    filtered_data_NoDup = []
    [filtered_data_NoDup.append(value) for value in filtered_data if value not in filtered_data_NoDup]

    data = find_maximal_sets(filtered_data_NoDup)
    return data

# -----------------------------------------------------------------------------
def count_sequences(sequences):
    from collections import Counter
    # Flatten the list of sequences
    flat_list = [seq for sublist in sequences for seq in sublist]
    # Count each sequence using Counter
    sequence_counts = Counter(flat_list)
    # Identify sequences that appear more than once
    multiple_occurrences = {seq: count for seq, count in sequence_counts.items() if count > 1}
    return multiple_occurrences

# -----------------------------------------------------------------------------
def numbering_dict(list_array):
    new_list = [name for element in list_array for name in element]
    new_list = sorted(new_list, key=lambda x: int(re.search(r'_(\d+)_', x).group(1)))
    counter = 1
    dict_naming = {}
    for name in new_list:
        dict_naming[name] = counter
        counter += 1
    return dict_naming

def generate_sequence():
    letters = string.ascii_uppercase
    single_letter_sequences = list(letters) # Generate single letter sequences
    double_letter_sequences = [a + b for a in letters for b in letters] # Generate double letter sequences
    triple_letter_sequences = [a + b + c for a in letters for b in letters for c in letters] # Generate triple letter sequences
    all_sequences = single_letter_sequences + double_letter_sequences + triple_letter_sequences # Combine all sequences
    return all_sequences

def subfamily_naming(chromosome, naming, sequences, naming_dict):
    abc_seq = generate_sequence()  # it will generate from A to ZZZ
    letter_slider = 0
    named_elements = []
    # subfamily_abc_counter_id = "A"
    for _, seqs in enumerate(sequences):  # selecting list
        if len(seqs) > 1:
            subfamily = []
            for _, element in enumerate(seqs):  # selecting element inside a list.
                chr_id = re.search(r"\d+", chromosome).group(0)
                number_id = naming_dict[element]  # get the correct number id
                number_id = str(number_id)[::-1].zfill(len(str(number_id)) + 1)[::-1]  # fill the number with 0
                subfamily.append(f"{naming}_c{chr_id}.{number_id}{abc_seq[letter_slider]}")
            letter_slider += 1
            named_elements.append(subfamily)
        else: # if len(seqs) == 1
            chr_id = re.search(r"\d+", chromosome).group(0)
            number_id = naming_dict[seqs[0]]  # get the correct number id
            number_id = str(number_id)[::-1].zfill(len(str(number_id)) + 1)[::-1]  # fill the number with 0
            named_elements.append([f"{naming}_c{chr_id}.{number_id}"])

    return named_elements
