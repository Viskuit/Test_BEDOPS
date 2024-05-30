# R. Pacheco #
# =============================================================================
# Importing needed modules
# =============================================================================
import pandas as pd
import subprocess
import os
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# =============================================================================
# Defining needed functions
# =============================================================================
def blastn_dic(path_input):
    os.system("makeblastdb -in " + path_input + " -dbtype nucl -parse_seqids")
    print("\nBlast Dictionary created in", path_input)

# -----------------------------------------------------------------------------
def blastn_blaster(query, db, perc_indentity):
    cmd = "blastn -word_size 11 " \
    + " -query " + query \
    + " -db "  + db \
    + " -perc_identity " + str(perc_indentity) \
    + " -outfmt '10 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp sstrand sseq'"
    data = subprocess.check_output(cmd, shell=True, universal_newlines=True)
    data = pd.DataFrame([x.split(",") for x in data.split("\n") if x])

    data.columns = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "sstrand", "sseq"]

    return data

# -----------------------------------------------------------------------------
# Get fasta f# Define fasta_creator file from the first chromosome
def fasta_creator(data, fasta_output_path):
    matrix = []
    for index, sequence in data.iterrows():
        rec = SeqRecord(Seq(sequence[5]),  # In the 5 position is the seq
                        id=f"Seq_{index}_{sequence[0]}",
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
    new_list = []
    skip_index_dict = {}  # Save the indexes that we have already joined.
    for index, seqs in enumerate(sequences):
        counter_no_coincidence = 0  # For when there is no coincidence between the two lists.

        for index2, seqs2 in enumerate(sequences):
            try:
                if index2 in skip_index_dict[index]:  # If the index2 is in the skip_index_dict, we skip it.
                    continue
            except:
                pass
            
            if index == index2:  # with this we avoid comparing the same list
                counter_no_coincidence += 1
                continue

            seqs_pdS = pd.Series(seqs)  # transform the first list into a pandas series, so we can use .apply() functions.

            try: # check if value True exist in the data
                count_coundicence = seqs_pdS.apply(lambda x: x in seqs2).value_counts().loc[True]  # Count the number of coincidences between the two lists.
            except:
                count_coundicence = 0  # If it doesn't exist, we will set it to 0.

            # For when both lists have coincidences, we will join them.
            if count_coundicence >= 1:  # If there is at least one coincidence, we will join the two lists.
                new_seqs = list(set(seqs + seqs2))  # Join the two lists and remove the repeated elements.
                new_seqs = sorted(new_seqs)  # Sort the new list.
                new_list.append(new_seqs)  # Append the new list to the new_list
                seqs = new_seqs  # Update the index to the new list.

                sequences[index2] = new_seqs  # Update the list in the sequences list. So when it's its turn, it will be updated in a more fiable way.

                if index2 in skip_index_dict.keys():  # If the index is in the skip_index_dict, we will append the index2 to the list.
                    skip_index_dict[index2].append(index)  # for index2, we save that it has been joined with index, so when its turn for index2 to be index in the loop, it skips the known index.
                else:
                    skip_index_dict[index2] = [index]
            

            if count_coundicence == 0:
                counter_no_coincidence += 1

        if counter_no_coincidence == len(sequences):
            new_list.append(seqs)
    
    # Remove duplicatse lists inside the lists of lists.
    new_list2 = []
    for _, value in enumerate(new_list):
        if value in new_list2:
            continue
        else:
            new_list2.append(value)

    # data1 = find_maximal_sets(sequences)
    data2 = find_maximal_sets(new_list2)

    return data2


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
def subfamily_naming(chromosome, naming, sequences):
    named_elements = []
    orphan_counter_id = 10
    subfamily_abc_counter_id = "A"
    for index, seqs in enumerate(sequences):  # selecting list
        if len(seqs) > 1:
            subfamily = []
            subfamily_counter_member_id = 10
            for index2, element in enumerate(seqs):  # selecting element inside a list.
                chr_id = re.search(r"\d+", chromosome).group(0)
                subfamily.append(f"{naming}_c{chr_id}.{subfamily_counter_member_id}{subfamily_abc_counter_id}")
                subfamily_counter_member_id += 10
            subfamily_abc_counter_id = chr(ord(subfamily_abc_counter_id) + 1)


            named_elements.append(subfamily)
        else: # if len(seqs) == 1
            chr_id = re.search(r"\d+", chromosome).group(0)
            named_elements.append([f"{naming}_c{chr_id}.{orphan_counter_id}"])
            orphan_counter_id += 10

    return named_elements

