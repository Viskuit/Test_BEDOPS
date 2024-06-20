import pandas as pd

def join_conflicted_sequences(sequences):
    filtered_data = []  # Save the data after all the filtering in this code
    joines_rel = {}  # Save the index relations that we have already joined.
    for i, seq1 in enumerate(sequences):
        counter_no_coincidence = 0  # When seq1 loops through every seq2 and there is no coincidence, we will increase the counter by 1. If the counter is equal to the length of the sequences, we will append the seq1 to the filtered_data, since there is no coincidence with any other seq2.

        for j, seq2 in enumerate(sequences):
            if i == j:  # Skip the same list
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
    filtered_data_NoDup = [value for value in filtered_data if value not in filtered_data]

    data = find_maximal_sets(filtered_data_NoDup)
    return data

