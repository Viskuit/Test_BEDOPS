import argparse
import os
import pandas as pd
import re


# =============================================================================
# Defining "argparse" arguments
# =============================================================================
def parse_arguments():
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description='From the relationship database CSV and the original FASTA file, divide the data in two subsets: related and unrelated')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to the NAMED relationship database CSV file')
    parser.add_argument('-csv2', '--csv_file2', type=str, required=True, help='Path to the original CSV file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data')
    return parser.parse_args()

# =============================================================================
# Defining needed functions
# =============================================================================
def pattern_search(input_string):
    
    # Define the regular expression pattern to extract the number
    pattern = r"Seq_(\d+)_LinJ"
    
    # Search for the pattern in the input string
    match = re.search(pattern, input_string)
    
    # Extract the number if the pattern is found
    if match:
        return int(match.group(1))
    else:
        return "Pattern not found"
    
def extract_number(df):
    main_list = []
    for i, row in df.iterrows():
        for elem in row:
            if pd.notna(elem):
                main_list.append(pattern_search(elem))
    return main_list

# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # Get the working folder
    path_working_folder = os.path.join(args.output, 'relationship_database')
    os.makedirs(path_working_folder, exist_ok=True)

    # Get the relationship database csv
    relationship_data = pd.read_csv(args.csv_file, sep=',', header=None)

    # Get the original csv
    og_data = pd.read_csv(args.csv_file2, sep=',', header=0)

    # Get lines with rows with only one element
    elem_counter = relationship_data.notna().sum(axis=1)

    # Get the rows with only one element
    one_elem_rows = relationship_data[elem_counter == 1]
    more_elem_rows = relationship_data[elem_counter > 1]

    # For each row like "Seq_10_LinJ.04" get the number "10"
    one_elem_index = extract_number(one_elem_rows)

    # The same for the "more_elem_rows" but for every column which is not NaN
    more_elem_index = extract_number(more_elem_rows)
    
    #

    # use "one_elem_rows" and "more_elem_rows" to get the related and unrelated data from the original csv since these numbers are the indices from the original csv
    related_data = og_data.iloc[more_elem_index]
    unrelated_data = og_data.iloc[one_elem_index]

    # Save the data
    related_data.to_csv(os.path.join(path_working_folder, 'negative_nomatchCDS_related_data.csv'), index=False, header=True)
    unrelated_data.to_csv(os.path.join(path_working_folder, 'negative_nomatchCDS_unrelated_data.csv'), index=False, header=True)