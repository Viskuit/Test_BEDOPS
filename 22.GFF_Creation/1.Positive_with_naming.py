import argparse
import pandas as pd
import os
import re

# =============================================================================
# Defining argparser arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Create new CSV of positive data with the new names')
    parser.add_argument('-p', '--main_path', type=str, required=True, help='Path to the main folder where the named and no named data is.')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to the CSV file with the positive SIDER elements.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data.')
    return parser.parse_args()

# =============================================================================
# Defining needed functions
# =============================================================================



# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # get the working folder
    os.makedirs(args.output, exist_ok=True)

    # Get the folders name
    folders = [folder for folder in os.listdir(args.main_path) if os.path.isdir(os.path.join(args.main_path, folder))]
    folders = sorted(folders, key=lambda x: int(x.split('.')[1]))

    # Get the files with the ID form the main data base and the NAMED data base. This gets the path
    files_dirs = []
    files_dirs_named = []
    for folder in folders:
        path_chr_folder = os.path.join(args.main_path, folder)
        file = [os.path.join(path_chr_folder, file) for file in os.listdir(path_chr_folder) if file.endswith(".csv") and file.startswith("families")]  # There will be only one file in the folder
        file_named = [os.path.join(path_chr_folder, file) for file in os.listdir(path_chr_folder) if file.endswith(".csv") and file.startswith("NAMED_families")]
        files_dirs.append(file[0])
        files_dirs_named.append(file_named[0])
    
    # Get dict with the IDs
    dict_ids = {}
    for file, file_named in zip(files_dirs, files_dirs_named):
        file_df = pd.read_csv(file, sep=',', header=None)
        file_named_df = pd.read_csv(file_named, sep=',', header=None)
        for (index, row), (index_named, row_named) in zip(file_df.iterrows(), file_named_df.iterrows()):  # Iterate over the rows of the dataframes
            row = [value for value in row if str(value) != 'nan']  # Remove the nan values
            row_named = [value for value in row_named if str(value) != 'nan']  # Remove the nan values
            for value, value_named in zip(row, row_named):
                slice_id = re.search(r'_(\d+)_', value).group(1)  # Get the slice ID
                slice_id = int(slice_id)  # Convert to integer
                dict_ids[slice_id] = [value, value_named]  # Add the ID and the named ID to the dictionary
    
    # Sort the dictionary by the key (its an integer) in ascending order
    dict_ids = dict(sorted(dict_ids.items()))

    # Now that the dict is sorted, create a list with the values[1], which are the named ones.
    named_values = [value[1] for value in dict_ids.values()]

    # Get the main data frame
    main_df = pd.read_csv(args.csv_file, sep=',', header=0)

    # Add a row with the new names
    main_df['new_IDs'] = named_values

    # Save the new CSV
    path_file_out = os.path.join(args.output, "positive_data_named.csv")
    main_df.to_csv(path_file_out, index=False)  # Save the new CSV
    print(f'New CSV saved in {path_file_out}')
            