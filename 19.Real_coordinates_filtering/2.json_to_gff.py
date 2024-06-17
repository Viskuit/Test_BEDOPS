# Needed functions
import argparse
import json
import pandas as pd
import os



# ======================================================================
# Defining parse arguments
# ======================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Convert python dictionary JSON file to a GFF file format')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input JSON python dictionary file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output GFF file.')
    return parser.parse_args()

# ======================================================================
# Main function
# ======================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # Read JSON file containing python dict
    with open(args.file, 'r') as file:
        data = json.load(file)
        # In the dict for each element:
            # element[0] = chromosome
            # element[1] = strand
            # element[2] = start coordinate
            # element[3] = end coordinate
            # element[4] = Accepted or Rejected
    counter = 0
    for key, value in data.items():
        for element in value:
            if element[4] == 'Accepted':
                counter += 1
    print(f'Analyzing {counter} elements.')

    # Loop through the python dict to create the dictionary
    pre_gff_list = []
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    for key, value in data.items():
        for i, element in enumerate(value):
            if element[4] == 'Accepted':  # Important, since I only want to print the accepted elements
                pre_gff_list.append(
                    {
                        'seqname': element[0],  # Need here the chromosome id
                        'source': 'CBM',
                        'feature': 'SIDER test',
                        'start': element[2],
                        'end': element[3],
                        'score': '.',
                        'strand': '+',  # Since my data all are in the positive strand
                        'frame': '.',
                        'attribute': f'ID={key}_element_{i}'
                    }
                )
    final_gff_df = pd.DataFrame(pre_gff_list, columns=column_names)
    path_output_folder = os.path.join(args.output, 'GFF_files')
    os.makedirs(path_output_folder, exist_ok=True)
    path_output_file = os.path.join(path_output_folder, 'SIDER_elements.gff')
    final_gff_df.to_csv(path_output_file, sep='\t', index=False, header=False)
    print(f'GFF file created at {path_output_file}')
