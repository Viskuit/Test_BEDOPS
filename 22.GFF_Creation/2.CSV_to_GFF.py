import argparse
import pandas as pd
import os

# =============================================================================
# Defining argparser arguments
# =============================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Create a GFF file from a CSV file.')
    parser.add_argument('-csv', '--csv_file', type=str, required=True, help='Path to the CSV file with the data.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data.')
    return parser.parse_args()

# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # Read the CSV file
    data = pd.read_csv(args.csv_file, sep=',', header=0)

    # Prepare GFF data
    pre_gff = []
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    for index, row in data.iterrows():
        pre_gff.append(
            {
                'seqname': row['sseqid'],
                'source': 'CBM',
                'feature': 'SIDER',
                'start': row['sstart'],
                'end': row['send'],
                'score': '.',
                'strand': '+',
                'frame': '.',
                'attribute': row['new_IDs']
            }
        )
    
    # Create the GFF file
    final_gff_df = pd.DataFrame(pre_gff, columns=column_names)

    # Folder path
    path_folder = os.path.join(args.output, 'GFF_files')
    os.makedirs(path_folder, exist_ok=True)

    # File path
    path_file = os.path.join(path_folder, 'SIDER_elements.gff')

    # Save the GFF file
    final_gff_df.to_csv(path_file, sep='\t', index=False, header=False)
    print(f'GFF file saved in {path_file}')