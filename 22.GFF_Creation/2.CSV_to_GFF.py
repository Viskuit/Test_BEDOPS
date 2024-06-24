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
    parser.add_argument('-n', '--name', type=str, required=True, help='Name of the GFF file.')
    parser.add_argument('--source', type=str, required=False, default='.', help='Source of the data.')
    parser.add_argument('--feature', type=str, required=False, default='.', help='Feature of the data.')
    parser.add_argument('--strand', action='store_true', help='Strand of the data.')
    parser.add_argument('--attribute', type=str, required=False, default='.', help='Column name of the attribute.')
    return parser.parse_args()

# =============================================================================
# Main function
# =============================================================================
if __name__ == '__main__':
    args = parse_arguments()

    # Read the CSV file
    data = pd.read_csv(args.csv_file, sep=',', header=0)
    if (data['sstart'] > data['send']).sum() > 0:
        print(f'There are {(data['sstart'] > data['send']).sum()}  elements with sstart > send. Correct this before continuing.')
        data.loc[data['sstrand'] == 'minus', ['sstart', 'send']] = data.loc[data['sstrand'] == 'minus', ['send', 'sstart']].values  # Swap the values
        print('Values swapped.')
        print(f'Now there are {(data['sstart'] > data['send']).sum()}  elements with sstart > send.')

    # Strand input
    if args.strand == True:
        data.loc[data['sstrand'] == 'minus', 'sstrand'] = '-'
        data.loc[data['sstrand'] == 'plus', 'sstrand'] = '+'
    else:
        pass


    # Prepare GFF data
    pre_gff = []
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    for index, row in data.iterrows():
        pre_gff.append(
            {
                'seqname': row['sseqid'],
                'source': args.source,
                'feature': args.feature,
                'start': row['sstart'],
                'end': row['send'],
                'score': '.',
                'strand': row['sstrand'] if args.strand == True else '.',  # If the strand is False, then the strand will be '.'
                'frame': '.',
                'attribute': f'ID={row[args.attribute]}' if args.attribute != '.' else '.'  # Checking args input.
            }
        )

    
    # Create the GFF file
    final_gff_df = pd.DataFrame(pre_gff, columns=column_names)

    # Folder path
    path_folder = os.path.join(args.output, 'GFF_files')
    os.makedirs(path_folder, exist_ok=True)

    # File path
    path_file = os.path.join(path_folder, args.name)

    # Save the GFF file
    final_gff_df.to_csv(path_file, sep='\t', index=False, header=False)
    print(f'GFF file saved in {path_file}')