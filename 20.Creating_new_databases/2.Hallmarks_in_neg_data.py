# Needed modules
import argparse
import pandas as pd
import os
import subprocess

# Importing from other .py files in this dict


# ======================================================================
# Defining parse arguments
# ======================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Check fasta hits in the a CSV data base.')
    parser.add_argument('-f1', '--fasta1', type=str, required=True, help='Path to the CSV file in fasta format where we want to search.')
    parser.add_argument('-f2', '--fasta2', type=str, required=True, help='Path to the fasta file with the Hallmarks to search.')
    parser.add_argument('-f3', '--csv', type=str, required=True, help='Path to the CSV file with the data base to search, same as --fasta1 but in csv format.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data.')
    # parser.add_argument('-db', '--database', type=str, required=True, help='Path to the genome fasta file.')
    return parser.parse_args()
# ======================================================================
# Defining neeeded functions
# ======================================================================
def blastn_dic(path_input, path_output):
    cmd = f'makeblastdb -in {path_input} -dbtype nucl -parse_seqids -out {path_output}'
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def blastn_blaster(query, path_genome, evalue):
    cmd = (
        f'blastn -word_size 11 '
        f'-query {query} '
        f'-db {path_genome} '
        f'-evalue {evalue} '
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


# ======================================================================
# Main function
# ======================================================================
if __name__ == '__main__':
    args = parse_arguments()  # Parse arguments

    # Create the blast database
    path_dict_folder = os.path.join(args.output, 'blast_dict')
    os.makedirs(path_dict_folder, exist_ok=True)
    path_blast_dict_file = os.path.join(path_dict_folder, os.path.basename(args.fasta1))
    blastn_dic(args.fasta1, path_blast_dict_file)

    # Run the blastn
    hallmark_match_df = blastn_blaster(query=args.fasta2, path_genome=path_blast_dict_file, evalue=0.001)
    # parser.add_argument('-db', '--database', type=str, required=True, help='Path to the genome fasta file.')

    # Save the data (just to analyze later; opcional)
    path_hallmark_df = os.path.join(args.output, 'hallmark_hits.csv')
    hallmark_match_df.to_csv(path_hallmark_df, index=False, header=True, sep=',')

    # Extract the index of the hits, because in 'sseqid' it will be sth like Seq_23_LinJ.32 and I want to extract the number 23.
    hallmark_match_df['index'] = hallmark_match_df['sseqid'].str.extract(r'Seq_(\d+)_')
    hallmark_match_df['index'] = hallmark_match_df['index'].astype(int)

    # Get the index list, sort it and remove duplicates
    index_list = hallmark_match_df['index'].sort_values().unique().tolist()

    # Check wich elements from that list are in the original CSV file
    csv_df = pd.read_csv(args.csv, sep=',', header=0)
    recaught_df = csv_df[csv_df.index.isin(index_list)]  # This will return the elements that are in the index list, and this elements are the evalue significative hits from the hallmark hits

    # Remove those elements from the original CSV file
    csv_df = csv_df[~csv_df.index.isin(index_list)]

    # Save the recaught data
    path_recaught_df = os.path.join(args.output, 'neg_data_recaught_hits.csv')
    recaught_df.to_csv(path_recaught_df, index=False, header=True, sep=',')

    # Save the new neg data
    path_new_neg_df = os.path.join(args.output, 'negative_data_after_recaught.csv')
    csv_df.to_csv(path_new_neg_df, index=False, header=True, sep=',')