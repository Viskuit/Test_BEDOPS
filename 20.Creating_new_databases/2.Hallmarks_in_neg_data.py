# Needed modules
import argparse
import pandas as pd
import os

# Importing from other .py files in this dict


# ======================================================================
# Defining parse arguments
# ======================================================================
def parse_arguments():
    parser = argparse.ArgumentParser(description='Check fasta hits in the a CSV data base.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the CSV file that contains the data.')
    parser.add_argument('--fasta', type=str, required=True, help='Path to the fasta file with the Hallmarks to search.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output directory to save the data.')
    return parser.parse_args()
# ======================================================================
# Defining neeeded functions
# ======================================================================




# ======================================================================
# Main function
# ======================================================================
if __name__ = '__main__':
