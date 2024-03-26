import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

import os
os.chdir("./Manual_Analysis")

def fasta_creator(path_input, fasta_output_path):

    data_df = pd.read_csv(path_input, sep=",", header=None)

    fasta_df = pd.DataFrame()
    for index, row in enumerate(data_df.iterrows()):
        rec = SeqRecord(
            Seq(row[5]),
            id="Seq_" + str(index) + "_" + row[0] + "_" + row[4],  # Que tenga aqui el sentido es esencial para luego filtrarlos
            description="Leishmania infantum " + row[4]
        )
        fasta_df.append(rec)

    SeqIO.write(fasta_df, fasta_output_path, "fasta")
    print("\nFasta created at:", fasta_output_path)

fasta_creator("RUN20_LinJ.01.csv", "RUN20_LinJ.01.fasta")