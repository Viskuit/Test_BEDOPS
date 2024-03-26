import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

import os
os.chdir("./Manual_Analysis")

plus_list = [[24093, 137458, 205558], [24758, 138019, 206213]]
minus_list = [[36163, 55775, 114388, 130854], [35313, 54697, 113357, 129878]]

def comparative_csv(path_input, coor, plus_list, minus_list):
    counter = 1
    for value in pd.read_csv(path_input, sep=",", header=None).iterrows():
        value = value[1]

        if coor == "start": 
            coor2 = value[10]
            plus_list = plus_list[0]
            minus_list = minus_list[0]
        elif coor == "end": 
            coor2 = value[11]
            plus_list = plus_list[1]
            minus_list = minus_list[1]

        if value[1] == "LinJ.01" and value[14] == "plus":
            if coor2 in plus_list:
                print(f"{counter} sequence in PLUS found {coor.upper()} in {coor2}")
                counter += 1
        
        elif value[1] == "LinJ.01" and value[14] == "minus":
            if coor2 in minus_list:
                print(f"{counter} sequence in MINUS found {coor.upper()} in {coor2}")
                counter += 1
comparative_csv("../Data/To_analyze/RUNS/run_19.csv", "start", plus_list, minus_list)