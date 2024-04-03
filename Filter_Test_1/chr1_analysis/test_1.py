import numpy as np
import pandas as pd
import os
import subprocess

os.chdir("/home/rfpacheco/Desktop/Projects/Testing_Leishmania_project/Filter_Test_1/chr1_analysis")

blast_data = pd.read_csv("./output3.csv", sep=",", header=None)

main_list = []
for index, element in blast_data.iterrows():  # element is a pandas.Series and we can't use boolean comparatives with that.
    main_statement = False
    print(f"===> Index: {index}")
    for index2, group in enumerate(main_list):
        # if element.iloc[0] not in saver and element.iloc[1] not in saver:
        if element.iloc[0] in group and element.iloc[1] not in group:
            group.append(element.iloc[1])
            main_statement = True
            print(f"{element.iloc[0]} added to group {index2}")
            break
        elif element.iloc[1] in group and element.iloc[0] in group:
            main_statement = True
            print(f"SKIPPED {element.iloc[0]} with row[0] {element.iloc[0]} since it's already in group {index2}")
            break

    if main_statement == False:
        main_list.append([element.iloc[0]])
        print(f"New group created with {element.iloc[0]}")
