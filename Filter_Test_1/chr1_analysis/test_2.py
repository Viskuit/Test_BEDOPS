import numpy as np
import pandas as pd
import os
import subprocess

os.chdir("/home/rfpacheco/Desktop/Projects/Testing_Leishmania_project/Filter_Test_1/chr1_analysis")

blast_data = pd.read_csv("./output3.csv", sep=",", header=None)

def grouping_elements(main_data):
    main_list = []
    for index, element in main_data.iterrows():
        main_statement = False
        element_counter = 0
        for index2, group in enumerate(main_list):
            if element.iloc[1] in group:
                element_counter += 1
            if element.iloc[0] in group and element_counter == 0:
                group.append(element.iloc[1])
                main_statement = True
                break
            elif element_counter > 0 and element.iloc[0] in group:
                main_statement = True
                break

        if main_statement == False:
            main_list.append([element.iloc[0]])
    return main_list

for i in range(1,101):
    shuffled_blast_data = blast_data.sample(frac=1).reset_index(drop=True)
    print(shuffled_blast_data.iloc[0][0])
    
    main_list = grouping_elements(shuffled_blast_data)

    with open(f"./results/sufamilies_test_{i}.txt", "w") as file:
        for group in main_list:
            file.write(f"{group}\n") 

