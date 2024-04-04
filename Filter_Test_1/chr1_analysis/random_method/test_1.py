import numpy as np
import pandas as pd
import os

# Count the number of times a group appears for each 1000 files
all_files = os.listdir("./results")
last_list = []
for index, file in enumerate(all_files, 1):
    print(f"Processing file {index}/{len(all_files)}")
    with open(f"./results/{file}", "r") as text_file:
        lines = text_file.readlines()
        for line in lines:
            counter = 0
            if line not in last_list:

                for file2 in all_files:
                    with open(f"./results/{file2}", "r") as text_file2:
                        lines2 = text_file2.readlines()
                        for line2 in lines2:
                            if line == line2:
                                counter += 1

                if counter >= 1:  # How many times it has to appear
                    last_list.append(line)
                    
last_list = [eval(line.strip()) for line in last_list]
len(last_list)

with open(f"./final_subfamilies_once.txt", "w") as file:  # name of the final file
    for group in last_list:
        file.write(f"{group}\n")