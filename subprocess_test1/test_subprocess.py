import subprocess

# subprocess.run(["bedops", "--merge", "df1_plus_sorted.bed", ">", "df1_plus_sorted_merged.bed"])

subprocess.run(["bedops", "--merge", "df1_plus_sorted.bed"], stdout=open("df1_plus_sorted_merged.bed", "w"))