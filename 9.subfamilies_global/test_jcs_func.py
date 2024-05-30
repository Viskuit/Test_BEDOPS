import pandas as pd
from subfamilies_global_functions import join_conflicted_sequences

test_data = [
    [807, 808, 809],
    [807, 809, 810]
]

test_data2 = [
    [807, 808, 809],
    [807, 809, 810],
    [700, 710, 720, 810]
]

test_data3 = [
    [807, 808, 809],
    [807, 809, 810],
    [10, 20, 30, 40],
    [700, 710, 720, 810],
    [40, 50, 60],
    [1000, 2000, 3000, 4000],
    [10000, 20000, 30000, 40000],
    [2000, 6000, 8000, 10000]
]

data, seqs = join_conflicted_sequences(test_data3)
print(f"Data is: {data}")
print("=" *20)
print(f"Second data is {seqs}")