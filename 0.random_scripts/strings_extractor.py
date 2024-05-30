import re
import sys

string = sys.argv[1]

# Select numbers between the caracters "_"
result = re.findall(r'_(\d+)_', string)

# print results without the "[]" and the " "" "
result2 = ", ".join(result)
print(result2)