import json
import sys


filename = sys.argv[1]
keys = sys.argv[2:]
with open(filename, "r") as fin:
    d = json.load(fin)

for k in keys:
    d = d[k]
    if isinstance(d, list):
        if len(d) > 1:
            sys.exit("SOMETHING WENT WRONG")
        else:
            d = d[0]

print(d)
