import json
import sys


filename = sys.argv[1]
with open(filename, "r") as fin:
    tasks = json.load(fin)
task = next(task for task in tasks if task["name"] == sys.argv[2])
print(task["id"])
