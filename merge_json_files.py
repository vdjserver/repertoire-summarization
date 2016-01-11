#!/usr/local/lib python
#VDJServer
#Tool to merge json files

import json
import glob

read_files = glob.glob("*.json")
output_list = []

for f in read_files:
    with open(f, "rb") as infile:
        output_list.append(json.load(infile))

with open("merged_file.json", "wb") as outfile:
    json.dump(output_list, outfile)
