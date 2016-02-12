#!/usr/bin/env python
#VDJServer UTSW
#Tool to combine several json output files from IgBlast analysis and provide summary results into one file


import json
import glob
import sys


input_dir = sys.argv[1]
out_file = sys.argv[2]

read_files = glob.glob("./"+input_dir+"/*.json")

def get_json_terminals_and_counts(j_data,existing_map):
    key = j_data['label']
    value = int(j_data['value'])
    existing_map['value'] = existing_map['value'] + value
    
    if('children' in j_data):
        for child in j_data['children']:
            key = child['label']
            for existing_child in existing_map['children']:
                if existing_child['label'] == key: break  
	    get_json_terminals_and_counts(child, existing_child)
    return None

existing_map = None
first = True
for f in read_files:
    with open(f, "rb") as infile:
        json_dict = json.load(infile)

	if (first):
	   first = False
           existing_map = json_dict
           continue
        get_json_terminals_and_counts(json_dict, existing_map)
     
with open("./output/"+out_file, "wb") as outfile:
    json.dump(existing_map, outfile)
