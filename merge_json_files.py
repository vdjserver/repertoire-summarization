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
    if(existing_map is None):
        existing_map=dict()
    if('children' in j_data):
        for child in j_data['children']:
            key = j_data['label']
            value = int(j_data['value'])
            existing_map[key] = value
            existing_map = get_json_terminals_and_counts(child, existing_map)
    elif(('label' in j_data) and ('value' in j_data)):
        key = j_data['label']
        value = int(j_data['value'])
        existing_map[key] = value
        return existing_map
    else:
        raise Exception("WARNING, either no children or no label/value found in ",j_data,"!")
    return existing_map

final_json_dict = {}
final_json_dict["gene_segments"] = {}
existing_maps = []
for f in read_files:
    with open(f, "rb") as infile:
        json_dict = json.load(infile)

        first_label = json_dict["label"]
        first_value = json_dict["value"]
        if first_label in final_json_dict.keys():
            final_json_dict[first_label] += first_value
        else:
            final_json_dict[first_label] = first_value

        existing_map = None
        for child in json_dict["children"]:
            existing_map = get_json_terminals_and_counts(child, existing_map)
            existing_maps.append(existing_map)
        # print(existing_map)
for em in existing_maps:
    for key in em.keys():
        if key in final_json_dict["gene_segments"].keys():
            final_json_dict["gene_segments"][key] += em[key]
        else:
            final_json_dict["gene_segments"][key] = em[key]


with open("./output/"+out_file, "wb") as outfile:
    json.dump(final_json_dict, outfile)