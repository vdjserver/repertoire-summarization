#!/usr/bin/env python
#VDJServer UTSW
#Tool to combine several json output files from IgBlast analysis and provide summary results into one file

import json
import glob
import sys

input_dir = sys.argv[1]
out_file = sys.argv[2]

read_files = glob.glob("./"+input_dir+"/*.json")

final_json_dict = {}
final_json_dict["gene_segments"] = {}
for f in read_files:
    with open(f, "rb") as infile:
        json_dict = json.load(infile)

        first_label = json_dict["label"]
        first_value = json_dict["value"]
        if first_label in final_json_dict.keys():
            final_json_dict[first_label] += first_value
        else:
            final_json_dict[first_label] = first_value
        while "children" in json_dict["children"][0].keys():
            first_label = json_dict["children"][0]["label"]
            first_value = json_dict["children"][0]["value"]
            if first_label in final_json_dict.keys():
                final_json_dict["gene_segments"][first_label] += first_value
            else:
                final_json_dict["gene_segments"][first_label] = first_value
            json_dict = json_dict["children"][0]
        last_label = json_dict["children"][0]["label"]
        last_value = json_dict["children"][0]["value"]
        if last_label in final_json_dict.keys():
            final_json_dict["gene_segments"][last_label] += last_value
        else:
            final_json_dict["gene_segments"][last_label] = last_value


with open("./output/"+out_file, "wb") as outfile:
    json.dump(final_json_dict, outfile)