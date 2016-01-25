#!/usr/bin/env python
#VDJServer UTSW
#Tool to combine several json output files from IgBlast analysis and provide summary results into one file

import json
import glob
import sys


input_files = sys.argv[1]
out_file = "json_out.json"

read_files = glob.glob("*.json")
final_json_dict = {}

file_counter = 0
for f in read_files:
    with open(f, "rb") as infile:
        if file_counter == 0:
            final_json_dict = json.load((infile))
        else:
            json_dict = json.load((infile))

            # add the current file's value to final file's value
            final_json_dict["value"] += json_dict["value"]

            children = json_dict["children"]
            final_children = final_json_dict["children"]
            for child in children:
                label = child["label"]

                dup_label = False
                for final_child in final_children:
                    if label in final_child.values():
                        # label is a duplicate
                        final_child["value"] += child["value"]
                        dup_label = True
                        break

                if not dup_label:
                    final_children.append(child)
    file_counter += 1

#print(final_json_dict)
with open(out_file, "wb") as outfile:
    json.dump(final_json_dict, outfile)
