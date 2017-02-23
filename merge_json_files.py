#!/usr/bin/env python
#VDJServer UTSW
#Tool to combine several json output files from IgBlast analysis and provide summary results into one file


import json
import glob
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infiles', nargs = '+', help='JSON files to combine')
parser.add_argument('-o', '--outfile', required = False, dest = 'outfile', help='Output filename. If not given, output is sent to stdout.')
args = parser.parse_args()

if args.outfile:
    fout = open(args.outfile,'wb')
else:
    fout = sys.stdout

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
for f in args.infiles:
    with open(f, "rb") as infile:
        json_dict = json.load(infile)

	if (first):
	   first = False
           existing_map = json_dict
           continue
        get_json_terminals_and_counts(json_dict, existing_map)
     
json.dump(existing_map, fout)
