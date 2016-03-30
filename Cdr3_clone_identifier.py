#!/usr/bin/env python
#VDJServer UTSW
#Tool to take a vdjserver igblast tsv file and parse it for cdr3 clonality

import csv
import sys
entries = []
duplicate_entries = []
dup_dictionary = {}
value=1
with open(sys.argv[1], 'r') as my_file:
    for line in my_file:
        columns = line.strip().split('\t')
        if columns[12] not in entries:
            entries.append(columns[12])
        else:
            if columns[12] in duplicate_entries: 
                continue
            else:
                duplicate_entries.append(columns[12]) 
                dup_dictionary[columns[12]]=value
                value+=1

with open(sys.argv[2], 'w') as out_file:
    with open(sys.argv[1], 'r') as my_file:
        unique_number=-1
        counter=0
        for line in my_file:
            if counter==0:
                line=line.strip()
                line+='\t'+ 'Clone Identifiers'+'\n'
            else:
                columns = line.strip().split('\t')
                if columns[12] in dup_dictionary.keys():
                    line=line.strip()
                    line+='\t'+ str(dup_dictionary[columns[12]])+'\n'
                else:
                    line=line.strip()
                    line+='\t' + str(unique_number)+'\n'
                    unique_number-=1
                out_file.write(line)
            counter+=1
