#!/usr/bin/env python
#VDJServer UTSW
#Tool to concatenate several igblast output tsv files

import glob
import csv
import os

IgblastOutput_files = glob.glob("*igblast.out.tsv") 

header_saved = False
with open('merged.tsv','wb') as fout:
    for filename in IgblastOutput_files:
        with open(filename) as fin:
            header = next(fin)
            if not header_saved:
                fout.write(header)
                header_saved = True
            for line in fin:
                fout.write(line)
input = open('merged.tsv', 'rb')
output = open('outfile.tsv', 'wb')
writer = csv.writer(output)
for row in csv.reader(input):
    if row:
        writer.writerow(row)
input.close()
output.close()	
os.remove('merged.tsv')
