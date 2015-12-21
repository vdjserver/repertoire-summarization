#!/usr/bin/env python
#VDJServer UTSW
#Tool to concatenate several igblast output tsv files

import sys
import os
import glob
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('infiles', nargs = '+', help='TSV files to combine')
parser.add_argument('-o', '--outfile', required = False, dest = 'outfile', help='Output filename. If not given, output is sent to stdout.')
args = parser.parse_args()

if args.outfile:
    fout = open(args.outfile,'wb')
else:
    fout = sys.stdout

header_saved = False
for filename in args.infiles:
    with open(filename, 'rb') as fin:
        header = fin.readline()
        if not header_saved:
            fout.write(header)
            header_saved = True
        for line in fin:
            if line:
                fout.write(line)
if fout is not sys.stdout:
    fout.close()

