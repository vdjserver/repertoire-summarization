#!/usr/bin/env python

import os
import re
import glob
import argparse
from utils import extractAsItemOrFirstFromList




#program to merge CDR3 histograms for kabat and imgt modes
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Merge multiple output tables, preserving ID')
	parser.add_argument('tables_in',type=str,nargs='+',help="path(s) to tables to merge.  At least one is requried!")
	args=parser.parse_args()
	if(args):
		accum=0
		inFirst=True
		#print "args.table_in ",args.tables_in
		#sys.exit(0)
		for table in args.tables_in:
			#print "Now reading ",table
			line_num=1
			reader=open(table,'r')
			for line in reader:
				line=line.strip()
				if(line_num==1 and inFirst):
					print line
				elif(line_num==1 and not(inFirst)):
					pass
				elif(line_num!=1):
					print "line_num=",line_num
					print "line=",line
					pieces=line.split('\t')					
					pieces[0]+=accum
					joiner="\t"
					line=joiner.join(pieces)
					print line
				line_num+=1
			accum+=(fileLineCount(table)-1)
			inFirst=False
	else:
		parser.print_help()
