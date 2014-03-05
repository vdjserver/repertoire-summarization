#!/usr/bin/env python

import argparse

#program to split IgBLAST output
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Split a single IgBLAST output into multiple files')
	#parser.add_argument('-ways',type=int,nargs=1,help="the number of ways to split the file (integer >=1)")
	parser.add_argument('way_size',type=int,nargs=1,help="the number of IgBLAST records to put in an output file")
	parser.add_argument('input_file',type=str,nargs=1,help="the path to the input file to split")
	parser.add_argument('split_base',type=str,nargs=1,help="the path to the output files (directory must exist, existing files are overwritten) .digits added as suffix")
	parser.print_help()
	args=parser.parse_args()
	if(args):
		pass
	else:
		parser.print_help()
	
