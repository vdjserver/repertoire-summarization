#!/usr/bin/env python

import argparse
import sys
import os

def split(infile,way_size,split_base):
	INPUT=open(infile,'r')
	o_num=(-1)#output file number
	c_num=(-1)#current record number
	writer=None
	for line in INPUT:
		if(line.startswith("# Query:")):
			c_num+=1
		if((c_num%way_size)==0 or c_num==0):
			if(writer is not None):
				writer.close()
			o_num+=1
			out_file_path=split_base+"/split."+str(o_num)
			print "Opening file #",str(o_num),"...",out_file_path
			writer=open(out_file_path,'w')
		if(writer is not None):
			writer.write(line)
	if(not writer is None):
		writer.close()



#program to split IgBLAST output
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Split a single IgBLAST output into multiple files')
	#parser.add_argument('-ways',type=int,nargs=1,help="the number of ways to split the file (integer >=1)")
	parser.add_argument('way_size',type=int,nargs=1,default=1000,help="the number of IgBLAST records to put in an output file")
	parser.add_argument('input_file',type=str,nargs=1,help="the path to the input file to split")
	parser.add_argument('split_base',type=str,nargs=1,help="the path to the output files (directory must exist, existing files are overwritten) .digits added as suffix")
	args=parser.parse_args()
	if(args):
		if(not(type(args.input_file)==list)):
			print "Invalid input file!"
			parser.print_help()
			sys.exit(1)
		infile=args.input_file[0]
		if(not(os.path.isfile(infile))):
			print "Invalid input file, it's not a file!"
			parser.print_help()
			sys.exit(1)
		if(not(type(args.way_size)==list)):
			print "Invalid way_size!"
			parser.print_help()
			sys.exit(1)
		way_size=args.way_size[0]
		if(way_size<0):
			print "Error, got negative way_size!  Require positive way_size! Exit !"
			parser.print_help()
			sys.exit(1)
		if(not(type(args.split_base)==list)):
			print "Error, invalid split_base ! Abort !"
			parser.print_help()
			sys.exit(1)
		split_base=args.split_base[0]
		if(not(os.path.isdir(split_base))):
			print "Error, require split_base to be a directory ! Abort !"
			parser.print_help()
			sys.exit(1)
		#print "got infile "+infile
		#print "got way size",way_size
		#print "got split_base",split_base
		split(infile,way_size,split_base)
	else:
		parser.print_help()











	
