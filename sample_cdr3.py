#!/usr/bin/env python

from utils import extractAsItemOrFirstFromList,canReadAccess
import argparse
import re

def cdr3LenStats(path,out_hist):
	reader=open(path,'r')
	line_num=1
	for line in reader:
		#print line
		if(line_num>1):
			pieces=line.split('\t')
			status=pieces[183]
			stat_counts[status]+=1
			if(status=="OK"):
				pass
			else:
				pass
		line_num+=1
	reader.close()



def cdr3DiverStats(path,out_diva):
	pass




def getCDR3MinMax(dir_base):
	max_cdr3=None
	min_cdr3=None
	cdr3_len_set=set()
	for root, dir, files in os.walk(top):
		for items in fnmatch.filter(files, "*out.tsv"):
			full_tsv_path=root+"/"+items
			reader=open(full_tsv_path,'r')
			line_num=1
			for line in reader:
				if(line_num>1):
					pass
					pieces=line.split('\t')
					status=pieces[183]
					if(status=="OK"):
						#19:CDR3 AA (imgt)
						#20:CDR3 AA length (imgt)
						#21:CDR3 NA (imgt)
						#22:CDR3 NA length (imgt)
						cdr3_na_len=pieces[21]
						if(cdr3_na_len!="None"):
							cdr3_na_len=int(cdr3_na_len)
							cdr3_len_set.add(cdr3_na_len)
				line_num+=1
			reader.close()
	max_cdr3=max(cdr3_len_set)
	min_cdr3=min(cdr3_len_set)
	return [min_cdr3,max_cdr3]



if (__name__=="__main__"):
	parser.add_argument('tsv_in',type=str,nargs=1,help="path to a rep_char TSV output file")
	parser.add_argument('bid',type=str,nargs=1,help="string for BATCH ID")
	parser.add_argument('sid',type=str,nargs=1,help="string for SAMPLE ID")
	parser.add_argument('out_hist',type=str,nargs=1,help="output file path for TSV of sample/CDR3 length table")
	parser.add_argument('out_diva',type=str,nargs=1,help="output file path for TSV of CDR3 length diversity analysis")
	parser.add_argument('dir_base',type=str,nargs=1,help="base dir for CDR3 length hist")
	args=parser.parse_args()
	if(not(args)):
		parser.print_help()
	input_file_path=extractAsItemOrFirstFromList(args.tsv_in)
	if(not(canReadAccess(input_file_path))):
		parser.print_help()
	bid=extractAsItemOrFirstFromList(args.bid)
	sid=extractAsItemOrFirstFromList(args.sid)
	out_hist=extractAsItemOrFirstFromList(args.out_hist)
	out_diva=extractAsItemOrFirstFromList(args.out_diva)
	







