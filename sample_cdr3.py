#!/usr/bin/env python

from utils import extractAsItemOrFirstFromList,canReadAccess
import argparse
import re

def getBatchIDAndSampleIDFromPath(tsv_path):
	return ["BID","SID"]


def cdr3LenStats(dir_base,out_hist,min_len,max_len):
	file_num=0
	for root, dir, files in os.walk(dir_base):
		for items in fnmatch.filter(files, "*out.tsv"):
			full_tsv_path=root+"/"+items
			reader=open(full_tsv_path,'r')
			line_num=1
			len_hist_this_file=dict()
			for line in reader:
				if(line_num>1):
					pass
					pieces=line.split('\t')
					status=pieces[183]
					if(status=="OK"):
						cdr3_na_len=pieces[21]
						if(cdr3_na_len!="None"):
							cdr3_na_len=int(cdr3_na_len)
							if(cdr3_na_len in len_hist_this_file):
								len_hist_this_file[cdr3_na_len]+=1
							else:
								len_hist_this_file[cdr3_na_len]=1
				line_num+=1
			reader.close()
			writer=open(out_hist,'a')
			if(file_num==0):
				#print headers
				pass
			for val in range(min_len,max_len+1):
								
			file_num+=1




def cdr3DiverStats(dir_base,out_diva):
	for root, dir, files in os.walk(dir_base):
		full_tsv_path=root+"/"+items
		reader=open(full_tsv_path,'r')
		line_num=1
		cdr3_na_this_file=dict()
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
					cdr3_na=pieces[20]
					if(cdr3_na!="None"):
						if(cdr3_na in cdr3_na_this_file):
							cdr3_na_this_file[cdr3_na]+=1
						else:
							cdr3_na_this_file[cdr3_na]=1
			line_num+=1
		reader.close()		
#coal_like.py-96-	#x = {1: 2, 3: 4, 4:3, 2:1, 0:0}
#coal_like.py:97:	#sorted_x = sorted(x.iteritems(), key=operator.itemgetter(1))
#coal_like.py-98-	import operator
#coal_like.py:99:	sorted_x=sorted(bases.iteritems(),key=operator.itemgetter(1))
#coal_like.py-100-	#most frequent base to appear at end of list
#coal_like.py-101-	i_count=0
#coal_like.py:102:	for c in range(len(sorted_x)-1):		
#coal_like.py:103:		i_count+=sorted_x[c][1]
#coal_like.py-104-	return i_count




def getCDR3MinMax(dir_base):
	max_cdr3=None
	min_cdr3=None
	cdr3_len_set=set()
	for root, dir, files in os.walk(dir_base):
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
	







