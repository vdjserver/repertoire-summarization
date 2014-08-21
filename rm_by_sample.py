#!/usr/bin/env python


from utils import extractAsItemOrFirstFromList,canReadAccess
import fnmatch
import argparse
import re
import os
import operator

def initNMCountMap():
	keys=getNumPosArr()
	cm=dict()
	for k in keys:
		cm[k]=0
	return cm

def getNumPosArr():
	np_list=[
	"31",
	"31A",
	"31B",
	"32",
	"33",
	"34",
	"35",
	"36",
	"37",
	"38",
	"39",
	"40",
	"41",
	"42",
	"43",
	"44",
	"45",
	"46",
	"47",
	"48",
	"49",
	"50",
	"51",
	"52",
	"52A",
	"52B",
	"52C",
	"53",
	"54",
	"55",
	"56",
	"57",
	"58",
	"59",
	"60",
	"61",
	"62",
	"63",
	"64",
	"65",
	"66",
	"67",
	"68",
	"69",
	"70",
	"71",
	"72",
	"73",
	"74",
	"75",
	"76",
	"77",
	"78",
	"79",
	"80",
	"81",
	"82",
	"82A",
	"82B",
	"82C",
	"83",
	"84",
	"85",
	"86",
	"87",
	"88",
	"89",
	"90",
	"91",
	"92"]
	return np_list


def writeRMInfoToTable(dir_base):
	file_num=0
	for root, dirs, files in os.walk(dir_base):
		for items in fnmatch.filter(files, "*out.tsv"):
			cdr1_lens=[0,0,0]
			full_tsv_path=root+"/"+items
			#print "Analyzing file ",full_tsv_path," for CDR3 min/max length...."
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			#print "The batch and sample :",bs
			reader=open(full_tsv_path,'r')
			line_num=1
			rm_count_map=initNMCountMap()
			for line in reader:
				if(line_num>1):
					pieces=line.split('\t')
					status=pieces[183]
					if(status=="OK"):
						CDR1_len=int(pieces[110])
						cdr1_amino_len=CDR1_len/3
						cdr1_lens[cdr1_amino_len-5]+=1
						rms=eval(pieces[185])
						for rm in rms:
							rm_str_len=len(rm)
							if(rm_str_len==4):
								numbered_pos=rm[1:3]
							elif(rm_str_len==5):
								numbered_pos=rm[1:4]
							else:
								raise Exception("Error, unknown number (extraction pos) in  numbered RM : "+rm)
							if(numbered_pos in rm_count_map):
								rm_count_map[numbered_pos]+=1
							else:
								raise Exception("Error, unseen number (extraction pos) in  numbered RM : "+rm)
								
				line_num+=1
			reader.close()
			if(file_num==0):
				#put header
				pass
				



if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Print DIOGENIX formatted rep_char CDR3 statistics')
	parser.add_argument('out_tbl',type=str,nargs=1,help="output file path for TSV of sample/CDR3 length table")
	parser.add_argument('dir_base',type=str,nargs=1,help="base dir for CDR3 length hist")
	args=parser.parse_args()
	if(not(args)):
		parser.print_help()
	else:
		out_hist=extractAsItemOrFirstFromList(args.out_hist)
		out_diva=extractAsItemOrFirstFromList(args.out_diva)
		dir_base=extractAsItemOrFirstFromList(args.dir_base)
		if(not(os.path.exists(dir_base))):
			parser.print_help()	





