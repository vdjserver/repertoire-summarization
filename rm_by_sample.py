#!/usr/bin/env python


from utils import extractAsItemOrFirstFromList,canReadAccess
import fnmatch
import argparse
import re
import os
import operator
from sample_cdr3 import getBatchIDAndSampleIDFromPath

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


def writeRMInfoToTable(dir_base,out_tbl,num_rm_filter_thresh="0+"):
	file_num=0
	validation_res=validate_filter_stat_arg(num_rm_filter_thresh)
	thresh_num=validation_res[0] #a non-negative integer
	num_plus=validation_res[1] #True or False
	for root, dirs, files in sorted(os.walk(dir_base)):
		for items in sorted(fnmatch.filter(files, "*out.tsv")):
			cdr1_lens=[0,0,0]
			num_pass_filter=0
			full_tsv_path=root+"/"+items
			print "Analyzing file ",full_tsv_path,"..."
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			print "The batch and sample :",bs
			reader=open(full_tsv_path,'r')
			line_num=1
			rm_count_map=initNMCountMap()
			for line in reader:
				if(line_num>1):
					pieces=line.split('\t')
					status=pieces[183]
					if(status=="OK"):
						num_pass_filter+=1
						CDR1_len=int(pieces[110])
						cdr1_amino_len=CDR1_len/3
						cdr1_lens[cdr1_amino_len-5]+=1
						rms=eval(pieces[185])
						num_rms_in_this_row=len(rms)
						examine_counts_for_this_row=False
						if(num_plus and num_rms_in_this_row>=thresh_num):
							examine_counts_for_this_row=True
						if(not(num_plus) and num_rms_in_this_row==thresh_num):
							examine_counts_for_this_row=True
						if(examine_counts_for_this_row):
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
						else:
							#undo these counts!
							cdr1_lens[cdr1_amino_len-5]-=1
							num_pass_filter-=1
				line_num+=1
			reader.close()
			nPos_arr=getNumPosArr()
			if(file_num==0):
				writer=open(out_tbl,'w')
				writer.write("BATCH\tSAMPLE\tCDR1_5_NUM\tCDR1_6_NUM\tCDR1_7_NUM\tNUM_PASS")
				for n in range(len(nPos_arr)):
					writer.write("\t"+nPos_arr[n])
				writer.write("\n")
				writer.close()
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			writer=open(out_tbl,'a')
			batch=bs[0]
			sample=bs[1]
			writer.write(batch+"\t")
			writer.write(sample)
			for c in range(len(cdr1_lens)):
				writer.write("\t"+str(cdr1_lens[c]))
			writer.write("\t"+str(num_pass_filter))
			for npi in range(len(nPos_arr)):
				writer.write("\t"+str(rm_count_map[nPos_arr[npi]]))
			writer.write("\n")
			writer.close()
			file_num+=1
			
			


def validate_filter_stat_arg(v):
	import re
	valid_re=re.compile(r'^(\d+)(\+?)$')
	search_res=re.search(valid_re,v)
	if(search_res):
		num=int(search_res.group(1))
		plus=search_res.group(2)
		print "num=",num," plus=",plus
		if(plus is None):
			return [num,False]
		else:
			if(plus=="+"):
				return [num,True]
			else:
				return [num,False]
	else:
		return None



if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Print DIOGENIX formatted RM frequency statistics (exact or exact+ modes)')
	parser.add_argument('out_tbl',type=str,nargs=1,help="output file path for TSV of sample/RM data table")
	parser.add_argument('dir_base',type=str,nargs=1,help="base dir for files")
	parser.add_argument('-exact_plus',type=str,nargs=1,default="0+",help="valid range to report over : may be : 0,0+,1,1+,2,2+, etc.... depending on range")
	args=parser.parse_args()
	if(not(args)):
		parser.print_help()
	else:
		exact_val=extractAsItemOrFirstFromList(args.exact_plus)
		exact_val_validated=validate_filter_stat_arg(exact_val)
		if(exact_val_validated is None):
			print "Invalid 'exact_plus' : '"+str(exact_val)+"' "
			parser.print_help()
			import sys
			sys.exit(0)
		out_tbl=extractAsItemOrFirstFromList(args.out_tbl)
		dir_base=extractAsItemOrFirstFromList(args.dir_base)
		if(not(os.path.exists(dir_base))):
			parser.print_help()
		else:
			writeRMInfoToTable(dir_base,out_tbl,exact_val)




