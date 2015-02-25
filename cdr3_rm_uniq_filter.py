#!/usr/bin/env python

from utils import glob_walk,extractAsItemOrFirstFromList,deAllelifyName,extractFromPosTo
import os
from sample_cdr3 import getBatchIDAndSampleIDFromPath
from Bio import SeqIO
from ags_mgr import ags_manager
import argparse
import re






#acquire signature from a line
def acquireSignatureFromTSVLine(l,use_na_cdr3=False):
	pieces=l.split('\t')
	ighv_allele=pieces[2]
	ighv_gene=deAllelifyName(ighv_allele)
	rm_text=pieces[185]
	sm_text=pieces[187]
	include_sm=False
	mut_arr=list()
	rms=eval(rm_text)
	for rm in rms:
		fromPosTo=extractFromPosTo(rm)
		rm_from=fromPosTo[0]
		rm_pos=fromPosTo[1]
		rm_to=fromPosTo[2]
		mut_arr.append(rm_pos+rm_to)
	if(include_sm):
		sms=eval(sm_text)
		for sm in sms:
			mut_arr.add(sm)
	if(use_na_cdr3):
		cdr3_index=20
	else:
		cdr3_index=18
	cdr3=pieces[cdr3_index]
	if(cdr3==None or cdr3=="None"):
		#no cdr3?!
		pass
	#mutation portion of signature
	mut_sig_port="".join(mut_arr)
	#cdr3 portion of signature
	cdr3_sig_port=cdr3
	#V call portion of signature
	v_sig_port=ighv_gene
	sig=v_sig_port+"."+mut_sig_port+"."+cdr3_sig_port
	return sig
		









#return uniq lines from TSV based on signature
def acquire_uniq_lines_from_path(p,include_header_in_output=True):
	reader=open(p,'r')
	line_num=0
	sig_to_data_map=dict()
	output=list()
	for line in reader:
		temp_line=line.strip()
		if(line_num==0):
			header=temp_line
			output.append(header)
		else:
			#pieces=temp_line.split('\t')
			#status=pieces[183]	
			#status presumed to ALWAYS be OKAY here????
			sig=acquireSignatureFromTSVLine(temp_line)
			if(sig in sig_to_data_map):
				pass
			else:
				output.append(temp_line)
			sig_to_data_map[sig]=temp_line
		line_num+=1
	return output




#batch process of files
def batch_sig_uniq(input_base,output_base):
	tsv_glob=input_base+"/*/*out.tsv"
	for tsv in glob_walk(tsv_glob):
		tsv_batch_sample=getBatchIDAndSampleIDFromPath(tsv)
		tsv_batch=tsv_batch_sample[0]
		tsv_sample=tsv_batch_sample[1]
		filtered_file_dir=output_base+"/"+tsv_batch+"/"
		filtered_file_path=filtered_file_dir+os.path.basename(tsv)
		filtered_lines=acquire_uniq_lines_from_path(tsv)
		print "Data read from ",tsv,"...now to write filtered data to",filtered_file_path,"..."
		if(not(os.path.isdir(filtered_file_dir)) and not(os.path.exists(filtered_file_dir))):
			os.makedirs(filtered_file_dir)
		writer=open(filtered_file_path,'w')
		for i in range(len(filtered_lines)):
			writer.write(filtered_lines[i]+"\n")
		writer.close()
	



#program to merge CDR3 histograms for kabat and imgt modes
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Perform new unique filtering (based on CDR3 and RM signature and VH4)')
	parser.add_argument('tsv_base',help="path to a directory containing samples and batches")
	parser.add_argument('output_base',help="path to a non-existent output directory")
	args=parser.parse_args()
	if(args):
		tsv_base=extractAsItemOrFirstFromList(args.tsv_base)
		#glob_walk(glob_str,sort_it=True):
		if(not(os.path.exists(tsv_base)) or not(os.path.isdir(tsv_base))):
			print "Error, tsv_base ",tsv_base," either not found or not a directory ! Abort !"
			parser.print_help()
		else:
			output_base=extractAsItemOrFirstFromList(args.output_base)
			if(os.path.exists(output_base)):
				print "Error, output path ",output_base," already exists! Abort !"
				parser.print_help()
			else:
				print "ready to process!"
				print "To read data from ",tsv_base,"..."
				batch_sig_uniq(tsv_base,output_base)
	else:
		parser.print_help()







