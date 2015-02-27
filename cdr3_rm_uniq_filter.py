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
	#gene/allele portion
	ighv_allele=pieces[2]
	ighv_gene=deAllelifyName(ighv_allele)
	rm_text=pieces[185]
	sm_text=pieces[187]
	include_sm=False
	mut_arr=list()
	rms=eval(rm_text)
	#RM portion
	for rm in rms:
		fromPosTo=extractFromPosTo(rm)
		rm_pos=fromPosTo[1]
		rm_to=fromPosTo[2]
		mut_arr.append(rm_pos+rm_to)
	if(include_sm):
		sms=eval(sm_text)
		for sm in sms:
			mut_arr.add(sm)
	#CDR3 portion
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
	#sig=v_sig_port+"."+mut_sig_port+"."+cdr3_sig_port
	sig=mut_sig_port+"."+cdr3_sig_port
	#print "Acquiring a signature...."
	#print "allele=",ighv_allele," gene=",ighv_gene," RM portion =",mut_sig_port," cdr3=",cdr3_sig_port
	#print "To return signature ",sig
	return sig
		









def acquire_sigs_to_max_dict(single_tsv_path):
	reader=open(single_tsv_path,'r')
	line_num=0
	sig_to_max_dict=dict()
	for line in reader:
		if(line_num>0):
			#data!
			temp_line=line.strip()
			sig=acquireSignatureFromTSVLine(temp_line)
			pieces=temp_line.strip()
			count=int(pieces[len(pieces)-1])
			print "For a line the sig is ",sig," and the count is ",count
			if(sig in sig_to_max_dict):
				#do a comparison
				existing_val=sig_to_max_dict[sig]
				if(existing_val<count):
					#count bigger than existing key?, then overwrite!
					sig_to_max_dict[sig]=count
				else:
					#if existing>=count, then leave the dict alone,
					#the first value is the one chosen if the two vals are the same
					pass
			else:
				#nothing? then store it!
				sig_to_max_dict[sig]=count
		else:
			#line_num==0?
			pass
		line_num+=1
	reader.close()
	return sig_to_max_dict






def perfrom_dgx_uniq(tsv_input_path,tsv_output_path,sig_to_max_dict):
	reader=open(tsv_input_path,'r')
	writer=open(tsv_output_path,'w')
	line_num=0
	for line in reader:
		temp_line=line.strip()
		if(line_num==0):
			#write the header!
			writer.write(temp_line+"\n")
		else:
			pieces=temp_line.split('\t')
			this_line_count=int(pieces[len(pieces)-1])
			this_line_sig=acquireSignatureFromTSVLine(temp_line)
			if(not(this_line_sig in sig_to_max_dict)):
				pass
			else:
				this_sig_max=sig_to_max_dict[this_line_sig]
				if(this_sig_max>=this_line_count):
					#output and remove!
					writer.write(temp_line+"\n")
					del sig_to_max_dict[this_line_sig]
				else:
					#skip!
					pass
		line_num+=1
	writer.close()
	reader.close()






#batch process of files
def batch_sig_uniq(input_base,output_base):
	tsv_glob=input_base+"/*/*out.tsv"
	for tsv in glob_walk(tsv_glob):
		#first acquire file level information
		tsv_batch_sample=getBatchIDAndSampleIDFromPath(tsv)
		tsv_batch=tsv_batch_sample[0]
		tsv_sample=tsv_batch_sample[1]
		filtered_file_dir=output_base+"/"+tsv_batch+"/"
		filtered_file_path=filtered_file_dir+os.path.basename(tsv)
		#acquire signature->max data/dict
		print "Now processing file ",tsv," to write to ",filtered_file_path,"..."
		sig_to_max_dict=acquire_sigs_to_max_dict(tsv)
		if(not(os.path.exists(filtered_file_dir))):
			os.makedirs(filtered_file_dir)
		#using the map perform the filtering/uniquing/outputting
		perfrom_dgx_uniq(tsv,filtered_file_path,sig_to_max_dict)



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







