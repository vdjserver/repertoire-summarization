#!/usr/bin/env python
import argparse
import os
import sys
import re
from utils import extractAsItemOrFirstFromList

#esalina2@eddiecomp:/home/data/big_sim/mouse/split_tr$ head  --lines=20  group_1.fasta.igblast.out.vdjml
#<?xml version="1.0" encoding="UTF-8"?>
#<vdjml xmlns="http://vdjserver.org/xml/schema/vdjml/" version="0.0">
#   <meta>
#      <generator name="libVDJML" version="0.0.9" time_gmt="2014-Apr-28 16:53:25"/>
#      <aligner aligner_id="1" name="IgBLAST" version="2.2.27+" uri="http://www.vdjserver.org" run_id="0">
#         <parameters>IG_BLAST_PARAMETERS</parameters>
#      </aligner>
#      <germline_db gl_db_id="1" name="DB_NAME" species="Mus_musculus" version="DB_VER" uri="http://www.vdjserver.org"/>
#   </meta>
#   <read_results>
#above keep only from first file
#here just copy
#      <read read_id="1|vKey=TRDV5*01|dKey=TRDD1*01|jKey=TRDJ1*01|vd_junc=AAGTAGGCCT|dj_junc=C|chain_type=heavy">
#         <vdj_alig


#keep only  from last file
#   </read_results>
#</vdjml>




def merge(input_list,output_file):
	only_first=["?xml","meta","generator","aligner","parameters","germline_db"]
	straddle=["read_results","vdjml"]
	if(len(input_list)==1):
		reader=open(input_list[0])
		writer=open(output_file,'w')
		for line in reader:
			writer.write(line)
		writer.close()
		reader.close()
	else:
		#print "input multi"
		firstIndex=0
		lastIndex=len(input_list)-1
		writer=open(output_file,'w')
		tagRE=re.compile(r'^\s*<\/?([^\s>]+)')
		closeRE=re.compile(r'^\s*<\s*\/')
		for i in range(len(input_list)):
			if(i==firstIndex or (i==firstIndex and len(input_list)==1)):
				inFirst=True
			else:
				inFirst=False
			if(i==lastIndex or (i==lastIndex and len(input_list)==1)):
				inLast=True
			else:
				inLast=False
			reader=open(input_list[i],'r')
			for line in reader:
				tsr=re.search(tagRE,line)
				passLine=True
				if(tsr):
					tag=tsr.group(1)
					#print "got tag "+str(tag)+" from line "+str(line[0:min(20,len(line))])
					if((tag in only_first) and (inFirst==True)):
						passLine=True
					elif((tag in only_first) and (inFirst==False)):
						passLine=False
					elif(tag in straddle and (inFirst==True or inLast==True)):
						passLine=True
					else:
						passLine=True
					closeREsearch=re.search(closeRE,line)
					if(closeREsearch):
						closeTag=True
					else:
						closeTag=False
					if(tag in straddle):
						#print "got tag ",tag
						if((tag in straddle) and not(inFirst) and not(inLast)):
							passLine=False
						elif((tag in straddle) and closeTag and inFirst):
							passLine=False
						elif((tag in straddle) and not(closeTag) and inLast):
							passLine=False
						elif((tag in straddle) and closeTag and inLast):
							passLine=True

				else:
					print "got no tag from ",line
					sys.exit(0)
				if(passLine):
					writer.write(line)
		writer.close()
			


if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Quick-and-dirty VDJML file merge')
	parser.add_argument('vdjmls',type=str,nargs='+',help="path(s) to VDJMLs to merge.  At least one is requried!")
	parser.add_argument('-output', type=str,nargs=1,default="/dev/stdout",help="Output file path ; default=/dev/stdout")
	args = parser.parse_args()
	if(args):
		#print "good"
		input_list=args.vdjmls
		#print "input is ",input_list
		output_file=extractAsItemOrFirstFromList(args.output)
		#print "output file : ",output_file
		if(len(input_list)>=1):
			merge(input_list,output_file)
			pass
		else:
			#nothing to do!
			pass
	else:
		parser.print_args()




