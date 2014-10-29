#!/usr/bin/env python

from sample_cdr3 import getBatchIDAndSampleIDFromPath
from utils import extractAsItemOrFirstFromList,canReadAccess
import argparse
import random


def createTempPathFromPath(tsv_path):
	proposed=tsv_path
	while(os.path.exists(proposed)):
		proposed=tsv_path+"."+str(random.random())
	return proposed



def addBatchAndSample(bs_arr,tsv_path,temp_path):
	import os
	if(os.path.exists(temp_path)):
		#don't overwrite the file !
		print "Error, file ",temp_path," exists! Abort!"
		import sys
		sys.exit(0)
	reader=open(tsv_path,'r')
	writer=open(temp_path,'w')
	ln=0
	for line in reader:
		temp_line=line
		temp_line=line.strip()
		pieces=temp_line.split('\t')
		if(ln==0):
			pieces.append("Batch")
			pieces.append("Sample")
		else:
			pieces.append(bs_arr[0])
			pieces.append(bs_arr[1])
		new_line="\t".join(pieces)
		writer.write(new_line+"\n")
		ln+=1
	reader.close()
	writer.close()
	#os.rename(src, dst)
	#overwrite the source with the modified/updated file
	os.rename(temp_path,tsv_path)
	#os.remove(temp_path)

	








if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Add batch/sample information to a TSV file.')
	parser.add_argument('in_file',type=str,nargs=1,help="path to a TSV file (whose name is interpreted as the sample and under a directory interpreted as a batch)")
	args=parser.parse_args()
	import os
	if(not(args)):
		parser.print_help()
		sys.exit(0)
	else:
		#print "Success!!"
		tsv_path=extractAsItemOrFirstFromList(args.in_file)
		if(not(os.path.exists(tsv_path))):
			print "File not found : ",tsv_path
			parser.print_help()
			sys.exit(0)
		bs=getBatchIDAndSampleIDFromPath(tsv_path)
		print bs
		temp=createTempPathFromPath(tsv_path)
		print "temp is ",temp
		addBatchAndSample(bs,tsv_path,temp)














