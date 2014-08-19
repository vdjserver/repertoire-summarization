#!/usr/bin/env python

from utils import extractAsItemOrFirstFromList,canReadAccess
import fnmatch
import argparse
import re
import os
import operator

def getBatchIDAndSampleIDFromPath(tsv_path):
	bsre=re.compile(r'/([^/]+)/(Sample[^\.]+)\.fasta')
	re_res=re.search(bsre,tsv_path)
	if(re_res):
		batch=re_res.group(1)
		sample=re_res.group(2)
		return [batch,sample]
	else:
		return ["BID_ERR","SID_ERR"]


def cdr3LenStats(dir_base,out_hist,min_len,max_len):
	file_num=0
	for root, dirs, files in sorted(os.walk(dir_base)):
		for items in sorted(fnmatch.filter(files, "*out.tsv")):
			full_tsv_path=root+"/"+items
			reader=open(full_tsv_path,'r')
			line_num=1
			len_hist_this_file=dict()
			for line in reader:
				if(line_num>1):
					pieces=line.split('\t')
					status=pieces[183]
					if(status=="OK"):
						cdr3_na_len=pieces[19]
						if(cdr3_na_len!="None" and cdr3_na_len!=""):
							cdr3_na_len=int(cdr3_na_len)
							if(cdr3_na_len in len_hist_this_file):
								len_hist_this_file[cdr3_na_len]+=1
							else:
								len_hist_this_file[cdr3_na_len]=1
				line_num+=1
			reader.close()
			if(file_num==0):
				#print headers
				writer=open(out_hist,'w')
				writer.write("batch\tsample")
				for val in range(min_len,max_len+1):
					writer.write("\tCDR3_LEN_"+str(val))
				writer.write("\n")
				writer.close()
			writer=open(out_hist,'a')
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			print "Writing CDR3 len/hist data for",bs,"to",out_hist
			writer.write(bs[0]+"\t"+bs[1])
			for val in range(min_len,max_len+1):
				if(val in len_hist_this_file):
					writer.write("\t"+str(len_hist_this_file[val]))
				else:
					writer.write("\t"+str(0))
			writer.write("\n")
			writer.close()
			file_num+=1




def cdr3DiverStats(dir_base,out_diva):
	file_num=0
	for root, dirs, files in sorted(os.walk(dir_base)):
		for items in sorted(fnmatch.filter(files, "*out.tsv")):
			full_tsv_path=root+"/"+items
			reader=open(full_tsv_path,'r')
			line_num=1
			cdr3_na_this_file=dict()
			denominator=0
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
							denominator+=1
							if(cdr3_na in cdr3_na_this_file):
								cdr3_na_this_file[cdr3_na]+=1
							else:
								cdr3_na_this_file[cdr3_na]=1
				line_num+=1
			reader.close()
			sorted_lens=sorted(cdr3_na_this_file.iteritems(),key=operator.itemgetter(1))
			sorted_lens.reverse()
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			#print "Analysis for",bs
			#print "sorted_lens :",sorted_lens
			#print "sorted_lens[0]=",sorted_lens[0]
			#dna=sorted_lens[0][0]
			#dna_count=sorted_lens[0][1]
			#print "dna and count are : ",dna,"<>",dna_count
			#print "\n\n\n\n\n"
			if(file_num==0):
				#print headers
				writer=open(out_diva,'w')
				writer.write("batch\tsample\tNUM_UNIQ_CDR3\tCDR3_DIVERSITY_VECTOR")
				writer.write("\n")
				writer.close()
			writer=open(out_diva,'a')
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			num_keys=len(sorted_lens)
			print "Writing CDR3 diversity data for",bs,"to ",out_diva
			writer.write(bs[0]+"\t"+bs[1]+"\t"+str(num_keys))
			for i in range(0,min(100,num_keys)):
				dna=sorted_lens[i][0]
				dna_count=sorted_lens[i][1]
				numerator=float(dna_count)
				if(denominator!=0):
					ratio=numerator/float(denominator)
				else:
					ratio=0.0
				writer.write("\t"+str(ratio))
			writer.write("\n")
			writer.close()			
			file_num+=1





def getCDR3MinMax(dir_base):
	max_cdr3=None
	min_cdr3=None
	cdr3_len_set=set()
	for root, dirs, files in os.walk(dir_base):
		for items in fnmatch.filter(files, "*out.tsv"):
			full_tsv_path=root+"/"+items
			#print "Analyzing file ",full_tsv_path," for CDR3 min/max length...."
			bs=getBatchIDAndSampleIDFromPath(full_tsv_path)
			#print "The batch and sample :",bs
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
						cdr3_na_len=pieces[19]
						if(cdr3_na_len!="None" and cdr3_na_len!=""):
							cdr3_na_len=int(cdr3_na_len)
							cdr3_len_set.add(cdr3_na_len)
				line_num+=1
			reader.close()
	max_cdr3=max(cdr3_len_set)
	min_cdr3=min(cdr3_len_set)
	return [min_cdr3,max_cdr3]



if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Print DIOGENIX formatted rep_char CDR3 statistics')
	parser.add_argument('out_hist',type=str,nargs=1,help="output file path for TSV of sample/CDR3 length table")
	parser.add_argument('out_diva',type=str,nargs=1,help="output file path for TSV of CDR3 length diversity analysis")
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
		else:
			print "using dir_base ",dir_base
			cdr3_min_max=getCDR3MinMax(dir_base)
			min_len=cdr3_min_max[0]
			max_len=cdr3_min_max[1]
			print "min/max are :",cdr3_min_max
			print "Writing histogram to",out_hist
			cdr3LenStats(dir_base,out_hist,min_len,max_len)
			cdr3DiverStats(dir_base,out_diva)





