#!/usr/bin/env python

from utils import extractAsItemOrFirstFromList,canReadAccess
import argparse


def readLog(logPath):
	reader=open(logPath,'r')
	total_re=re.compile(r'total\sof\s(\d+)\sreads')
	total_reads=0
	AGS=0
	AGS_RM=0
	TOT_RM=0
	AGS5=0
	AGS5_RM=0
	NMO=0
	NMO_Q_RM=0
	NMO_NUC_TOT=0
	for line in reader:
		#Processed total of 53 reads!
		#Wrote CDR3 lengths histogram to  human.IG.fna.igblast.kabat.out.cdr3_hist.tsv
		#AGS_SCORE	11.2532123961	AGS6_TOT_RM	29	TOT_RM	147
		#AGS5_SCRE	12.2751322751	AGS5_TOT_RM	28	TOT_RM	147
		#NMO_SCORE	109.523809524	NMO_QUERIES_RM	21	NMO_SAMP_NUC_TOT	23
		re_res=re.search(total_re,line)
		if(re_res):
			total_reads=int(re_res.group(1))
		pieces=line.split('\t')
		if(len(pieces)==6):
			if(pieces[0]=="AGS_SCORE")
				AGS=float(pieces[1])
				AGS_RM=int(pieces[3])
				TOT_RM=float(pieces[5])
			elif(pieces[0]=="AGS5_SCRE"):
				AGS5_SCRE=float(pieces[1])
			elif(pieces[0]=="NMO_SCORE"):
				
		
	reader.close()


def printStats(path):
	reader=open(path,'r')
	stat_counts=dict()
	#stat counts
	stat_counts["Found a stop codon"]=0
	stat_counts["Had Indels!"]=0
	stat_counts["Incomplete regions"]=0
	stat_counts["NoVHit"]=0
	stat_counts["OK"]=0
	stat_counts["VJ Out of Frame"]=0
	stat_counts["Not an IGHV4 hit"]
	#V4 counts
	vCounts=dict()
	vCounts["IGHV4-4"]=0
	vCounts["IGHV4-28"]=0
	vCounts["IGHV4-30-1"]=0
	vCounts["IGHV4-30-2"]=0
	vCounts["IGHV4-30-4"]=0
	vCounts["IGHV4-31"]=0
	vCounts["IGHV4-34"]=0
	vCounts["IGHV4-38-2"]=0
	vCounts["IGHV4-39"]=0
	vCounts["IGHV4-55"]=0
	vCounts["IGHV4-59"]=0
	vCounts["IGHV4-61"]=0
	vCounts["IGHV4-80"]=0
	vCounts["IGHV4/OR15-8"]=0
	#J counts
	jCounts=dict()
	jCounts["IGHJ1"]=0
	jCounts["IGHJ2"]=0
	jCounts["IGHJ3"]=0
	jCounts["IGHJ4"]=0
	jCounts["IGHJ5"]=0
	jCounts["IGHJ6"]=0


	for line in reader:
		#print line
		pieces=line.split('\t')
		status=pieces[184]
		stat_counts[status]+=1
		if(status=="OK"):
			#proceed
		else:
			pass
		
	
	reader.close()


if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Print DIOGENIX formatted rep_char run statistics')
	parser.add_argument('tsv_in',type=str,nargs=1,help="path to a rep_char TSV output file")
	parser.add_argument('log_in',type=str,nargs=1,help="path to the corresponding rep_char output log file")
	args=parser.parse_args()
	if(not(args)):
		parser.print_help()
	input_file_path=extractAsItemOrFirstFromList(args.tsv_in)
	if(not(canReadAccess(input_file_path))):
		parser.print_help()
	printStats(input_file_path)









