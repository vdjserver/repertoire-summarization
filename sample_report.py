#!/usr/bin/env python

from utils import extractAsItemOrFirstFromList,canReadAccess
from segment_utils import deAllelifyName
import argparse
import re

def getRegion(note):
	note_re=re.compile(r'^[A-Z]+(\d+)[A-Z]+$',re.IGNORECASE)
	sr=re.search(note_re,note)
	if(sr):
		num=int(sr.group(1))
		if(num>=31 and num<=35):
			return "CDR1"
		elif(num>=36 and num<=51):
			return "FR2"
		elif(num>=52 and num<=56):
			return "CDR2"
		elif(num>=57 and num<=92):
			return "FR3"
		else:
			raise Exception("bad number "+note)			
	else:
		raise Exception("bad number "+note)


def isCDR(s):
	if(s.startswith("CDR")):
		return True
	else:
		return False

def isFR(s):
	if(s.startswith("FR")):
		return True
	else:
		return False

def isNeither(s):
	if(not(isFR(s)) and not(isCDR(s))):
		return True
	else:
		return False



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
			if(pieces[0]=="AGS_SCORE"):
				AGS=float(pieces[1])
				AGS_RM=int(pieces[3])
				TOT_RM=float(pieces[5])
			elif(pieces[0]=="AGS5_SCRE"):
				AGS5=float(pieces[1])
				AGS5_RM=float(pieces[3])
			elif(pieces[0]=="NMO_SCORE"):
				NMO=float(pieces[1])
				NMO_Q_RM=int(pieces[3])
				NMO_NUC_TOT=int(pieces[5])
	reader.close()
	ret_map=dict()
	ret_map['total_reads']=total_reads
	ret_map['AGS']=AGS
	ret_map['AGS_RM']=AGS_RM
	ret_map['TOT_RM']=TOT_RM
	ret_map['AGS5']=AGS5
	ret_map['AGS5_RM']=AGS5_RM
	ret_map['NMO']=NMO
	ret_map['NMO_Q_RM']=NMO_Q_RM
	ret_map['NMU_NUC_TOT']=NMO_NUC_TOT
	return ret_map



def printStats(path,logPath):
	reader=open(path,'r')
	stat_counts=dict()
	#stat counts
	stat_list=list()
	stat_list.append("NoVHit")
	stat_list.append("Not an IGHV4 hit")
	stat_list.append("Had Indels!")
	stat_list.append("Found a stop codon")
	stat_list.append("VJ Out of Frame")
	stat_list.append("Incomplete regions")
	stat_list.append("OK")
	for stat in stat_list:
		stat_counts[stat]=0
	#V4 counts
	vList=list()
	vList.append("IGHV4-4")
	vList.append("IGHV4-28")
	vList.append("IGHV4-30-1")
	vList.append("IGHV4-30-2")
	vList.append("IGHV4-30-4")
	vList.append("IGHV4-31")
	vList.append("IGHV4-34")
	vList.append("IGHV4-38-2")
	vList.append("IGHV4-39")
	vList.append("IGHV4-55")
	vList.append("IGHV4-59")
	vList.append("IGHV4-61")
	vList.append("IGHV4-80")
	vList.append("IGHV4/OR15-8")
	vCounts=dict()
	for v in vList:
		vCounts[v]=0
	#J counts
	jList=list()
	jList.append("IGHJ1")
	jList.append("IGHJ2")
	jList.append("IGHJ3")
	jList.append("IGHJ4")
	jList.append("IGHJ5")
	jList.append("IGHJ6")
	jCounts=dict()
	for j in jList:
		jCounts[j]=0
	base_bp_count=195
	base_CDR_nuc_len=63
	line_num=1
	num_silent_mut=0
	tot_seq_len_by_bp=0
	tot_seq_len_by_codon=0
	tot_CDR_len_by_nuc=0
	tot_CDR_len_by_codon=0
	tot_RM_in_CDR=0
	tot_RM_in_FR=0
	tot_SM_in_CDR=0
	tot_SM_in_FR=0
	tot_bp_mut_in_CDR=0
	tot_bp_mut_in_FR=0
	for line in reader:
		#print line
		if(line_num>1):
			pieces=line.split('\t')
			status=pieces[183]
			stat_counts[status]+=1
			if(status=="OK"):
				#proceed
				vHit=pieces[2]
				vHitNA=deAllelifyName(vHit)
				if(vHitNA in vList):
					vCounts[vHitNA]+=1
				jHit=pieces[3]
				jHitNA=deAllelifyName(jHit)
				if(jHitNA in jList):
					jCounts[jHitNA]+=1
				#cdr3_len=pieces[21]
				#if(not(cdr3_len=="None")):
				#	cdr3_len=int(cdr3_len)
				#print "cdr3 len for ",pieces[1]," is ",cdr3_len
				rms=eval(pieces[185])
				sms=eval(pieces[187])
				bprms=eval(pieces[184])
				bpsms=eval(pieces[186])
				for silent_mutation in sms:
					num_silent_mut+=1
					sm_reg=getRegion(silent_mutation)
					if(isCDR(sm_reg)):
						tot_SM_in_CDR+=1
					elif(isFR(sm_reg)):
						tot_SM_in_FR+=1
					else:
						raise Exception("Error, unplacable silent mutation ",rm," for read="+pieces[1])
				for rm in rms:
					rm_reg=getRegion(rm)
					if(isCDR(rm_reg)):
						tot_RM_in_CDR+=1
					elif(isFR(rm_reg)):
						tot_RM_in_FR+=1
					else:
						raise Exception("Error, unplacable mutation ",rm," for read="+pieces[1])
				codon_change_set=set(bprms)
				codon_change_set=codon_change_set.union(set(bpsms))
				for cc in codon_change_set:
					codon_from=cc[0:3]
					codon_to=cc[len(cc)-3:]
					reg=getRegion(cc)
					for bpi in range(len(codon_from)):
						if(codon_from[bpi]!=codon_to[bpi]):
							#mutation detected!
							if(isCDR(reg)):
								tot_bp_mut_in_CDR+=1
							elif(isFR(reg)):
								tot_bp_mut_in_FR+=1
				#print "The codon change set for ",pieces[1]," is ",codon_change_set,"\n\n\n"
				#print "THE CDR1 LENGTH FOR "+pieces[1]+" is "+CDR1_len
				CDR1_len=int(pieces[110])
				if(CDR1_len in [15,18,21]):
					num_extra_bp=CDR1_len-15
					num_extra_aa=num_extra_bp/3
					tot_seq_len_by_bp+=base_bp_count+num_extra_bp
					tot_seq_len_by_codon+=(3*(base_bp_count+num_extra_bp))
					tot_CDR_len_by_nuc+=base_CDR_nuc_len+num_extra_bp
					tot_CDR_len_by_codon+=3*(base_CDR_nuc_len)+num_extra_aa
				else:
					raise Exception("Error on CDR1 length in file "+path+" for read="+pieces[1])
			else:
				pass
		line_num+=1
	reader.close()
	log_data=readLog(logPath)
	out_pieces=list()
	out_pieces_header=list()
	out_pieces_header.append("BATCH ID")
	out_pieces.append("BID")
	out_pieces_header.append("SAMP ID")
	out_pieces.append("SID")
	out_pieces_header.append("TOTAL READS")
	out_pieces.append(log_data['total_reads'])
	for stat in stat_list:
		#print stat_counts[stat]+"\t"
		out_pieces_header.append("FILTER "+stat+" COUNT")
		out_pieces.append(stat_counts[stat])
	for v in vList:
		out_pieces_header.append("V4 GENE "+v+" COUNT")
		out_pieces.append(vCounts[v])
	for j in jList:
		out_pieces_header.append("J GENE "+j+" COUNT")
		out_pieces.append(jCounts[j])
	#AGS5 and AGS5RM
	out_pieces_header.append("AGS5")
	out_pieces.append(log_data['AGS5'])
	out_pieces_header.append("AGS5 RM")
	out_pieces.append(log_data['AGS5_RM'])
	out_pieces_header.append("AGS6")
	out_pieces.append(log_data['AGS'])
	out_pieces_header.append("AGS6 RM")
	out_pieces.append(log_data['AGS_RM'])

	#TOT RM SM
	out_pieces_header.append("TOT RM")
	out_pieces.append(log_data['TOT_RM'])
	out_pieces_header.append("TOT SM")
	out_pieces.append(num_silent_mut)

	#items dependent on length
	out_pieces_header.append("Total seq length by codon")
	out_pieces.append(tot_seq_len_by_codon)
	
	



	#stringify
	for op in range(len(out_pieces)):
		out_pieces[op]=str(out_pieces[op])
	joiner="\t"
	if(len(out_pieces)!=len(out_pieces_header)):
		print "header len=",len(out_pieces_header)
		print "value  len=",len(out_pieces)
		raise Exception("Error, header, value length mismatch!")
	print joiner.join(out_pieces_header)
	print joiner.join(out_pieces)
	


if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Print DIOGENIX formatted rep_char run statistics')
	parser.add_argument('tsv_in',type=str,nargs=1,help="path to a rep_char TSV output file")
	parser.add_argument('log_in',type=str,nargs=1,help="path to the corresponding rep_char output log file")
	args=parser.parse_args()
	if(not(args)):
		parser.print_help()
	input_file_path=extractAsItemOrFirstFromList(args.tsv_in)
	log_file_path=extractAsItemOrFirstFromList(args.log_in)
	if(not(canReadAccess(input_file_path))):
		parser.print_help()
	if(not(canReadAccess(log_file_path))):
		parser.print_help()
	printStats(input_file_path,log_file_path)









