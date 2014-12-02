#!/usr/bin/env python

from utils import extractAsItemOrFirstFromList,canReadAccess,glob_walk
from segment_utils import deAllelifyName
from sample_cdr3 import getBatchIDAndSampleIDFromPath
import argparse
import re
import os
from ags_mgr import ags_manager

def getLogPathGivenTSVPath(tsv_path):
	tsv_path_len=len(tsv_path)
	tsv_log_path=tsv_path[:(tsv_path_len-3)]+"log"
	#print "given ",tsv_path," returning log ",tsv_log_path
	return tsv_log_path
	


def getRegion(note):
	note_re=re.compile(r'^[A-Z\*]+(\d+)[A-Z\*]+$',re.IGNORECASE)
	sr=re.search(note_re,note)
	if(sr):
		num=int(sr.group(1))
		if(num>=31 and num<=35):
			return "CDR1"
		elif(num>=36 and num<=49):
			return "FR2"
		elif(num>=50 and num<=65):
			return "CDR2"
		elif(num>=66 and num<=92):
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
	if(not(os.path.isfile(logPath))):
		logPath="/dev/null"
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
				if(pieces[1]!="None"):
					AGS=float(pieces[1])
				else:
					AGS=0.0
				AGS_RM=int(pieces[3])
				TOT_RM=float(pieces[5])
			elif(pieces[0]=="AGS5_SCRE"):
				if(pieces[1]!="None"):
					AGS5=float(pieces[1])
				else:
					AGS5=0.0
				AGS5_RM=float(pieces[3])
			elif(pieces[0]=="NMO_SCORE"):
				if(pieces[1]!="None"):
					NMO=float(pieces[1])
				else:
					NMO=0.0
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
	ret_map['NMO_NUC_TOT']=NMO_NUC_TOT
	return ret_map



def printStats(path,logPath,bid,sid,minCountNum,returnInsteadOfPrint=False):
	reader=open(path,'r')
	stat_counts=dict()
	#stat counts
	stat_list=list()
	stat_list.append("MinCount_"+str(minCountNum)+"_Fail")
	stat_list.append("NoVHit")
	stat_list.append("Not an IGHV4 hit")
	stat_list.append("Had Indels!")
	stat_list.append("Found a stop codon")
	stat_list.append("VJ Out of Frame")
	stat_list.append("Incomplete regions")
	stat_list.append("Fail Homology Filter 0.85")
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
	base_FR_AA_len=44
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
	tot_bp_mut_all_reg=0
	tot_len_AA_in_FR=0
	tot_len_bp_in_FR=0
	num_seqs_with_equal_or_greater_rm=dict()
	rm_gq_max=6
	for y in range(0,rm_gq_max+1):
		num_seqs_with_equal_or_greater_rm[y]=0
	num_seqs_with_equal_or_greater_nuc_muts=dict()
	num_mut_max=6
	for y in range(0,num_mut_max+1):
		num_seqs_with_equal_or_greater_nuc_muts[y]=0
	homology_over_v_thresh=85.00
	num_homology_over_v_less_than_thresh=0
	mySampleAGSMGR=ags_manager(str(bid+"."+sid))

	#per GENE AGS managers
	per_gene_ags_init_list=list()
	per_gene_ags_init_list.append("IGHJ1")
	per_gene_ags_init_list.append("IGHJ2")
	per_gene_ags_init_list.append("IGHJ3")
	per_gene_ags_init_list.append("IGHJ4")
	per_gene_ags_init_list.append("IGHJ5")
	per_gene_ags_init_list.append("IGHJ6")
	per_gene_ags_init_list.append("IGHV4-4")
	per_gene_ags_init_list.append("IGHV4-61")
	per_gene_ags_init_list.append("IGHV4-38")
	per_gene_ags_init_list.append("IGHV4-39")
	per_gene_ags_init_list.append("IGHV4-34")
	per_gene_ags_init_list.append("IGHV4-59")
	per_gene_ags_init_list.append("IGHV4-30")
	per_gene_ags_init_list.append("IGHV4-31")
	ags_by_gene=dict()
	for gene in per_gene_ags_init_list:
		ags_by_gene[gene]=ags_manager(gene)


	#per RM+ mgr
	ags_by_rmp=dict()
	ags_by_rmp[2]=ags_manager(str(2))
	ags_by_rmp[3]=ags_manager(str(3))

	#counts for weighted filtering
	weightedNoIGHV4Hit=0
	weightedNoVHit=0
	weightedSumAll=0
	weightedSumOKOnly=0


	for line in reader:
		#print line
		if(line_num>1):
			pieces=line.split('\t')
			status=pieces[183]
			count_val=int(pieces[len(pieces)-1])
			#print "for ",pieces[1]," count is ",str(count_val)
			#llsys.exit(0)
			if(count_val>=minCountNum):
				stat_counts[status]+=1
				weightedSumAll+=count_val
				if(status=="NoVHit"):
					weightedNoVHit+=count_val
				if(status=="Not an IGHV4 hit"):
					weightedNoIGHV4Hit+=count_val
				if(status=="OK"):
					#proceed
					weightedSumOKOnly+=count_val
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
					rms=eval(pieces[185]) #replacement mutations amino
					sms=eval(pieces[187]) #silent mutations amino
					bprms=eval(pieces[184]) #replacment mutation codons
					bpsms=eval(pieces[186]) #silent mutation codons
					homology_over_v=float(pieces[5])
					if(homology_over_v<homology_over_v_thresh):
						num_homology_over_v_less_than_thresh+=1
					for silent_mutation in sms:
						num_silent_mut+=1
						sm_reg=getRegion(silent_mutation)
						if(isCDR(sm_reg)):
							tot_SM_in_CDR+=1
							#print "SM",silent_mutation," detected in CDR for ",pieces[1]
						elif(isFR(sm_reg)):
							#print "SM",silent_mutation," detected in FR for ",pieces[1]
							tot_SM_in_FR+=1
						else:
							raise Exception("Error, unplacable silent mutation ",rm," for read="+pieces[1])
					num_rms_in_this_row=0
					for rm in rms:
						mySampleAGSMGR.receive_numbered_mut(rm)
						for gene in per_gene_ags_init_list:
							if(jHit.startswith(gene)):
								ags_by_gene[gene].receive_numbered_mut(rm)
							if(vHit.startswith(gene)):
								ags_by_gene[gene].receive_numbered_mut(rm)
						if(len(rms)>=2):
							ags_by_rmp[2].receive_numbered_mut(rm)
							if(len(rms)>=3):
								ags_by_rmp[3].receive_numbered_mut(rm)
						rm_reg=getRegion(rm)
						num_rms_in_this_row+=1
						if(isCDR(rm_reg)):
							#print "RM",rm," detected in CDR for ",pieces[1]
							tot_RM_in_CDR+=1
						elif(isFR(rm_reg)):
							#print "RM",rm," detected in FR for ",pieces[1]
							tot_RM_in_FR+=1
						else:
							raise Exception("Error, unplacable mutation ",rm," for read="+pieces[1])
					num_seqs_with_equal_or_greater_rm[min(num_rms_in_this_row,rm_gq_max)]+=1
					#here MERGE the two sets (both silent and non-silent by UNION)
					codon_change_set=set(bprms)
					codon_change_set=codon_change_set.union(set(bpsms))
					num_nuc_muts_this_row=0
					for cc in codon_change_set:
						codon_from=cc[0:3]
						codon_to=cc[len(cc)-3:]
						reg=getRegion(cc)
						for bpi in range(len(codon_from)):
							if(codon_from[bpi]!=codon_to[bpi]):
								#mutation detected!
								num_nuc_muts_this_row+=1
								tot_bp_mut_all_reg+=1
								if(isCDR(reg)):
									tot_bp_mut_in_CDR+=1
								elif(isFR(reg)):
									tot_bp_mut_in_FR+=1
					num_seqs_with_equal_or_greater_nuc_muts[min(num_mut_max,num_nuc_muts_this_row)]+=1
					#print "The codon change set for ",pieces[1]," is ",codon_change_set,"\n\n\n"
					#print "THE CDR1 LENGTH FOR "+pieces[1]+" is "+CDR1_len
					CDR1_len=int(pieces[110])
					if(CDR1_len in [15,18,21]):
						num_extra_bp=CDR1_len-15
						num_extra_aa=num_extra_bp/3
						tot_seq_len_by_bp+=base_bp_count+num_extra_bp
						tot_seq_len_by_codon+=((base_bp_count+num_extra_bp)/3)
						tot_CDR_len_by_nuc+=base_CDR_nuc_len+num_extra_bp
						tot_CDR_len_by_codon+=(((base_CDR_nuc_len)/3)+num_extra_aa)
					else:
						raise Exception("Error on CDR1 length in file "+path+" for read="+pieces[1])
					tot_len_AA_in_FR+=base_FR_AA_len
					tot_len_bp_in_FR+=(base_FR_AA_len*3)
				else:
					pass
			else:
				#count_val<min_count_num
				stat_counts["MinCount_"+str(minCountNum)+"_Fail"]+=1
				pass
		line_num+=1
	reader.close()
	log_data=readLog(logPath)
	out_pieces=list()
	out_pieces_header=list()
	out_pieces_header.append("BATCH ID")
	out_pieces.append(bid)
	out_pieces_header.append("SAMP ID")
	out_pieces.append(sid)
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

	#PER GENE AGS data
	for gene in per_gene_ags_init_list:
		types=[5,6]
		for ags_type in types:
			if(gene.startswith("IGHV")):
				gene_ags_info=ags_by_gene[gene].compute_ags("AGS"+str(ags_type))
				#AGS score
				out_pieces_header.append(gene+"_AGS"+str(ags_type))
				out_pieces.append(gene_ags_info[0])
				#AGS AGS TOT	
				out_pieces_header.append(gene+"_AGS"+str(ags_type)+"_TOT")
				out_pieces.append(gene_ags_info[1])
				#TOT	
				out_pieces_header.append(gene+"_AGS"+str(ags_type)+"_SAMPTOT")
				out_pieces.append(gene_ags_info[2])
			else:
				gene_ags_info=ags_by_gene[gene].compute_ags("AGS"+str(ags_type))
				#AGS score
				out_pieces_header.append(gene+"_AGS"+str(ags_type))
				out_pieces.append(gene_ags_info[0])
				#AGS AGS TOT	
				out_pieces_header.append(gene+"_AGS"+str(ags_type)+"_AGSTOT")
				out_pieces.append(gene_ags_info[1])
				#TOT	
				out_pieces_header.append(gene+"_AGS"+str(ags_type)+"_TOT")
				out_pieces.append(gene_ags_info[2])				


	#weighted filtering counts
	out_pieces_header.append("Filter Weighted Not an IGHV4 hit")
	out_pieces.append(weightedNoIGHV4Hit)
	out_pieces_header.append("Filter Weighted NoVHit")
	out_pieces.append(weightedNoVHit)
	out_pieces_header.append("Weighed Sum ALL")
	out_pieces.append(weightedSumAll)
	out_pieces_header.append("Weighted Sum OK Only")
	out_pieces.append(weightedSumOKOnly)

	#REPLACE LOG AGS/RM data with AGS_MGR dynamically computed data
	ags_5_info=mySampleAGSMGR.compute_ags("AGS5")#[ags5_score,sampAGSTot,sampTot]
	ags_6_info=mySampleAGSMGR.compute_ags("AGS6")#[ags6_score,sampAGSTot,sampTot]
	#print "ags_5_info : ",ags_5_info
	#print "ags_6_info : ",ags_6_info
	log_data['AGS5']=ags_5_info[0]
	log_data['AGS5_RM']=ags_5_info[1]
	log_data['AGS']=ags_6_info[0]
	log_data['AGS_RM']=ags_6_info[1]
	log_data['TOT_RM']=ags_6_info[2]

	#AGS5 and AGS5RM
	out_pieces_header.append("AGS5")
	out_pieces.append(log_data['AGS5'])
	out_pieces_header.append("AGS5 RM")
	out_pieces.append(log_data['AGS5_RM'])
	#do AGS5 RM2+,RM3+ scores
	for rmp in range(2,4):
		agss_rm_info=ags_by_rmp[rmp].compute_ags("AGS5")
		out_pieces_header.append("RM"+str(rmp)+"+ AGS5")
		out_pieces.append(agss_rm_info[0])
		out_pieces_header.append("RM"+str(rmp)+"+ AGS5 RM")
		out_pieces.append(agss_rm_info[1])
		out_pieces_header.append("RM"+str(rmp)+"+ AGS5 TOT RM")
		out_pieces.append(agss_rm_info[2])


	out_pieces_header.append("AGS6")
	out_pieces.append(log_data['AGS'])
	out_pieces_header.append("AGS6 RM")
	out_pieces.append(log_data['AGS_RM'])
	#do AGS5 RM2+,RM3+ scores
	for rmp in range(2,4):
		agss_rm_info=ags_by_rmp[rmp].compute_ags("AGS6")
		out_pieces_header.append("RM"+str(rmp)+"+ AGS6")
		out_pieces.append(agss_rm_info[0])
		out_pieces_header.append("RM"+str(rmp)+"+ AGS6 RM")
		out_pieces.append(agss_rm_info[1])
		out_pieces_header.append("RM"+str(rmp)+"+ AGS6 TOT RM")
		out_pieces.append(agss_rm_info[2])


	#TOT RM SM
	out_pieces_header.append("TOT RM")
	out_pieces.append(log_data['TOT_RM'])
	out_pieces_header.append("TOT SM")
	out_pieces.append(num_silent_mut)

	#items dependent on length
	out_pieces_header.append("Total seq length by codon")
	out_pieces.append(tot_seq_len_by_codon)
	out_pieces_header.append("Tot RM in CDR")
	out_pieces.append(tot_RM_in_CDR)
	out_pieces_header.append("Tot SM in CDR")
	out_pieces.append(tot_SM_in_CDR)
	out_pieces_header.append("Tot CDR Length by Codon")
	out_pieces.append(tot_CDR_len_by_codon)
	out_pieces_header.append("Tot RM in FR")
	out_pieces.append(tot_RM_in_FR)	
	out_pieces_header.append("Tot SM in FR")
	out_pieces.append(tot_SM_in_FR)		
	out_pieces_header.append("Tot FR Length by Codon")
	out_pieces.append(tot_len_AA_in_FR)
	out_pieces_header.append("Tot nucleotide mutations")
	out_pieces.append(tot_bp_mut_all_reg)
	out_pieces_header.append("Tot seq. len by nucleotide")
	out_pieces.append(tot_seq_len_by_bp)
	out_pieces_header.append("Tot nucleotide CDR mutations")
	out_pieces.append(tot_bp_mut_in_CDR)
	out_pieces_header.append("Tot seq. len by CDR nucleotide")
	out_pieces.append(tot_CDR_len_by_nuc)
	out_pieces_header.append("Tot. nucleotide FR mutations")
	out_pieces.append(tot_bp_mut_in_FR)
	out_pieces_header.append("Tot FR length by nucleotide")
	out_pieces.append(tot_len_bp_in_FR)

	
	#NMO
	out_pieces_header.append("NMO Score (no norm by 100)")
	out_pieces.append(log_data['NMO']/100.0)
	out_pieces_header.append("NMO Nuc Muts (sum nuc mutations in NMO codons)")
	out_pieces.append(log_data['NMO_NUC_TOT'])
	out_pieces_header.append("# Seqs with RM>=1")
	out_pieces.append(log_data['NMO_Q_RM'])


	#RM FREQUENCY SPECTRA
	for rm_count in range(0,rm_gq_max+1):
		if(rm_count==rm_gq_max):
			out_pieces_header.append("NUM_SEQS_WITH_RM_GEQ_"+str(rm_count))
		else:
			out_pieces_header.append("NUM_SEQS_WITH_RM_"+str(rm_count))
		out_pieces.append(str(num_seqs_with_equal_or_greater_rm[rm_count]))

	#HOMOLOGY THRESHOLD
	out_pieces_header.append("NUM_SEQS_WITH_V_HOMOLOGY_LT_"+str(homology_over_v_thresh))
	out_pieces.append(str(num_homology_over_v_less_than_thresh))


	#MUT FREQ SPECTRA
	for nuc_mut_count in range(0,num_mut_max+1):
		if(nuc_mut_count==num_mut_max):
			out_pieces_header.append("NUM_SEQS_WITH_NUC_MUT_GEQ_"+str(nuc_mut_count))		
		else:
			out_pieces_header.append("NUM_SEQS_WITH_NUC_MUT_"+str(nuc_mut_count))
		out_pieces.append(str(num_seqs_with_equal_or_greater_nuc_muts[nuc_mut_count]))


	#stringify
	for op in range(len(out_pieces)):
		out_pieces[op]=str(out_pieces[op])
	joiner="\t"
	if(len(out_pieces)!=len(out_pieces_header)):
		print "header len=",len(out_pieces_header)
		print "value  len=",len(out_pieces)
		raise Exception("Error, header, value length mismatch!")
	if(not(returnInsteadOfPrint)):
		print joiner.join(out_pieces_header)
		print joiner.join(out_pieces)
	else:
		header=joiner.join(out_pieces_header)
		values=joiner.join(out_pieces)
		return [header,values]


if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Print DIOGENIX formatted rep_char run statistics')
	#parser.add_argument('tsv_in',type=str,nargs=1,help="path to a rep_char TSV output file")
	#parser.add_argument('log_in',type=str,nargs=1,help="path to the corresponding rep_char output log file")
	parser.add_argument('tsv_base',type=str,nargs=1,help="path to a directory itself holding BATCHES, each BATCH directory holding samples")
	parser.add_argument('-min_count_num',type=int,default=1,nargs=1,help="required minimum count value that triggers row into aggregtation (default 1)")
	args=parser.parse_args()
	if(not(args)):
		parser.print_help()
		sys.exit(0)
	#input_file_path=extractAsItemOrFirstFromList(args.tsv_in)
	#log_file_path=extractAsItemOrFirstFromList(args.log_in)
	input_base=extractAsItemOrFirstFromList(args.tsv_base)
	min_count_num=extractAsItemOrFirstFromList(args.min_count_num)
	input_glob=input_base+"/*/*.rc_out.tsv"
	tsv_num=1
	#print "using min_count_num=",min_count_num
	#sys.exit(0)
	for TSV in glob_walk(input_glob):
		LOG=getLogPathGivenTSVPath(TSV)
		bid_sid=getBatchIDAndSampleIDFromPath(TSV)
		#print "FROM TSV=",TSV
		#print "\tLOG=",LOG
		#print "\tBID,SID=",bid_sid
		BID=bid_sid[0]
		SID=bid_sid[1]
		returnInsteadOfPrint=True
		report_lines=printStats(TSV,LOG,BID,SID,min_count_num,returnInsteadOfPrint)
		if(tsv_num==1):
			print report_lines[0]
		print report_lines[1]
		tsv_num+=1










