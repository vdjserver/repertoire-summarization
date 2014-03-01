#!/usr/bin/env python

import os
import re
import glob
from os.path import basename
from segment_utils import getFastaListOfDescs,getQueryIndexGivenSubjectIndexAndAlignment,getAdjustedCDR3StartFromRefDirSetAllele,getADJCDR3EndFromJAllele,looksLikeAlleleStr
from imgt_utils import imgt_db
from utils import read_fasta_file_into_map,biopythonTranslate,rev_comp_dna,printMap
from igblast_utils import printNiceAlignment,buildAlignmentWholeSeqs,buildAlignmentWholeSeqsDirect
import vdjml
from vdjml_utils import getTopVDJItems,getHitInfo



#test

#define inputs for run
vfasta="/home/esalina2/round1_imgt/human_IG_V.fna"
dfasta="/home/esalina2/round1_imgt/human_IG_D.fna"
jfasta="/home/esalina2/round1_imgt/human_IG_J.fna"
fasta_file_list=[vfasta,dfasta,jfasta]
blast_out="/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna.igblastn.imgt.out"

#define inputs for database
base_dir="/home/data/DATABASE/01_22_2014/"
organism="human"

#create DB OBJ
imgtdb_obj=imgt_db(base_dir)
#getAdjustedCDR3StartFromRefDirSetAllele(allele,imgtdb_obj,organism="human"):

vlist=getFastaListOfDescs(fasta_file_list[0])
dlist=getFastaListOfDescs(fasta_file_list[1])
jlist=getFastaListOfDescs(fasta_file_list[2])

blast_data=open(blast_out,'r')
line_num=0
counts_map=dict()
capture_flag=False
is_heavy=True
VHitLine=None
JHitLine=None
currentV=None
currentJ=None
currentD=None

trans=0
tot=0
nona=0
max_len=0
cdr3_hist=dict()
query_seq_map=read_fasta_file_into_map("/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna")



def get_domain_modes():
	domain_list=["imgt","kabat"]
	#domain_list=["kabat"]
	return domain_list


rsmap=dict()
remap=dict()
domain_list=get_domain_modes()
for d in domain_list:
	rsmap[d]=dict()
	remap[d]=dict()
	cdr3_hist[d]=dict()
	for i in range(0,1000):
		cdr3_hist[d][i]=0
v_alleles_review=dict()




vqmap=dict()


files=glob.glob("/home/esalina2/Downloads/HSEQS/IMGT_HighV-QUEST_individual_files_folder/*")
files=[]
for f in files:
	reader=open(f,'r')
	num=0
	fbase=basename(f)
	fbpieces=fbase.split("_")
	readname=str(fbpieces[0].strip())
	vqmap[readname]=dict()
	for line in reader:
		segs=['V','D','J']
		for seg in segs:
			if(not(seg in vqmap[readname])):
				vqmap[readname][seg]=None
			if(line.startswith(seg+"-GENE")):
				pieces=line.split(";")
				if(len(pieces)>=2):
					info=pieces[1]
					segre=re.compile(r'homsap\s+(\S+)\s+',re.IGNORECASE)
					mr=re.search(segre,info)
					if(mr):
						seg_res=mr.group(1)
						#print "readname is ",readname
						vqmap[readname][seg]=seg_res
						print "\n\n\n*************************"
						print "readname="+readname
						print "seg="+seg
						print "storage="+seg_res
						print "*******************************\n\n\n"
							
					else:
						print "FAIL re MATCH ON ",readname," : ",info," for seg=",seg
		if(num>=40):
			break;
		num+=1
	reader.close()


def getNumDashInStr(s):
	orig_len=len(s)
	s_no_dash=re.sub(r'\-','',s)
	nodash_len=len(s_no_dash)
	return orig_len-nodash_len


def repStr(s,n):
	if(n<=0):
		return ""
	else:
		return s+repStr(s,n-1)



def annotatedCDR3(s_aln,q_aln,s_c,q_c,s_s,q_s):
	#get CDR3 starts relative to the alignment starts
	aln_rel_s=s_c-s_s
	aln_rel_q=q_c-q_s
	print "subject from is ",s_s
	print "subject cdr3 is ",s_c
	print "aln rel sub is ",aln_rel_s
	print "query from is ",q_s
	print "query CDR3 is ",q_c
	print "aln rel q is ",aln_rel_q
	#sys.exit(0)


	#find position including gaps in the alignment for the SUBJECT STAR
	#all this code is 1-based indices
	#except for actually doing string comparison and string-printing which is 0-based
	s_temp=0
	if(s_aln[0]=="-"):
		s_pos=0
	else:
		s_pos=1
	while(s_pos<=aln_rel_s):
		if(s_aln[s_temp]!="-"):
			s_pos+=1
		s_temp+=1
	s_actual=s_temp
	#0-based indices!
	q_temp=0
	if(q_aln[0]=="-"):
		q_pos=0
	else:
		q_pos=1
	while(q_pos<=aln_rel_q):
		if(q_aln[q_temp]!="-"):
			q_pos+=1
		q_temp+=1
	q_actual=q_temp 
	#0-based indices
	#if(q_actual!=s_actual):
	#	print "***WARNING***"
	annLines=["","","",""]
	annLines[0]=repStr(" ",s_actual)
	annLines[0]+="X"
	annLines[1]=s_aln
	annLines[2]=q_aln
	annLines[3]=repStr(" ",q_actual)
	annLines[3]+="X"
	annotated=""
	for i in range(len(annLines)):
		#print annLines[i]
		annotated+=str(i)+" : "+annLines[i]
		if(i<len(annLines)-1):
			annotated+="\n"
	return annotated
	

	


def getOtherMode(m):
	if(m=="imgt"):
		return "kabat"
	else:
		return "imgt"



def addAlignmentsPreCDR3(dataMap,alleleName,imgtdb_obj,organism,query_record):
	btop=dataMap['btop']
	bopyname=str(query_record.id)
	#print "the biopython name is ",bopyname
	query_seq=str(query_record.seq)
	if(dataMap['is_inverted']):
		#print "inversion is necessary...."
		query_seq=rev_comp_dna(query_seq)
	query_seq=query_seq[int(dataMap['q. start']-1):int(dataMap['q. end'])]
	subject_seq=imgtdb_obj.getRefDirSetFNASeqGivenOrgAndAllele(alleleName,organism)
	subject_seq=subject_seq[int(dataMap['s. start']-1):int(dataMap['s. end'])]
	qsAlgn=buildAlignmentWholeSeqs(btop,query_seq,subject_seq)
	dataMap['query seq']=qsAlgn[0]
	dataMap['subject seq']=qsAlgn[2]
	return dataMap

def CDR3LengthAnalysisVDMLOBJ(read_result_obj,meta,organism,imgtdb_obj,query_record):
	#for V and J require:
	# 1) allele names
	# 2) alignment from BTOP reconstruction both Q and S
	# 3) q. start and q. end and s. start and s. end
	# 4) read inversion flags
	#print "got into cdr3 hist wrapper"
	read_name=read_result_obj.id()
	#print "got read name ",read_name
	segment_combinations=read_result_obj.segment_combinations()
	#if(len(segment_combinations)
	#print "the length is ",len(segment_combinations)
	VDJMap=getTopVDJItems(read_result_obj,meta)
	vAllele=VDJMap['V']
	jAllele=VDJMap['J']
	empty_map=dict()
	empty_map['kabat']=(-1)
	empty_map['imgt']=(-1)
	if(vAllele!=None and jAllele!=None):
		vData=getHitInfo(read_result_obj,meta,vAllele)
		#print "Vdata IS "
		vData=addAlignmentsPreCDR3(vData,vAllele,imgtdb_obj,organism,query_record)
		jData=getHitInfo(read_result_obj,meta,jAllele)
		jData=addAlignmentsPreCDR3(jData,jAllele,imgtdb_obj,organism,query_record)
		#printMap(vData)
		#printMap(jData)
		cdr3_analysis_map=CDR3LengthAnalysis(vData,jData,vData['query id'],organism,imgtdb_obj)
		return cdr3_analysis_map
		#printMap(cdr3_analysis_map)
	else:
		#print "insufficient data!"
		return empty_map
	pass



class histoMapClass:

	#counter
	count_map=dict()
	modes=None

	#constructor
	def __init__(self,modes):
		self.count_map=dict()
		self.modes=modes
		for mode in modes:
			self.count_map[mode]=dict()

	#increment
	def inc(self,mode,val):
		val=int(val)
		if(not val in self.count_map[mode]):
			self.count_map[mode][val]=1
		else:
			self.count_map[mode][val]+=1

	#get min val over all vals
	def gminVal(self):
		min_val=float("inf")
		for mode in self.modes:
			for val in self.count_map[mode]:
				if(val<min_val):
					min_val=val
		return min_val

	#get max val over all vals
	def gmaxVal(self):
		max_val=float("-inf")
		for mode in self.modes:
			for val in self.count_map[mode]:
				if(val>max_val):
					max_val=val
		return max_val
							


				

	#print maps basic
	def printMaps(self):
		for mode in self.modes:
			print "SHOWING MAP FOR ",mode
			for val in self.count_map[mode]:
				print val,"\t",self.count_map[mode][val]

	#write a basic histogram to a file!
	def writeToFile(self,ofile):
		writer=open(ofile,'w')
		min_val=self.gminVal()
		max_val=self.gmaxVal()
		print "min and max are ",min_val,max_val
		round_num=0
		for v in range(min_val,max_val+1):
			if(round_num==0):
				writer.write("VALUE")
				for mode in self.modes:
					writer.write("\t"+mode)
				writer.write("\n")
			writer.write(str(v))
			for mode in self.modes:
				actual_val=0
				if(v in self.count_map[mode]):
					actual_val=self.count_map[mode][v]
				writer.write("\t"+str(actual_val))
			writer.write("\n")
			round_num+=1
		writer.close()




def CDR3LengthAnalysis(vMap,jMap,currentQueryName,organism,imgtdb_obj):
	currentQueryName=str(currentQueryName.strip())
	currentV=vMap['subject ids']
	currentJ=jMap['subject ids']
	cdr3_hist=dict()
	if(looksLikeAlleleStr(currentV) and looksLikeAlleleStr(currentJ)):
		#print "WE'RE IN BUSINESS!"
		domain_modes=["kabat","imgt"]
		for dm in domain_modes:
			cdr3_hist[dm]=(-1)
			#print "processing in ",dm
			if(not currentV in rsmap[dm]):
				#print currentV,"not in lookup for dm=",dm
				ref_cdr3_start=getAdjustedCDR3StartFromRefDirSetAllele(currentV,imgtdb_obj,organism,dm)
				rsmap[dm][currentV]=ref_cdr3_start
			else:
				ref_cdr3_start=rsmap[dm][currentV]
			if(not currentJ in remap[dm]):
				#print currentJ,"not in lookup for dm=",dm
				ref_cdr3_end=getADJCDR3EndFromJAllele(currentJ,imgtdb_obj,organism,dm)
				remap[dm][currentJ]=ref_cdr3_end
			else:
				ref_cdr3_end=remap[dm][currentJ]
			if(ref_cdr3_start!=(-1) and ref_cdr3_end!=(-1)):
				vq_aln=vMap['query seq']
				vs_aln=vMap['subject seq']
				vq_f=int(vMap['q. start'])
				vq_t=int(vMap['q. end'])
				vs_f=int(vMap['s. start'])
				vs_t=int(vMap['s. end'])
				jq_aln=jMap['query seq']
				js_aln=jMap['subject seq']
				jq_f=int(jMap['q. start'])
				jq_t=int(jMap['q. end'])
				js_f=int(jMap['s. start'])
				js_t=int(jMap['s. end'])
				qry_cdr3_start=getQueryIndexGivenSubjectIndexAndAlignment(vq_aln,vs_aln,vq_f,vq_t,vs_f,vs_t,ref_cdr3_start)
				qry_cdr3_end=getQueryIndexGivenSubjectIndexAndAlignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t,ref_cdr3_end)
				if(qry_cdr3_start!=(-1) and qry_cdr3_end!=(-1)):
					query_coding_seq=query_seq_map[currentQueryName]
					#if(vMap['query id'].find("reversed|"+currentQueryName)==(-1)):
					#	query_coding_seq=rev_comp_dna(query_coding_seq)
					#coding_seq=query_coding_seq[(qry_cdr3_start-1):(qry_cdr3_end-1)]
					cdr3_len=qry_cdr3_end-qry_cdr3_start+1
					#translation=biopythonTranslate(coding_seq)
					#print "the coding seq ("+dm+") is : ",coding_seq
					#print "The translation ("+dm+") is : ",translation
					#print "CDR3_LEN ("+dm+") ="+str(cdr3_len)
					cdr3_hist[dm]=cdr3_len
				else:
					#print "BADQRYMAP Failure to map to query for mode=",dm," V=",currentV," J=",currentJ," read=",currentQueryName," REFSTART=",ref_cdr3_start,"QRYSTART=",qry_cdr3_start,"REFEND=",ref_cdr3_end,"QRYEND=",qry_cdr3_end
					pass
			else:
				#print "BADREFMAP mode=",dm," refVCDR3=(-1) for ",currentV," = ",ref_cdr3_start," or refJCDR3 ",currentJ," = ",ref_cdr3_end
				pass
	else:
		#print "Ref names ",currentV," and ",currentJ," don't appear alleleic!"
		pass
	return cdr3_hist

