#!/usr/bin/env python

import os
import re
import glob
from os.path import basename
from segment_utils import getFastaListOfDescs,getQueryIndexGivenSubjectIndexAndAlignment,getAdjustedCDR3StartFromRefDirSetAllele,getADJCDR3EndFromJAllele,looksLikeAlleleStr,getEmptyRegCharMap
from imgt_utils import imgt_db
from utils import read_fasta_file_into_map,biopythonTranslate,rev_comp_dna,printMap,get_domain_modes,repStr
from igblast_utils import printNiceAlignment,buildAlignmentWholeSeqs,buildAlignmentWholeSeqsDirect
import vdjml
from vdjml_utils import getTopVDJItems,getHitInfo
import argparse




#how many deletions???
def getNumDashInStr(s):
	orig_len=len(s)
	s_no_dash=re.sub(r'\-','',s)
	nodash_len=len(s_no_dash)
	return orig_len-nodash_len


#from an alignment mark where CDR3 is with X
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


#utility to reconstruct alignments from a BTOP and data map and add it to the 
#data map and return the data map
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



#wrapper for getting CDR3 length using pyVDJML objects as input
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
	modes=get_domain_modes()
	for mode in modes:
		empty_map[mode]=(-1)
		empty_map[mode+'_from']=(-1)
		empty_map[mode+'_to']=(-1)
		empty_map[mode+'_to']=(-1)
		empty_map['qry_rev']=False
	if(vAllele!=None and jAllele!=None):
		vData=getHitInfo(read_result_obj,meta,vAllele)
		vData=addAlignmentsPreCDR3(vData,vAllele,imgtdb_obj,organism,query_record)
		jData=getHitInfo(read_result_obj,meta,jAllele)
		jData=addAlignmentsPreCDR3(jData,jAllele,imgtdb_obj,organism,query_record)
		#printMap(vData)
		#printMap(jData)
		cdr3_analysis_map=CDR3LengthAnalysis(vData,jData,organism,imgtdb_obj)
		return cdr3_analysis_map
		#printMap(cdr3_analysis_map)
	else:
		#print "insufficient data!"
		return empty_map
	pass


#class for CDR3 histogram
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

	#verify integer string
	def appearsAsInt(self,s):
		ire=re.compile(r'^\-?\d+$')
		if(re.match(ire,s)):
			return True
		else:
			return False


	#read a file into the object!
	def read_from_file(self,infile):
		line_num=1
		try:
			reader=open(infile,'r')
			modes=list()
			for line in reader:
				line=line.strip()
				line_pieces=line.split('\t')
				#print "got line=",line," and line_pieces=",line_pieces
				if(line_num==1 and len(line_pieces)>=2):
					for m in range(1,len(line_pieces)):
						modes.append(line_pieces[m].strip())
				elif(line_num==1 and len(line_pieces)<2):
					raise Exception('Error, invalid header in histogram file '+str(infile))
				elif(len(line_pieces)==len(modes)+1):
					#got data!
					for lp in line_pieces:
						if(not(self.appearsAsInt(lp))):
							raise Exception("Error, data ",lp," on line # ",line_num," appears as non-integral!  Integer values are expected!")
					val=line_pieces[0]
					for m in range(0,len(modes)):
						#print "for val=",val," and mode=",modes[m]," need to inc count ",line_pieces[m+1]
						for i in range(int(line_pieces[m+1])):
							mode=modes[m]
							#print "calling inc with mode=",mode," val=",val
							self.inc(mode,val)
				else:
					print "warning, bad data in file ",infile," line ",line_num
				line_num+=1				
			reader.close()
		except Exception as my_exception:
			print my_exception," line #",line_num


	#increment
	def inc(self,mode,val):
		val=int(val)
		#print "INC CALLED"
		if(not val in self.count_map[mode]):
			self.count_map[mode][val]=1
		else:
			self.count_map[mode][val]+=1

	#get min val over all vals
	def gminVal(self):
		min_val=1000
		for mode in self.modes:
			for val in self.count_map[mode]:
				if(val<min_val):
					min_val=val
		return min_val

	#get max val over all vals
	def gmaxVal(self):
		max_val=(-1)
		for mode in self.modes:
			for val in self.count_map[mode]:
				if(val>max_val):
					max_val=val
		return max_val
							


	#merge from other hist
	def merge(self,other):
		#print "MERGE IS CALLED"
		other_map=other.count_map
		#print "other_map is ",other_map
		for mode in other_map:
			if(not mode in self.count_map):
				self.count_map[mode]=dict()
			for val in other_map[mode]:
				#print "for mode=",mode," at val=",val
				other_count=other_map[mode][val]
				#print "the count is ",other_count
				for i in range(other_count):
					self.inc(mode,val)

	#print maps basic
	def printMaps(self):
		for mode in self.modes:
			#print "SHOWING MAP FOR ",mode
			for val in self.count_map[mode]:
				print val,"\t",self.count_map[mode][val]

	#write a basic histogram to a file!
	def writeToFile(self,ofile):
		writer=open(ofile,'w')
		min_val=self.gminVal()
		max_val=self.gmaxVal()
		#print "min and max are ",min_val,max_val
		round_num=0
		if((type(min_val)!=int) or (type(max_val)!=int)):
			min_val=(-1)
			max_val=(-1)
			print "Warning, could not determine min/max values for CDR3 lengths!  Output may not be defined!"
		for v in range(min(min_val,max_val),max(min_val,max_val)+1):
			if(round_num==0):
				writer.write("CDR3_LENGTH")
				for mode in self.modes:
					writer.write("\t"+mode)
				writer.write("\n")
			if(v!=(0)):
				writer.write(str(v))
				for mode in self.modes:
					actual_val=0
					if(v in self.count_map[mode]):
						actual_val=self.count_map[mode][v]
					writer.write("\t"+str(actual_val))
				writer.write("\n")
			round_num+=1
		writer.close()



#cache for CDR3 start/end
rsmap=dict()
remap=dict()
domain_list=get_domain_modes()
for d in domain_list:
       rsmap[d]=dict()
       remap[d]=dict()



#given info maps for V and J and the returnd a 
#dictionary with kabat and imgt lengths
def CDR3LengthAnalysis(vMap,jMap,organism,imgtdb_obj):
	#currentQueryName=str(currentQueryName.strip())
	currentV=vMap['subject ids']
	currentJ=jMap['subject ids']
	cdr3_hist=dict()
	if(vMap['query id'].find("reversed|")==0):
		cdr3_hist['qry_rev']=True
	else:
		cdr3_hist['qry_rev']=False
	if(looksLikeAlleleStr(currentV) and looksLikeAlleleStr(currentJ)):
		#print "WE'RE IN BUSINESS!"
		domain_modes=get_domain_modes()
		for dm in domain_modes:
			cdr3_hist[dm]=(-1)
			cdr3_hist[dm+'_from']=(-1)
			cdr3_hist[dm+'_to']=(-1)
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
				qry_cdr3_end=getQueryIndexGivenSubjectIndexAndAlignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t,ref_cdr3_end,"left")
				if(qry_cdr3_start!=(-1) and qry_cdr3_end!=(-1)):
					#query_coding_seq=query_seq_map[currentQueryName]
					#coding_seq=query_coding_seq[(qry_cdr3_start-1):(qry_cdr3_end-1)]
					if(qry_cdr3_start<qry_cdr3_end):
						#good
						cdr3_len=qry_cdr3_end-qry_cdr3_start+1
						#translation=biopythonTranslate(coding_seq)
						#print "the coding seq ("+dm+") is : ",coding_seq
						#print "The translation ("+dm+") is : ",translation
						#print "CDR3_LEN ("+dm+") ="+str(cdr3_len)
						cdr3_hist[dm]=cdr3_len
						cdr3_hist[dm+'_from']=qry_cdr3_start
						cdr3_hist[dm+'_to']=qry_cdr3_end
					else:
						#messed up alignment presumably due to overlap!
						cdr3_hist[dm]=(-1)
						#cdr3_hist[dm+'_from']=(-1)
						#cdr3_hist[dm+'_to']=(-1)
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


#program to merge CDR3 histograms for kabat and imgt modes
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Merge multiple CDR3 length histograms into a single histogram.  Write the merged result to stdout')
	parser.add_argument('cdr3_hist_in',type=str,nargs='+',help="path(s) to CDR3 histograms to merge.  At least one is requried!")
	args=parser.parse_args()
	if(args):
		#print "in args"
		cdr3_hist_in_list=args.cdr3_hist_in
		main_hist=histoMapClass(get_domain_modes())
		#print "got ",cdr3_hist_in_list
		for infile in cdr3_hist_in_list:
			#print "now to analyze for ",infile
			temp_hist=histoMapClass(get_domain_modes())
			temp_hist.read_from_file(infile)
			main_hist.merge(temp_hist)
		main_hist.writeToFile("/dev/stdout")
	else:
		parser.print_help()



