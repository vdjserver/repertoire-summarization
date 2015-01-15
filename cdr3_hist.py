#!/usr/bin/env python

import os
import re
import glob
from os.path import basename
from segment_utils import getFastaListOfDescs,getQueryIndexGivenSubjectIndexAndAlignment,getAdjustedCDR3StartFromRefDirSetAllele,getADJCDR3EndFromJAllele,getEmptyRegCharMap,alleleIsTR,getTheFrameForThisReferenceAtThisPosition,getTheFrameForThisJReferenceAtThisPosition
from imgt_utils import imgt_db
from utils import *
import vdjml
from vdjml_utils import getTopVDJItems,getHitInfo
from igblast_parse import rev_comp_dna
import argparse
from alignment import alignment,CodonAnalysis,codonAnalyzer








#utility to reconstruct alignments from a BTOP and data map and add it to the 
#data map and return the data map
def addAlignmentsPreCDR3(dataMap,alleleName,imgtdb_obj,organism,query_record):
	btop=dataMap['btop']
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
			self.modes=list()
			for line in reader:
				line=line.strip()
				line_pieces=line.split('\t')
				#print "got line=",line," and line_pieces=",line_pieces
				if(line_num==1 and len(line_pieces)>=2):
					for m in range(1,len(line_pieces)):
						self.modes.append(line_pieces[m].strip())
				elif(line_num==1 and len(line_pieces)<2):
					raise Exception('Error, invalid header in histogram file '+str(infile))
				elif(len(line_pieces)==len(self.modes)+1):
					#got data!
					for lp in line_pieces:
						if(not(self.appearsAsInt(lp))):
							raise Exception("Error, data ",lp," on line # ",line_num," appears as non-integral!  Integer values are expected!")
					val=line_pieces[0]
					for m in range(0,len(self.modes)):
						#print "for val=",val," and mode=",self.modes[m]," need to inc count ",line_pieces[m+1]
						for i in range(int(line_pieces[m+1])):
							mode=self.modes[m]
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


	#write as STREAM JSON
	def JSONIFYStreamAsString(self):
		self.modes.sort()
		mode_stream_list=list()
		for mode in self.modes:
			#print "Now JSONIFying for MODE=",mode
			min_val=self.gminVal()
			max_val=self.gmaxVal()
			#print "min and max are ",min_val,max_val
			round_num=0
			if((type(min_val)!=int) or (type(max_val)!=int)):
				min_val=(-1)
				max_val=(-1)
				print "Warning, could not determine min/max values for CDR3 lengths!  Output may not be defined!"			
			x_y_vals_array=list()
			for v in range(min(min_val,max_val),max(min_val,max_val)+1):
				x_y_str=""
				x_y_str+="\"x\": "+str(v)
				x_y_str+=","
				actual_val=0
				if(v in self.count_map[mode]):
					actual_val=self.count_map[mode][v]
				x_y_str+="\"y\": "+str(actual_val)
				x_y_vals_array.append(x_y_str)
			obj_sep="}\n,\n{"
			x_y_vals_array_str=obj_sep.join(x_y_vals_array)
			stream_str="{\n\"key\":\""+str(mode)+"\",\n"
			stream_str+="\"values\": [ \n{"+x_y_vals_array_str+"}\n]\n}"
			mode_stream_list.append(stream_str)
		stream_sep=","
		JSON="["+stream_sep.join(mode_stream_list)+"]"
		return JSON
#The sample data is 
#[
# {
#  "key": "KABAT",
#  "values": [
#   {
#    "x": 0,
#    "y": 0.14202716750712038
#   },
#   {
#    "x": 1,
#    "y": 0.17686738381528422
#   },
#   {
#    "x": 2,
#    "y": 0.2520439938176909
#   },


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
		self.modes.sort()
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




def getEmptyCDR3Map():
	domain_modes=get_domain_modes()
	cdr3_hist=dict()
	for dm in domain_modes:
		cdr3_hist[dm+'_from']=(-1)
		cdr3_hist[dm+'_to']=(-1)
		cdr3_hist[dm]=(-1)
	return cdr3_hist
		
		







def VJRearrangementInFrameTest(vInfo,jInfo,imgtdb_obj,organism):
	if(vInfo==None or jInfo==None):
		#need valid data to test. return false in this case
		return False
	#use V and J frame
	if(jInfo==None):
		#no J means no prod. rearrangment???
		return False
	s_end=int(vInfo['s. end'])
	q_end=int(vInfo['q. end'])
	refName=vInfo['subject ids']
	s_end_frame=getTheFrameForThisReferenceAtThisPosition(refName,organism,imgtdb_obj,s_end)
	q_end_frame=s_end_frame
	print "========================================================="
	print "\n\n\nLooking at ",vInfo['query id']," for FRAME test...."
	#print "q_start_frame (spos=",s_end,") ",q_end_frame
	q_end_j=int(jInfo['q. end'])
	q_bgn_j=int(jInfo['q. start'])
	print "First test...."
	if(q_end_j<=q_end):
		print "WAY TOO SHORT or BAD ALIGNMENT"
		print "IF Q END IN J IS LESS THAN OR EQUAL TO Q END IN V\n\n\n"
		return False
	else:
		print "Second test"
		j_bgn_j=int(jInfo['s. start'])
		j_bgn_frame=getTheFrameForThisJReferenceAtThisPosition(jInfo['subject ids'],organism,imgtdb_obj,j_bgn_j)
		print "v_end_frame is ",s_end_frame
		print "j_bgn_frame is ",j_bgn_frame
		if(j_bgn_frame!=None):
			print "passed into final...."
			print "V sub and qry end are : ",s_end," and ",q_end
			print "J sub and query bgn are : ",j_bgn_j
			expected_frame_based_on_V=(q_end_frame+(q_bgn_j-q_end))%3
			print "expected_frame_based_on_V=",expected_frame_based_on_V
			expected_frame_based_on_J=j_bgn_frame
			print "expected_frame_based_on_J=",expected_frame_based_on_J,"\n\n\n\n\n"
			if(expected_frame_based_on_J==expected_frame_based_on_V):
				#if both V and J impose the same frame, then the read should NOT be filtered
				return True
			else:
				return False
		return None





#given info maps for V and J and the return a 
#dictionary with kabat and imgt lengths
#return other related information as well
def CDR3LengthAnalysis(vMap,jMap,organism,imgtdb_obj):
	#currentQueryName=str(currentQueryName.strip())
	currentV=vMap['subject ids']
	currentJ=jMap['subject ids']
	cdr3_hist=getEmptyCDR3Map()
	cdr3_hist['Missing CYS']=True
	cdr3_hist['Missing TRP/PHE']=True
	cdr3_hist['Out-of-frame junction']=None
	cdr3_hist['Out-of-frame CDR3']=None
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
			if(alleleIsTR(currentV) and dm=="kabat"):
				#no kabat analysis for CDR3 TR!
				return cdr3_hist
			#print "processing in ",dm
			if(not currentV in rsmap[dm]):
				#print currentV,"not in lookup for dm=",dm
				ref_cdr3_start=getAdjustedCDR3StartFromRefDirSetAllele(currentV,imgtdb_obj,organism,dm)
				rsmap[dm][currentV]=ref_cdr3_start
			else:
				ref_cdr3_start=rsmap[dm][currentV]
			#print "After initial retrieval, the ref CDR3 start for ",currentV," is ",ref_cdr3_start
			if(not currentJ in remap[dm]):
				#print currentJ,"not in lookup for dm=",dm
				ref_cdr3_end=getADJCDR3EndFromJAllele(currentJ,imgtdb_obj,organism,dm)
				remap[dm][currentJ]=ref_cdr3_end
			else:
				ref_cdr3_end=remap[dm][currentJ]
			if(ref_cdr3_start!=(-1) and ref_cdr3_end!=(-1)):
				if(dm=="imgt"):
					cdr3_hist['Out-of-frame junction']=not(VJRearrangementInFrameTest(vMap,jMap,imgtdb_obj,organism))
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
				ref_cdr3_start-=1
				qry_cdr3_start=getQueryIndexGivenSubjectIndexAndAlignment(vq_aln,vs_aln,vq_f,vq_t,vs_f,vs_t,ref_cdr3_start)
				if(qry_cdr3_start!=(-1)):
					qry_cdr3_start+=1
				if(qry_cdr3_start!=(-1) and dm=='imgt'):
					#code for test for CYS-104 found
					#print "For ",currentV," the ref CDR3 start is ",ref_cdr3_start
					ref_cys_start=ref_cdr3_start-2
					ref_cys_end=ref_cdr3_start
					v_aln_obj=alignment(vq_aln,vs_aln,vq_f,vq_t,vs_f,vs_t)
					cys_aln=v_aln_obj.getSubAlnInc(ref_cys_start,ref_cys_end,"subject")
					#print "CYS ALN="
					#print "\n",cys_aln.getNiceString()
					query_cys_na=cys_aln.q_aln
					ref_cys_na=cys_aln.s_aln
					if(len(query_cys_na)==3  and len(ref_cys_na)==3):
						if(codonAnalyzer.is_unambiguous_codon(query_cys_na)==True):
							query_cys_trx=codonAnalyzer.fastTrans(query_cys_na)
							#print "The fast trans (with query=",vMap['query id']," and ref=",currentV,")=",query_cys_trx
							#print "\n\n\n\n\n\n\n===========================================================\n\n\n\n\n\n\n"
							#cdr3_hist['CYS']=true
							if(query_cys_trx=='C'):
								cdr3_hist['Missing CYS']=False
							else:
								cdr3_hist['Missing CYS']=True
				ref_cdr3_end+=1
				qry_cdr3_end=getQueryIndexGivenSubjectIndexAndAlignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t,ref_cdr3_end,"left")
				if(qry_cdr3_end!=(-1)):
					if(dm=='imgt'):
						#code for test for J-TRP/J-PHE found
						ref_trp_start=ref_cdr3_end
						ref_trp_end=ref_trp_start+2
						j_aln=alignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t)
						trp_aln=j_aln.getSubAlnInc(ref_trp_start,ref_trp_end,"subject")
						trp_q_na=trp_aln.q_aln
						trp_s_na=trp_aln.s_aln
						if(len(trp_q_na)==3 and len(trp_s_na)==3):
							#print "For ",vMap['query id']," got TRP-ALN : "
							#print "\n",trp_aln.getNiceString()
							qry_trp_trx=codonAnalyzer.fastTrans(trp_q_na)
							#print "Fast trans=",qry_trp_trx
							#print 
							#print "\n"
							if(qry_trp_trx=='W' or qry_trp_trx=='F'):
								cdr3_hist['Missing TRP/PHE']=False
							else:
								cdr3_hist['Missing TRP/PHE']=True
					qry_cdr3_end-=1
				if(qry_cdr3_start!=(-1) and qry_cdr3_end!=(-1)):
					if(dm=="imgt"):
						if(((qry_cdr3_end-qry_cdr3_start+1)%3)==0):
							cdr3_hist['Out-of-frame CDR3']=False
						else:
							cdr3_hist['Out-of-frame CDR3']=True
					#query_coding_seq=query_seq_map[currentQueryName]
					#coding_seq=query_coding_seq[(qry_cdr3_start-1):(qry_cdr3_end-1)]
					if(qry_cdr3_start<=qry_cdr3_end):
						#good
						cdr3_len=qry_cdr3_end-qry_cdr3_start+1
						#print "mode=",dm,"read=",vMap['query id']
						#print "For ",currentV," the ref_cdr3_bgn is ",ref_cdr3_start,"!"
						#print "For ",currentJ," the ref_cdr3_end is ",ref_cdr3_end,"!"
						#print "read start and end are ",qry_cdr3_start," and ",qry_cdr3_end,"\n\n\n\n"
						cdr3_hist[dm]=cdr3_len
						cdr3_hist[dm+'_from']=qry_cdr3_start
						cdr3_hist[dm+'_to']=qry_cdr3_end
					else:
						#messed up alignment presumably due to overlap! or in the RARE case where J aligns before V in the read!
						#print "messed up alignment presumably due to overlap! or in the RARE case where J aligns before V in the read!"
						pass
				else:
					#print "BADQRYMAP Failure to map to query for mode=",dm," V=",currentV," J=",currentJ," read=",vMap['query id'],"  REFSTART=",ref_cdr3_start,"QRYSTART=",qry_cdr3_start,"REFEND=",ref_cdr3_end,"QRYEND=",qry_cdr3_end
					pass
			else:
				#print "BADREFMAP mode=",dm," VALLELE=",currentV," and JALLELE=",currentJ," refVCDR3=(-1) for ",currentV," = ",ref_cdr3_start," or refJCDR3 ",currentJ," = ",ref_cdr3_end
				pass
	else:
		print "Ref names ",currentV," and ",currentJ," don't appear alleleic!"
		pass
	#print "***************************\n"
	#print "RETURNING THIS CDR3_HIST for ",vMap['query id'],":"
	#printMap(cdr3_hist)
	#print "***************************\n"	
	return cdr3_hist


#program to merge CDR3 histograms for kabat and imgt modes
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Merge multiple CDR3 length histograms into a single histogram.  Write the merged result to stdout')
	parser.add_argument('cdr3_hist_in',type=str,nargs='+',help="path(s) to CDR3 histograms to merge.  At least one is requried!")
	parser.add_argument('-j', action='store_true',help="Write the output as a 'stream-layer' formatted JSON string")
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
		if(args.j):
			#print "GOT LAYER!"
			print main_hist.JSONIFYStreamAsString()
		else:
			main_hist.writeToFile("/dev/stdout")
	else:
		parser.print_help()



