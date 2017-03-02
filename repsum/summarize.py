"""
Repertoire summarization
"""

import vdjml
#from igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs,rev_comp_dna
from vdjml.igblast_parse import rev_comp_dna
#from utils import printMap,get_domain_modes,biopythonTranslate,makeAllMapValuesVal,getNumberBpInAlnStr,repStr
from utils import *
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
#from vdjml_utils import getTopVDJItems,getHitInfo,getVDJServerRegionAlignmentFromLargerVAlignmentPyObj,getAlignmentString,getTypeNameScoreMap,getVDJTieMap
from vdjml_utils import *
from Bio import SeqIO
from segment_utils import getVRegionsList,getEmptyRegCharMap,getAdjustedCDR3StartFromRefDirSetAllele,getTheFrameForThisReferenceAtThisPosition,getVRegionStartAndStopGivenRefData,getADJCDR3EndFromJAllele,alleleIsTR
from segment_utils import recombFreqManager
from char_utils import getNumberBaseSubsFromBTOP,getNumberIndelsFromBTOP,getIndelMapFromBTOP
from alignment import alignment,CodonAnalysis,codonAnalyzer
from CharacterizationThread import CharacterizationThread
from codon_analysis import *
import re
import Queue
import threading
import time
import multiprocessing
import math

global_key_base="vdj_server_ann_"
myCodonCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")
release_info_tag=None
release_info_hash=None




#initialize the sums to 0 for accumulation during execution
#and output at program termination

#whole read
num_seqs_pf=0
whole_read_b2b_aln=0
whole_read_b2b_sbst=0

#b2b region values
framework_b2b_aln_imgt=0
framework_b2b_sbst_imgt=0
cdr_b2b_aln_imgt=0
cdr_b2b_subst_imgt=0
framework_b2b_aln_kabat=0
framework_b2b_sbst_kabat=0
cdr_b2b_aln_kabat=0
cdr_b2b_subst_kabat=0

#r:s for region variables
framework_r_sum_imgt=0
framework_s_sum_imgt=0
cdr_r_sum_imgt=0
cdr_s_sum_imgt=0
framework_r_sum_kabat=0
framework_s_sum_kabat=0
cdr_r_sum_kabat=0
cdr_s_sum_kabat=0	



def loadReleaseInfoTagAndHash():
	global release_info_tag
	global release_info_hash
	release_info_tag=getVersionInfo("desc")
	release_info_hash=getVersionInfo("hash")
	


# use a read_result_ojbect return several things:
# 1) the segment counts (a list of 1, 2, or 3 items) , 
# 2) the cdr3_lengths (kabat and imgt) ( a dict() with two keys)
# 3) a clone of the object, but with additional (with additional read information)
def rep_char_read(read_result_obj,meta,organism,imgtdb_obj,read_rec,skip_char=False):
	#print "In loop analyzing ",read_result_obj.id()


	#prepare an object to return
	return_obj=dict()

	#retrieve the top VDJ
	topVDJ=getTopVDJItems(read_result_obj,meta)
	return_obj['VDJ']=topVDJ

	#retrieve the CDR3 lengths
	cdr3_length_results=CDR3LengthAnalysisVDMLOBJ(read_result_obj,meta,organism,imgtdb_obj,read_rec)
	return_obj['cdr3_length_results']=cdr3_length_results

	#perform own annotation
	read_ann_map=readAnnotate(read_result_obj,meta,organism,imgtdb_obj,read_rec,cdr3_length_results,skip_char)
	return_obj['ann_map']=read_ann_map

	return return_obj
	





#specialized annotation for CDR3 region!
def readAnnotate_cdr3(read_result_obj,meta,organism,imgtdb_obj,read_rec,read_ann_map,cdr3_length_results):
	global global_key_base
	#print "SOME CDR3LEN RES"
	#printMap(cdr3_length_results)
	#print "\n\n\n\n\n"
	if('Out-of-frame junction' in cdr3_length_results):
		read_ann_map['Out-of-frame junction']=cdr3_length_results['Out-of-frame junction']
	else:
		read_ann_map['Out-of-frame junction']=None
	if('Out-of-frame CDR3' in cdr3_length_results):
		read_ann_map['Out-of-frame CDR3']=cdr3_length_results['Out-of-frame CDR3']
	else:
		read_ann_map['Out-of-frame CDR3']=None
	if('Missing CYS' in cdr3_length_results):
		read_ann_map["Missing CYS"]=cdr3_length_results['Missing CYS']
	else:
		read_ann_map["Missing CYS"]=None
	if('Missing TRP/PHE' in cdr3_length_results):
		read_ann_map["Missing TRP/PHE"]=cdr3_length_results['Missing TRP/PHE']
	else:
		read_ann_map["Missing TRP/PHE"]=None
	for mode in get_domain_modes():
		f=cdr3_length_results[mode+'_from']
		t=cdr3_length_results[mode+'_to']
		len_key_aa="CDR3 AA length ("+mode+")"
		aa_key="CDR3 AA ("+mode+")"
		na_key="CDR3 NA ("+mode+")"
		len_key_na="CDR3 NA length ("+mode+")"
		read_ann_map[mode+'_from']=f
		read_ann_map[mode+'_to']=t
		if(f!=(-1) and t!=(-1) and cdr3_length_results[mode]!=(-1) ):
			qw=str(read_rec.seq)
			if(cdr3_length_results['qry_rev']):
				qw=rev_comp_dna(qw)
			query_cdr3=qw[f-1:t]
			read_ann_map[global_key_base+mode+'_cdr3_na']=query_cdr3
			read_ann_map[na_key]=query_cdr3.upper()
			read_ann_map[aa_key]=biopythonTranslate(read_ann_map[global_key_base+mode+'_cdr3_na'])
			#add in lengths!
			read_ann_map[global_key_base+mode+'_cdr3_na_len']=len(query_cdr3)
			read_ann_map[len_key_aa]=len(read_ann_map[aa_key])
			read_ann_map[len_key_na]=len(read_ann_map[na_key])
			#this stuff below was once in here becuase we might've done some CDR3 specific annotation...but
			#it is commented out for now because such goals are either postponed or canceled
			#vMap=getHitInfo(read_result_obj,meta,read_ann_map['top_V'],read_rec,imgtdb_obj,organism)
			#dMap=None
			#if(read_ann_map['top_D'] is not None):
			#	dMap=getHitInfo(read_result_obj,meta,read_ann_map['top_D'],read_rec,imgtdb_obj,organism)
			#jMap=getHitInfo(read_result_obj,meta,read_ann_map['top_J'],read_rec,imgtdb_obj,organism)
			#cdr3RegCharAnalysis(vMap,dMap,jMap,mode,cdr3_length_results,organism,imgtdb_obj)
			#getCDR3RegionSpecificCharacterizationSubAln(vMap,dMap,jMap,organism,imgtdb_obj,mode,read_rec)
			#sys.exit(0)
		else:
			read_ann_map[global_key_base+mode+'_cdr3_na']=""
			read_ann_map[aa_key]=""
			#add in lengths!
			read_ann_map[global_key_base+mode+'_cdr3_na_len']=""
			read_ann_map[len_key_aa]=""
			pass
	return read_ann_map





def getIndelsFlags(vInfo,dInfo,jInfo):
	strs_to_analyze=list()
	basicIndelsFlags=False
	framePreservingIndelsFlag=False
	#printMap(vInfo)
	if(vInfo!=None):
		v_q_seq=vInfo['query seq']
		v_s_seq=vInfo['subject seq']
		strs_to_analyze.append(v_q_seq)
		strs_to_analyze.append(v_s_seq)
	if(dInfo!=None):
		d_q_seq=dInfo['query seq']
		d_s_seq=dInfo['subject seq']
		strs_to_analyze.append(d_q_seq)
		strs_to_analyze.append(d_s_seq)
	if(jInfo!=None):
		j_q_seq=jInfo['query seq']
		j_s_seq=jInfo['subject seq']
		strs_to_analyze.append(j_q_seq)
		strs_to_analyze.append(j_s_seq)
	gap_re=re.compile(r'([\-]{1,})')	
	#first look for >=1 indels
	for s in strs_to_analyze:
		#print "Looking for gap in ",s
		if(s.find("-")!=(-1)):
			#found a gap!
			#print "FOUND A GAP"
			basicIndelsFlags=True
		else:
			#print "didn't find a gap!"
			pass
	#if at least one was found then look at the frame-preserving status
	if(basicIndelsFlags):
		#set the p flag to TRUE unless at least on frame-destroying indel is found
		framePreservingIndelsFlag=True
		for s in strs_to_analyze:
			gap_groups=re.findall(gap_re,s)
			for gap_group in gap_groups:
				gap_len=len(gap_group)
				if(gap_len%3!=0):
					#it's not frame preserving!
					framePreservingIndelsFlag=False
				else:
					#found a frame preserving gap!
					pass
	gap_ret_package=[basicIndelsFlags,framePreservingIndelsFlag]
	return gap_ret_package



def getValnWholeSeqStopFlag(vInfo,dInfo,jInfo,imgtdb_obj,organism,annMap,seq_rec):
	if(vInfo==None):
		#if no V alignment how can there be a frame, a J alignment, or a stop codon?
		return None
	else:
		vAlnObj=alignment(vInfo['query seq'],vInfo['subject seq'],vInfo['q. start'],vInfo['q. end'],vInfo['s. start'],vInfo['s. end'])
		#print "Looking at STOP log for ",str(seq_rec.id)
		#print "V allele is ",vInfo['subject ids']
		if(jInfo!=None and 'subject ids' in jInfo):
			#print "J allele is ",jInfo['subject ids']
			pass
		if(dInfo!=None and 'subject ids' in dInfo):
			#print "D allele is ",dInfo['subject ids']
			pass
			
		#print "The V alignment plain is \n",vAlnObj.getNiceString()
		vAlnNoGapStart=vAlnObj.getGEQAlignmentFirstNonGap()
		#print "The V alignment non-gap ('de-gapped' begin') start is ",vAlnNoGapStart.getNiceString()
		s_start=int(vInfo['s. start'])
		q_start=int(vInfo['q. start'])
		s_start_frame=getTheFrameForThisReferenceAtThisPosition(vInfo['subject ids'],organism,imgtdb_obj,s_start)
		info_to_use=None
		if(jInfo==None):
			#see if D is available
			if(dInfo==None):
				info_to_use=vInfo
			else:
				info_to_use=dInfo
		else:
			info_to_use=jInfo
		#use max in case there is germline alignment overlap
		#print "S start frame is ",s_start_frame
		q_end=max(int(info_to_use['q. end']),int(vInfo['q. end']))
		if(not(vInfo['is_inverted'])):
			#print "Computing to_trans regular block...."
			to_trans=str(seq_rec.seq)[(q_start-1):q_end].upper()
		else:
			proper=rev_comp_dna(str(seq_rec.seq))
			to_trans=str(proper[q_start-1:q_end]).upper()
			#print "computing to_trans inverted block"
		#print "s_start_frame is ",s_start_frame
		if(s_start_frame!=0):
			to_trans=to_trans[(3-s_start_frame):]
		#print "To be translated for (adjusted for start frame) (pre-trim by alignment) read=",str(seq_rec.id)," is\n",to_trans
		to_trans_len=len(to_trans)
		if((to_trans_len%3)!=0):
			necessary_trim=(to_trans_len%3)
			trim_to_trans=to_trans[:to_trans_len-necessary_trim]
		else:
			trim_to_trans=to_trans
		#print "The trim_to_trans is (trimmed by V alignment and the more extensive of D or J) is\n",trim_to_trans
		translation=codonAnalyzer.fastTransStr(trim_to_trans)
		stopFlag=stringContainsAStar(translation)
		#print "The translation : ",translation
		#print "The stop flag : ",stopFlag
		#print "\n\n\n\n"
		return stopFlag

		
		




#get pre and post CDR3 alignment objects.
#this can be use in some whole_seq stuff
#or in just pre or post CDR3 stuff
def getPrePostCDR3AlnObjs(vInfo,jInfo,imgtdb_obj,organism,annMap):
	if(not(vInfo==None)):
		#at least work with V
		v_cdr3_start=getAdjustedCDR3StartFromRefDirSetAllele(vInfo['subject ids'],imgtdb_obj,organism,"imgt")
		if(v_cdr3_start==(-1)):
			#if CDR3 start is not found, just use the entire thing!
			v_cdr3_start=vInfo['s. end']+1
		vAlnObj=alignment(vInfo['query seq'],vInfo['subject seq'],vInfo['q. start'],vInfo['q. end'],vInfo['s. start'],vInfo['s. end'])
		preCDR3Aln=vAlnObj.getAlnAtAndCond(v_cdr3_start-1,"subject","leq")
		preCDR3Aln.setSFM(getTheFrameForThisReferenceAtThisPosition(vInfo['subject ids'],organism,imgtdb_obj,preCDR3Aln.s_start))
		if(jInfo==None):
			#if no info is found, return the value from just V
			return 	[preCDR3Aln,None]
		else:
			#work on V and J
			jAlnObj=alignment(jInfo['query seq'],jInfo['subject seq'],jInfo['q. start'],jInfo['q. end'],jInfo['s. start'],jInfo['s. end'])
			J_cdr3_end=getADJCDR3EndFromJAllele(jInfo['subject ids'],imgtdb_obj,organism,"imgt")
			#print "jAlnObj"
			#print jAlnObj.getNiceString()
			if(J_cdr3_end==(-1)):
				#could not detect CDR3 end! so just return V map data
				#return preCDR3Aln.characterize()
				return [preCDR3Aln,None]
			elif(jAlnObj.q_end<=preCDR3Aln.q_start or jAlnObj.q_end<=preCDR3Aln.q_end):
				#j alignment TOTALLY within V so just return V ; avoid double countings
				return [preCDR3Aln,None]
			else:
				postCDR3AlnStart=J_cdr3_end+1
				#okay, got J alignment post-CDR3 position
			if(postCDR3AlnStart>jAlnObj.s_end):
				#CDR3 start is after the s_aln end which means J aligns only in CDR3, ....so just return V map
				return [preCDR3Aln,None]
			#extract the subalignment at CDR3 end
			postCDR3AlnObj=jAlnObj.getAlnAtAndCond(postCDR3AlnStart,"subject","geq")
			if(postCDR3AlnObj.q_start<=preCDR3Aln.q_end):
				#this is in case V and J overlap in the CDR3 area alignment
				new_q_start=preCDR3Aln.q_end+1
				#print "It is deteced that this is in case V and J overlap in the CDR3 area alignment"
				#print "new q_start is ",new_q_start
				#adjust the J sub-alignment to start after the V ends
				postCDR3AlnObj=jAlnObj.getAlnAtAndCond(new_q_start,"query","geq")
			if(postCDR3AlnObj.s_start==postCDR3AlnStart):
				#great, set the FRAME mask to start with zero cause the alignment starts where it ideally should start
				#print "postCDR3AlnObj.s_start=",postCDR3AlnObj.s_start,", and postCDR3AlnStart=", postCDR3AlnStart," so normal case for SFM"
				postCDR3AlnObj.setSFM(0)
			elif(postCDR3AlnStart<postCDR3AlnObj.s_start):
				#if the alignment starts later, set the proper frame mask start
				#print "postCDR3AlnStart=",postCDR3AlnStart," and postCDR3AlnObj.s_start=",postCDR3AlnObj.s_start," so setting SFM as ",str(int(int(postCDR3AlnObj.s_start-postCDR3AlnStart)%3))
				postCDR3AlnObj.setSFM((postCDR3AlnObj.s_start-postCDR3AlnStart)%3)
			else:
				raise Exception('postCDR3AlnStart>postCDR3AlnObj.s_start ?????  Error in alignment class???')
			return [preCDR3Aln,postCDR3AlnObj]
	else:
		#nothing to work with!? :(
		return [None,None]

	
#return the char map for the V alignment
def returnVAlnAllCharMap(vInfo,imgtdb_obj,organism,read_result_obj=None,meta=None):
	global global_key_base
	emptyMap=getEmptyRegCharMap()
	if(vInfo==None):
		return emptyMap
	vAlnObj=alignment(vInfo['query seq'],vInfo['subject seq'],vInfo['q. start'],vInfo['q. end'],vInfo['s. start'],vInfo['s. end'])
	vAlnObj.setSFM(getTheFrameForThisReferenceAtThisPosition(vInfo['subject ids'],organism,imgtdb_obj,int(vInfo['s. start'])))
	if(vAlnObj==None):
		return emptyMap
	if(vInfo['q. start']>vInfo['q. end'] or vInfo['q. start']<=0 or vInfo['q. end']<=0):
		return emptyMap
	else:
		vCharMap=vAlnObj.characterize()
		if((read_result_obj is not None) and (meta is not None)):
			aa_mut_list=vAlnObj.AA_mut_list
			target_seg_match_name=vInfo['subject ids']
			print "Adding muts at ",read_result_obj.id()
			addAAMuts(read_result_obj,aa_mut_list,target_seg_match_name,meta)
			showAAMuts(read_result_obj)
		new_map=dict()
		for k in vCharMap:
			newKey=k
			newKey=newKey[0].upper()+newKey[1:]
			new_map[newKey+ " (over V)"]=vCharMap[k]
			if(k=="homology%"):
				new_map['V Sequence similarity']=vCharMap[k]
		new_map['V length (nucleotides)']=abs(vInfo['q. start']-vInfo['q. end'])
		return new_map




#return characterization for V and J (EXCLUDING CDR3 REGION)
#if CDR3 is unknown, just return an empty map! :(
def returnWholeSeqCharMap(vInfo,jInfo,imgtdb_obj,organism,annMap):
	global global_key_base
	cdr3_surround_aln=getPrePostCDR3AlnObjs(vInfo,jInfo,imgtdb_obj,organism,annMap)
	preCDR3Aln=cdr3_surround_aln[0]
	postCDR3AlnObj=cdr3_surround_aln[1]
	emptyMap=getEmptyRegCharMap()
	if(preCDR3Aln==None and postCDR3AlnObj==None):
		return emptyMap
	elif(preCDR3Aln!=None and postCDR3AlnObj==None):
		return preCDR3Aln.characterize()
	elif(preCDR3Aln==None and postCDR3AlnObj!=None):
		return emptyMap


	vCharMap=preCDR3Aln.characterize()
	#print "the vCharMap annotation map and alignment (for whole seq named ",vInfo['query id'],") is "
	#printMap(vCharMap)
	#print preCDR3Aln.getNiceString()
	#print "the jCharMap annotation map and alignment (for whole seq name ",jInfo['query id'],") is "
	##print postCDR3AlnObj.getNiceString()
	jCharMap=postCDR3AlnObj.characterize()
	#printMap(jCharMap)
	#to get/compute the 'whole' map, add up the "non-special" values
	#otherwise, compute them specially
	#toSkip=['bsb_freq','pct_id','ns_rto','indel_freq']
	toSkip=[]
	totCharMap=dict()
	for k in vCharMap:
		if(not(k in toSkip)):
			if(type(vCharMap[k])==type(jCharMap[k]) and type(vCharMap[k])!=str and vCharMap[k]!=None and jCharMap[k]!=None):
				totCharMap[k]=vCharMap[k]+jCharMap[k]
			else:
				totCharMap[k]=str(vCharMap[k])+" & "+str(jCharMap[k])
	if(vCharMap['Stop codons?']==None and jCharMap['Stop codons?']==None):
		totCharMap['Stop codons?']=None
	elif(vCharMap['Stop codons?']==True or jCharMap['Stop codons?']==True):
		totCharMap['Stop codons?']=True
	else:
		totCharMap['Stop codons?']=False
	
	if("CDR3 AA (imgt)" in annMap):
		starPos=annMap["CDR3 AA (imgt)"].find("*")
		if(starPos!=(-1)):
			totCharMap['CDR3 Stop codon?']=True
		else:
			totCharMap['CDR3 Stop codon?']=False
	else:
		totCharMap['CDR3 Stop codon?']=False
		#handle pct id , bsb_freq , indel_freq specially
	totQ=getNumberBpInAlnStr(postCDR3AlnObj.q_aln)+getNumberBpInAlnStr(preCDR3Aln.q_aln)
	totBS=totCharMap['base substitutions']
	if(totQ!=0):
		totCharMap['homology%']=float((float(totQ)-float(totBS))/float(totQ))
		totCharMap['base substitution freq%']=float(totCharMap['base substitutions'])/float(totQ)
		totCharMap['indel frequency']=float(totCharMap['insertion count']+totCharMap['deletion count'])/float(totQ+totCharMap['insertion count']+totCharMap['deletion count'])
	else:
		totCharMap['homology%']=0.0
		totCharMap['base substitution freq%']=0
		totCharMap['indel frequency']=0
	#handle ns_ratio
	if(totCharMap['synonymous base substitutions']!=0):
		totCharMap['ns_rto']=float(totCharMap['nonsynonymous base substitutions'])/float(totCharMap['synonymous base substitutions'])
	else:
		totCharMap['ns_rto']=0
	return totCharMap


	


#characterization_queue=multiprocessing.JoinableQueue()
characterization_queue=Queue.Queue()
characterization_queue_results=Queue.Queue()
alignment_output_queue=Queue.Queue()
characterization_thread_set=None
			
#given a read object, meta data, organism and database data and the read record and CDR3 data
#make a characterization map using the database/lookups
#of regions FR1, CDR1, FR2, CDR2, FR3
def readAnnotate(read_result_obj,meta,organism,imgtdb_obj,read_rec,cdr3_map,skip_char=False):
	#print "To annotate a read...."
	global global_key_base
	topVDJ=getTopVDJItems(read_result_obj,meta)
	annMap=dict()
	for seg in topVDJ:
		topkey=seg+" gene"
		if(topVDJ[seg] is not None):
			annMap[topkey]=topVDJ[seg]
		else:
			annMap[topkey]="None"

	vInfo=None
	dInfo=None
	jInfo=None
	noneSeg_flag=False
	for seg in topVDJ:
		#print "in second loop...."
		#print "seg=",seg
		if(topVDJ[seg] is not None):
			#print "to get info...."
			hit_info=getHitInfo(read_result_obj,meta,topVDJ[seg],read_rec,imgtdb_obj,organism)
			if(seg=="V" and not(hit_info==None)):
				vInfo=hit_info
			elif(seg=="J" and not(hit_info==None)):
				jInfo=hit_info
			elif(seg=="D" and not(hit_info==None)):
				dInfo=hit_info
		else:
			if(seg=='V'):
				noneSeg_flag=True


	#annotations such as NA and AA from CDR3 from cdr3_map
	annMap=readAnnotate_cdr3(read_result_obj,meta,organism,imgtdb_obj,read_rec,annMap,cdr3_map)

	gap_ret_package=getIndelsFlags(vInfo,dInfo,jInfo)
	annMap['Indels Found']=gap_ret_package[0]
	annMap['Only Frame-Preserving Indels Found']=gap_ret_package[1]


	#characterization beside segments/CDR3
	global characterization_thread_set
	global characterization_queue
	global characterization_queue_results
	global alignment_output_queue
	num_submitted_jobs=0
	if(not(skip_char)):
		#V seq characterization (OVER v) only
		#V_map=returnVAlnAllCharMap(vInfo,imgtdb_obj,read_result_obj,meta)
		V_map=returnVAlnAllCharMap(vInfo,imgtdb_obj,organism)
		#print "A V_MAP\n"
		#printMap(V_map)
		if(V_map is not None):
			for k in V_map:
				annMap[k]=V_map[k]

		#whole seq characterization (over V and J)
		#whole_char_map=returnWholeSeqCharMap(vInfo,jInfo,imgtdb_obj,organism,annMap)
		#for w in whole_char_map:
		#	#new_key=global_key_base+'whole_seq_'+w
		#	new_key=w+" (over V and J)"
		#	new_key=new_key[0].upper()+new_key[1:]
		#	annMap[new_key]=whole_char_map[w]
		whole_seq_stp_cdn_Tot_flag=getValnWholeSeqStopFlag(vInfo,dInfo,jInfo,imgtdb_obj,organism,annMap,read_rec)
		annMap['Stop Codon?']=whole_seq_stp_cdn_Tot_flag


		#productive rearrangement 
		#if IMGT CDR3 length is a multiple of 3 and it's not (-1), then consider it productive
		prodRearrangeKey="Productive CDR3 rearrangement (T/F)"
		if(not(noneSeg_flag)):
			if(cdr3_map['imgt']!=(-1) and cdr3_map['imgt']>0):
				if((cdr3_map['imgt'])%3==0):
					annMap[prodRearrangeKey]=True
				else:
					annMap[prodRearrangeKey]=False
			else:
				#annMap[prodRearrangeKey]="N/A"
				
				pass

		#regions characterization FR1, FR2, FR3, CDR1, CDR2 (imgt and kabat)
		mode_list=get_domain_modes()

		if(characterization_thread_set==None):
			num_cpus=multiprocessing.cpu_count()
			num_threads=max(1,num_cpus)
			#num_threads=1
			characterization_thread_set=set()
			for c in range(num_threads):
				#print "ABOUT TO CALL CONST...."
				temp_thread=CharacterizationThread(characterization_queue,characterization_queue_results,alignment_output_queue)
				#print "RETURNED FROM CONSTRUCTOR"
				temp_thread.setDaemon(True)
				#sys.exit(0)
				temp_thread.start()
				characterization_thread_set.add(temp_thread)
		if(not(noneSeg_flag)):
			for temp_thread in characterization_thread_set:
				temp_thread.clear_result()
			for mode in mode_list:
				if(not(mode=="kabat" and alleleIsTR(topVDJ[seg]))):
					for region in getVRegionsList():
						key_base=global_key_base+mode+"_"
						char_job_dict=dict()
						char_job_dict['refName']=topVDJ[seg]
						char_job_dict['read_result_obj']=read_result_obj
						char_job_dict['meta']=meta
						char_job_dict['organism']=organism
						char_job_dict['mode']=mode
						char_job_dict['region']=region
						char_job_dict['imgtdb_obj']=imgtdb_obj
						char_job_dict['wholeOnly']=False
						char_job_dict['read_rec']=read_rec
						char_job_dict['key_base']=key_base
						char_job_dict['region']=region
						char_job_dict['noneSeg_flag']=noneSeg_flag
						char_job_dict['read_name']=read_rec.id
						characterization_queue.put(char_job_dict)
						num_submitted_jobs+=1
				else:
					#if mode=kabat and type is TR, then skip the regions!
					pass
		else:
			#no annotation possible if V is empty!
			pass

	#here, get results from the job, putting them in the annMap
	if(num_submitted_jobs>0):
		#print "CHAR_QUEUE pre-join...."
		characterization_queue.join()
		get_res=0
		while(get_res<num_submitted_jobs):
			temp_results=characterization_queue_results.get()
			if(temp_results is not None):
				#print "NOW LOOKING AT TEMP_RESULTS:"
				#printMap(temp_results)
				#print "\n\n\n\n\n"
				#sys.exit(0)
				for temp_key in temp_results:
					annMap[temp_key]=temp_results[temp_key]			
			get_res+=1

	#get mutation map
	#mut_map_info=annotationMutationMap(vInfo,dInfo,jInfo,alignment_output_queue,num_submitted_jobs,imgtdb_obj,myCodonCounter,organism,read_rec,cdr3_map,0.85)
	mut_map_info=annotationMutationMap(vInfo,dInfo,jInfo,alignment_output_queue,num_submitted_jobs,imgtdb_obj,myCodonCounter,organism,read_rec,cdr3_map,0)
	qualifNote=mut_map_info[0]
	maps=mut_map_info[1]
	amino_map=maps['aminos']
	codon_map=maps['codons']
	amino_silent_map=maps['aminos_silent']
	codon_silent_map=maps['codons_silent']
	annMap["AGS Q Note"]=qualifNote
	annMap["Codon Muts (kabat)"]=codon_map
	annMap["Amino Muts (kabat)"]=amino_map
	annMap["Silent Codon Muts (kabat)"]=codon_silent_map
	annMap["Silent Amino Muts (kabat)"]=amino_silent_map



	#print "***************************************\n\n  A READ\n"
	#printMap(annMap)
	#print "***************************************\n\n\n"



	
	#here, compute the region total REPLACEMENTS, SILENTS
	#for computation of CDR R:S ratio and FR R:S ratio
	#Out-of-frame junction	Missing CYS	Missing TRP/PHE	Stop Codon?	Indels Found	Only Frame-Preserving Indels Found
	if(
		annMap['Out-of-frame junction']==False and 
		annMap['Missing CYS']==False and 
		annMap['Missing TRP/PHE']==False 
		and annMap['Stop Codon?']==False
		):
		indelComposite=True
		if(annMap['Indels Found']==False):
			#no indels so pass filter
			indelComposite=True
		else:
			#indels found
			if(annMap['Only Frame-Preserving Indels Found']==True):
				#they're frame preserving!
				indelComposite=True
			else:
				indelComposite=False
		if(indelComposite):

			global num_seqs_pf
			global whole_read_b2b_aln
			global whole_read_b2b_sbst

			#whole read
			num_seqs_pf+=1
			whole_read_num_key="Base substitutions (over V)"
			whole_read_dnm_key="Number base-to-base aligned (over V)"
			if(whole_read_num_key in annMap and whole_read_dnm_key in annMap):
				whole_num=annMap[whole_read_num_key]
				whole_dnm=annMap[whole_read_dnm_key]
				if(whole_num!=None and whole_dnm!=None):
					if(whole_dnm>=0):
						#whole read
						whole_read_b2b_aln+=whole_dnm
						whole_read_b2b_sbst+=whole_num



			#b2b region values
			global framework_b2b_aln_imgt
			global framework_b2b_sbst_imgt
			global cdr_b2b_aln_imgt
			global cdr_b2b_subst_imgt
			global framework_b2b_aln_kabat
			global framework_b2b_sbst_kabat
			global cdr_b2b_aln_kabat
			global cdr_b2b_subst_kabat

			#r:s for region variables
			global framework_r_sum_imgt
			global framework_s_sum_imgt
			global cdr_r_sum_imgt
			global cdr_s_sum_imgt
			global framework_r_sum_kabat
			global framework_s_sum_kabat
			global cdr_r_sum_kabat
			global cdr_s_sum_kabat

			for mode in mode_list:
				for region in getVRegionsList(False):
					#print "Got a region ",region,"\n\n"
					r_key=region+" AA subst. ("+mode+")"
					s_key=region+" codons with silent mut. ("+mode+")"
					if((r_key in annMap) and (s_key in annMap)):
						if(region.startswith("F")):
							#is a framework
							if(mode=="imgt"):
								framework_r_sum_imgt+=annMap[r_key]
								framework_s_sum_imgt+=annMap[s_key]
							else:
								framework_r_sum_kabat+=annMap[r_key]
								framework_s_sum_kabat+=annMap[s_key]
						elif(region.startswith("C")):
							#is a CDR
							if(mode=="imgt"):
								cdr_r_sum_imgt+=annMap[r_key]
								cdr_s_sum_imgt+=annMap[s_key]
							else:
								cdr_r_sum_kabat+=annMap[r_key]
								cdr_s_sum_kabat+=annMap[s_key]
					b_num_key=region+" base substitutions ("+mode+")"
					b_dnm_key=region+" Number base-to-base aligned ("+mode+")"
					if(b_dnm_key in annMap and b_num_key in annMap):
						#make sure the keys are there
						num=annMap[b_num_key]
						dnm=annMap[b_dnm_key]
						if(num!=None and dnm!=None):
							if(region.startswith("F")):
								#is a framework
								if(mode=="imgt"):
									framework_b2b_aln_imgt+=dnm
									framework_b2b_sbst_imgt+=num
								else:
									framework_b2b_aln_kabat+=dnm
									framework_b2b_sbst_kabat+=num									
							else:
								#is a CDR
								if(mode=="imgt"):
									cdr_b2b_aln_imgt+=dnm
									cdr_b2b_subst_imgt+=num
								else:
									cdr_b2b_aln_kabat+=dnm
									cdr_b2b_subst_kabat+=num	

			

	#print "SHOWING MAP WITH OVER REGION TOTALS"
	#printMap(annMap)
	#print "\n\n\n"	

	
	annMap["Read identifier"]=read_rec.id


	t_map=getTypeNameScoreMap(read_result_obj,meta)

	alt_map=getVDJAltMap(t_map,topVDJ)
	segments=get_segment_list()
	for segment in segments:
		alt_key="Alternate "+segment+" gene"
		alt_list=alt_map[segment]
		if(len(alt_list)>0):
			alt_str=",".join(alt_list)
			annMap[alt_key]=alt_str
		else:
			annMap[alt_key]=None
			

	return annMap



def computeSmartRatio(num,denom):
	if(denom!=None):
		if(num!=None):
			s_num=float(num)
			s_dnm=float(denom)
			if(s_dnm!=0):
				ratio=s_num/s_dnm
				return ratio
			else:
				return None
		else:
			return None
	else:
		return None






def prodRearrangmentVJ(vInfo,jInfo,imgtdb_obj,organism):
	if(vInfo==None or jInfo==None):
		return  False
	query_max_v=vInfo['q. end']
	query_min_j=jInfo['q. start']
	if(query_max_v>=query_min_j):
		#junction ovrelap or J comes before V
		query_max_j=jInfo['q. end']
		query_min_v=vInfo['q. start']
		if(query_max_j<=query_min_v):
			#J aligns before V! Nooooooooo!
			return False
		else:
			pass	
	else:
		v_end_aln=vInfo['s. end']
		v_end_frm=getTheFrameForThisReferenceAtThisPosition(vInfo['subject ids'],organism,imgtdb_obj,v_end_aln)
		sys.exit(0)
		
	







def appendAnnToFileWithMap(fHandle,m,rid,read_name,desiredKeys=None,defaultValue="None",logHandle=None):


	#apapend intro data
	keys_to_append=[

	#PART 1 : basic alignment/germline data
	"Read identifier",
	"V gene",
	"J gene",
	"D gene",
	"V Sequence similarity",

	#PART 2 : post-alignment filter reports for standard use
	"Out-of-frame junction",
	"Missing CYS",
	"Missing TRP/PHE",
	"Stop Codon?",
	"Indels Found",
	"Only Frame-Preserving Indels Found"
	]	
	
	#PART 3 IMGT/KABAT
	#append IMGT/KABAT region data
	modes=get_domain_modes()
	modes.sort()
	for mode in modes:
		keys_to_append.append("CDR3 AA ("+mode+")")
		#keys_to_append.append("CDR3 AA length ("+mode+")")
		keys_to_append.append("CDR3 NA ("+mode+")")
		#keys_to_append.append("CDR3 NA length ("+mode+")")
		#Base substition totals
		#keys_to_append.append("Total base substitution freq(%) over regions ("+mode+")") 	#a/b
		#keys_to_append.append("Total base substitutions over regions ("+mode+")") 		#a
		#keys_to_append.append("Total germline length (bp) over regions ("+mode+")") 		#b
		#keys_to_append.append("Total R:S ratio over all regions ("+mode+")") 			#c/e
		#keys_to_append.append("Total replacement (codon) freq(%) over regions ("+mode+")") 	#c/d
		#keys_to_append.append("Total replacements (codon) over regions ("+mode+")")		#c
		##keys_to_append.append("Total codons over regions ("+mode+")")				#d,f
		#keys_to_append.append("Total silent (codon) freq(%) over regions ("+mode+")")		#e/f
		#keys_to_append.append("Total silent (codon) over regions ("+mode+")")			#e		
		regions=getVRegionsList()
		for region in regions:
			keys_to_append.append(region.upper()+" aligned bases ("+mode+")")
			keys_to_append.append(region.upper()+" base subst. ("+mode+")")
			keys_to_append.append(region.upper()+" AA subst. ("+mode+")")
			keys_to_append.append(region.upper()+" codons with silent mut. ("+mode+")")
			#keys_to_append.append(region+" base substitution freq% ("+mode+")")
			#keys_to_append.append(region+" length ("+mode+")")
			#keys_to_append.append(region+" base substitutions ("+mode+")")
			#keys_to_append.append(region+" R:S ratio ("+mode+")")
			#keys_to_append.append(region+" replacement mutation freq% ("+mode+")")
			#keys_to_append.append(region+" replacement mutations (codons) ("+mode+")")
			#keys_to_append.append(region+" silent mutation freq% ("+mode+")")
			#keys_to_append.append(region+" silent mutations (codons) ("+mode+")")
			#keys_to_append.append(region+" indel frequency ("+mode+")")
			keys_to_append.append(region+" insertion count ("+mode+")")
			keys_to_append.append(region+" deletion count ("+mode+")")



	
	#keys_to_append.append("Homology% (over V and J)")
	#keys_to_append.append("Indel frequency (over V)")
	#keys_to_append.append("Insertion count (over V)")
	#keys_to_append.append("Length (over V and J)")
	#keys_to_append.append("Mutations (over V and J)")
	#keys_to_append.append("Mutations (over V)")
	#keys_to_append.append("Nonsynonymous base substitutions (over V and J)")
	#keys_to_append.append("Nonsynonymous base substitutions (over V)")
	#keys_to_append.append("AA (over V and J)")
	keys_to_append.append("AA (over V)")
	#keys_to_append.append("Nucleotide read (over V and J)")
	keys_to_append.append("Nucleotide read (over V)")
	#keys_to_append.append("R:S ratio (over V and J)")
	#keys_to_append.append("R:S ratio (over V)")
	#keys_to_append.append("Replacement mutation freq% (over V and J)")
	#keys_to_append.append("Replacement mutations (codons) (over V and J)")
	#keys_to_append.append("Base substitutions (over V and J)")
	#keys_to_append.append("CDR3 Stop codon? (over V and J)")
	#keys_to_append.append("Deletion count (over V)")
	#keys_to_append.append("Silent mutation freq% (over V and J)")
	#keys_to_append.append("Silent mutations (codons) (over V and J)")
	#keys_to_append.append("Base substitution freq% (over V and J)")
	keys_to_append.append("Stop codons? (over V)")
	#keys_to_append.append("Synonymous base substitutions (over V and J)")
	#keys_to_append.append("Synonymous base substitutions (over V)")
	segments=get_segment_list()
	for segment in sorted(segments,reverse=True):
		alt_key="Alternate "+segment+" gene"
		keys_to_append.append(alt_key)

	#keys_to_append.append("AGS Q Note")
	keys_to_append.append("Codon Muts (kabat)")
	keys_to_append.append("Amino Muts (kabat)")
	keys_to_append.append("Silent Codon Muts (kabat)")
	keys_to_append.append("Silent Amino Muts (kabat)")
	#keys_to_append.append("Release Version Tag")
	#keys_to_append.append("Release Version Hash")






###############################33
#extra stuff "reject section"
#	"FR1 Stop codons? (imgt)",
#	"FR1 homology% (imgt)",
#	"FR1 nucleotide read (imgt)"

	#use the keys in the to_append list and add them as keys
	keys=list()
	for ki in range(len(keys_to_append)):
		keys.append(keys_to_append[ki])
	#set the READ ID in the map
	m[keys[0]]=read_name
	global release_info_tag
	global release_info_hash	
	#m[keys[len(keys)-2]]=release_info_tag
	#m[keys[len(keys)-1]]=release_info_hash




	#now, append everything from the map passed in that hasn't been added
	#for map_key in m:
	#	if(not(map_key in keys)):
	#		keys.append(map_key)


	

	if(rid==1):
		#if read id is 1, write out column headers
		for k in range(len(keys)):
			if(k<len(keys)-1):
				fHandle.write(keys[k]+"\t")
			else:
				fHandle.write(keys[k])
		fHandle.write("\n")
	for k in range(len(keys)):
		if(k<len(keys)-1):
			sep="\t"
		else:
			sep=""
		if(keys[k] in m):
			if( ( keys[k].find("%")!=(-1)  or keys[k].upper().find("FREQ")!=(-1))    and type(m[keys[k]])==float):
				#if the column header contains a % and the type is float, multiply by 100
				fHandle.write(str( m[keys[k]]*100.00   )+sep)
			else:
				fHandle.write(str(m[keys[k]])+sep)
		else:
			fHandle.write(defaultValue+sep)
	fHandle.write("\n")
	if(not(logHandle==None)):
		n=0
		logHandle.write("\n\n\n********************************************\n")
		mkeys=m.keys()
		mkeys.sort()
		for k in mkeys:
			#logHandle.write(str(n)+"\t"+str(k)+"\t"+str(m[k])+"\n")	
			logHandle.write(str(m[k])+"\t"+str(k)+"\t"+str(n)+"\n")
			n+=1
		
		



def generateSampleLevelStatsJSON():

	#IMGT_BASED		
	imgt_json="{\n"
	imgt_json+="F_R_SUM : "+jsonOptWrapVal(framework_r_sum_imgt)+",\n"
	imgt_json+="F_S_SUM : "+jsonOptWrapVal(framework_s_sum_imgt)+",\n"
	imgt_json+="C_R_SUM : "+jsonOptWrapVal(cdr_r_sum_imgt)+",\n"
	imgt_json+="C_S_SUM : "+jsonOptWrapVal(cdr_s_sum_imgt)+",\n"
	imgt_json+="F_RS : "+jsonOptWrapVal(computeSmartRatio(framework_r_sum_imgt,framework_s_sum_imgt))+",\n"
	imgt_json+="C_RS : "+jsonOptWrapVal(computeSmartRatio(cdr_r_sum_imgt,cdr_s_sum_imgt))+",\n"
	imgt_json+="F_B2B : "+jsonOptWrapVal(framework_b2b_aln_imgt)+",\n"
	imgt_json+="F_BSB : "+jsonOptWrapVal(framework_b2b_sbst_imgt)+",\n"
	imgt_json+="F_BSB_F : "+jsonOptWrapVal(computeSmartRatio(framework_b2b_sbst_imgt,framework_b2b_aln_imgt))+",\n"
	imgt_json+="C_B2B : "+jsonOptWrapVal(cdr_b2b_aln_imgt)+",\n"
	imgt_json+="C_BSB : "+jsonOptWrapVal(cdr_b2b_subst_imgt)+",\n"
	imgt_json+="C_BSB_F : "+jsonOptWrapVal(computeSmartRatio(cdr_b2b_subst_imgt,cdr_b2b_aln_imgt))+"\n"
	imgt_json+="}"

	#KABAT_BASED
	kabat_json="{\n";
	kabat_json+="F_R_SUM : "+jsonOptWrapVal(framework_r_sum_kabat)+",\n"
	kabat_json+="F_S_SUM : "+jsonOptWrapVal(framework_s_sum_kabat)+",\n"
	kabat_json+="C_R_SUM : "+jsonOptWrapVal(cdr_r_sum_kabat)+",\n"
	kabat_json+="C_S_SUM : "+jsonOptWrapVal(cdr_s_sum_kabat)+",\n"
	kabat_json+="F_RS : "+jsonOptWrapVal(computeSmartRatio(framework_r_sum_kabat,framework_s_sum_kabat))+",\n"
	kabat_json+="C_RS : "+jsonOptWrapVal(computeSmartRatio(cdr_r_sum_kabat,cdr_s_sum_kabat))+",\n"
	kabat_json+="F_B2B : "+jsonOptWrapVal(framework_b2b_aln_kabat)+",\n"
	kabat_json+="F_BSB : "+jsonOptWrapVal(framework_b2b_sbst_kabat)+",\n"
	kabat_json+="F_BSB_F : "+jsonOptWrapVal(computeSmartRatio(framework_b2b_sbst_kabat,framework_b2b_aln_kabat))+",\n"
	kabat_json+="C_B2B : "+jsonOptWrapVal(cdr_b2b_aln_kabat)+",\n"
	kabat_json+="C_BSB : "+jsonOptWrapVal(cdr_b2b_subst_kabat)+",\n"
	kabat_json+="C_BSB_F : "+jsonOptWrapVal(computeSmartRatio(cdr_b2b_subst_kabat,cdr_b2b_aln_kabat))+"\n"
	kabat_json+="}"
	
	#WHOLE json
	whole_json="{\n"
	whole_json+="NUM_SEQS : "+jsonOptWrapVal(num_seqs_pf)+",\n"
	whole_json+="WHOLE_BSB : "+jsonOptWrapVal(whole_read_b2b_sbst)+",\n"
	whole_json+="WHOLE_B2B : "+jsonOptWrapVal(whole_read_b2b_aln)+",\n"
	whole_json+="WHOLE_BSB_F : "+jsonOptWrapVal(computeSmartRatio(whole_read_b2b_sbst,whole_read_b2b_aln))+",\n"
	whole_json+="IMGT : "+imgt_json+",\n"
	whole_json+="KABAT : "+kabat_json+"\n"
	whole_json+="}\n"

	return whole_json
		
def summary_file_header_mappings(header):
        """Generate mapping from header names to columns given the header line for a summary file"""
        headerDict = {}
        headers = header.split('\t')
        # TODO: this should use an established set of names
        headerDict = {h: headers.index(h) for h in headers}
        return headerDict