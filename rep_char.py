#!/usr/bin/env python

import vdjml
from igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs,rev_comp_dna
#from utils import printMap,get_domain_modes,biopythonTranslate,makeAllMapValuesVal,getNumberBpInAlnStr,repStr
from utils import *
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
#from vdjml_utils import getTopVDJItems,getHitInfo,getVDJServerRegionAlignmentFromLargerVAlignmentPyObj,getAlignmentString,getTypeNameScoreMap,getVDJTieMap
from vdjml_utils import *
from Bio import SeqIO
from segment_utils import IncrementMapWrapper,getVRegionsList,getEmptyRegCharMap,getAdjustedCDR3StartFromRefDirSetAllele,getTheFrameForThisReferenceAtThisPosition,getVRegionStartAndStopGivenRefData,getADJCDR3EndFromJAllele,alleleIsTR
from segment_utils import recombFreqManager
from char_utils import getNumberBaseSubsFromBTOP,getNumberIndelsFromBTOP,getIndelMapFromBTOP
from alignment import alignment,CodonAnalysis,codonAnalyzer
from CharacterizationThread import CharacterizationThread
from char_utils import getRegPosFromInvertedPos
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
	framePreservingIndelsFlag=True
	#printMap(vInfo)
	if(vInfo!=None):
		v_q_seq=vInfo['query seq']
		v_s_seq=vInfo['subject seq']
		strs_to_analyze.append(v_q_seq)
		strs_to_analyze.append(v_s_seq)
	if(dInfo!=None):
		d_q_seq=vInfo['query seq']
		d_s_seq=vInfo['subject seq']
		strs_to_analyze.append(d_q_seq)
		strs_to_analyze.append(d_s_seq)
	if(jInfo!=None):
		j_q_seq=vInfo['query seq']
		j_s_seq=vInfo['subject seq']
		strs_to_analyze.append(j_q_seq)
		strs_to_analyze.append(j_s_seq)
	gap_re=re.compile(r'([\-]{1,})')	
	for s in strs_to_analyze:
		if(s.find("-")!=(-1)):
			#found a gap!
			basicIndelsFlags=True
		gap_groups=re.findall(gap_re,s)
		for gap_group in gap_groups:
			gap_len=len(gap_group)
			if(gap_len%3==0):
				#it's frame preserving!
				pass
			else:
				#found at least one non-frame preserving gap!
				framePreservingIndelsFlag=False
	gap_ret_package=[basicIndelsFlags,framePreservingIndelsFlag]
	return gap_ret_package








def getValnWholeSeqStopFlag(vInfo,dInfo,jInfo,imgtdb_obj,organism,annMap,seq_rec):
	if(vInfo==None):
		#if no V alignment how can there be a frame, a J alignment, or a stop codon?
		return None
	else:
		#vAlnObj=alignment(vInfo['query seq'],vInfo['subject seq'],vInfo['q. start'],vInfo['q. end'],vInfo['s. start'],vInfo['s. end'])
		#vAlnNoGapStart=vAlnObj.getGEQAlignmentFirstNonGap()
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
		q_end=max(int(info_to_use['q. end']),int(vInfo['q. end']))
		if(not(vInfo['is_inverted'])):
			to_trans=str(seq_rec.seq)[(q_start-1):q_end].upper()
		else:
			proper=rev_comp_dna(str(seq_rec.seq))
			to_trans=str(proper[q_start-1:q_end]).upper()
		#print "s_start_frame is ",s_start_frame
		if(s_start_frame!=0):
			to_trans=to_trans[(3-s_start_frame):]
		#print "To be translated for read=",str(seq_rec.id)," is ",to_trans
		to_trans_len=len(to_trans)
		if((to_trans_len%3)!=0):
			necessary_trim=(to_trans_len%3)
			trim_to_trans=to_trans[:to_trans_len-necessary_trim]
		else:
			trim_to_trans=to_trans
		#print "The trim_to_trans is ",trim_to_trans
		translation=codonAnalyzer.fastTransStr(trim_to_trans)
		stopFlag=stringContainsAStar(translation)
		#print "The translation : ",translation
		#print "\n\n\n\n"
		return stopFlag
		#return False
		
		




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
def returnVAlnAllCharMap(vInfo,imgtdb_obj,read_result_obj=None,meta=None):
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
				new_map['Sequence similarity']=vCharMap[k]
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
		V_map=returnVAlnAllCharMap(vInfo,imgtdb_obj)
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
	mut_map_info=annotationMutationMap(vInfo,dInfo,jInfo,alignment_output_queue,num_submitted_jobs,imgtdb_obj,myCodonCounter,organism,read_rec,cdr3_map,0.85)
	qualifNote=mut_map_info[0]
	maps=mut_map_info[1]
	amino_map=maps['aminos']
	codon_map=maps['codons']
	amino_silent_map=maps['aminos_silent']
	codon_silent_map=maps['codons_silent']
	annMap["AGS Q Note"]=qualifNote
	annMap["AGS Codon Muts"]=codon_map
	annMap["AGS Amino Muts"]=amino_map
	annMap["AGS Silent Codon Muts"]=codon_silent_map
	annMap["AGS Silent Amino Muts"]=amino_silent_map


	#here, compute the region total BSB, REPLACEMENTS, SILENTS
	#retrieve the mode/region results from the annMap
	for mode in mode_list:
		bsb_reg_tot=0
		rplcmt_reg_tot=0
		silent_reg_tot=0
		ref_reg_len_tot=0
		ref_reg_len_cdn_tot=0
		temp_len_codons=0
		for region in getVRegionsList():
			bsb_key=region+" base substitutions ("+mode+")"
			ref_reg_srt_key=global_key_base+mode+"_"+region+"_ref_srt"
			ref_reg_end_key=global_key_base+mode+"_"+region+"_ref_end"
			rplcmt_key=region+" replacement mutations (codons) ("+mode+")"
			silent_key=region+" silent mutations (codons) ("+mode+")"
			#replacement/silent ref germline increment
			if((ref_reg_srt_key in annMap) and (ref_reg_end_key in annMap)):
				#base pairs length
				this_region_na_len=int(annMap[ref_reg_end_key])-int(annMap[ref_reg_srt_key])+1
				ref_reg_len_tot=ref_reg_len_tot+(this_region_na_len)
				#codon length (first na, then divide by 3, then increment total)
				temp_len_codons=int(math.floor(float(this_region_na_len)/float(3.0)))
				#print "temp_len_codons=",temp_len_codons,"mode=",mode,"region=",region,"read=",read_rec.id
				#print "this_region_na_len=",this_region_na_len
				ref_reg_len_cdn_tot+=temp_len_codons
				#print "ref_reg_len_cdn_tot=",ref_reg_len_cdn_tot,"mode=",mode
			else:
				#print "at least one of "+ref_reg_srt_key+" or "+ref_reg_end_key+" missing in "+read_rec.id
				pass
			#BSB totals incrementing
			if(bsb_key in annMap):
				bsb_reg_tot+=annMap[bsb_key]
			else:
				#print bsb_key+" not found for "+read_rec.id
				pass
			#replacement totals incrementing
			if(rplcmt_key in annMap):
				rplcmt_reg_tot+=annMap[rplcmt_key]
			else:
				#print rplcmt_key+" not found "+read_rec.id
				pass
			#silent totals incrementing
			if(silent_key in annMap):
				silent_reg_tot+=annMap[silent_key]
			else:
				#print silent_key+" not found for "+read_rec.id
				pass
				
		#having accumulated totals for the regions, compute mode values
		#Base substitutions
		sum_bsb_overregion_key_num="Total base substitutions over regions ("+mode+")"
		sum_bsb_overregion_key_dnm="Total germline length (bp) over regions ("+mode+")"
		annMap[sum_bsb_overregion_key_dnm]=ref_reg_len_tot
		annMap[sum_bsb_overregion_key_num]=bsb_reg_tot
		if(ref_reg_len_tot!=0):
			
			reg_bsb_tot=float(bsb_reg_tot)/float(ref_reg_len_tot)
		else:
			reg_bsb_tot=None
		bsb_overregion_tot_key="Total base substitution freq(%) over regions ("+mode+")"
		annMap[bsb_overregion_tot_key]=reg_bsb_tot
		#Codon/AA replacements (and silents)
		sum_repl_overregion_key_num="Total replacements (codon) over regions ("+mode+")"
		sum_slnt_overregion_key_num="Total silent (codon) over regions ("+mode+")"
		sum_gl_overregion_key_dnm="Total codons over regions ("+mode+")"
		annMap[sum_repl_overregion_key_num]=rplcmt_reg_tot
		annMap[sum_slnt_overregion_key_num]=silent_reg_tot
		annMap[sum_gl_overregion_key_dnm]=ref_reg_len_cdn_tot
		#print "ref_reg_len_cdn_tot=",ref_reg_len_cdn_tot
		if(ref_reg_len_cdn_tot!=0):
			rep_over_region=float(rplcmt_reg_tot)/float(ref_reg_len_cdn_tot)
			silent_over_region=float(silent_reg_tot)/float(ref_reg_len_cdn_tot)
		else:
			rep_over_region=None
			silent_over_region=None
		rep_overregion_tot_key="Total replacement (codon) freq(%) over regions ("+mode+")"
		annMap[rep_overregion_tot_key]=rep_over_region
		silent_overregion_tot_key="Total silent (codon) freq(%) over regions ("+mode+")"
		annMap[silent_overregion_tot_key]=silent_over_region
		total_rs_key="Total R:S ratio over all regions ("+mode+")"
		if(annMap[sum_slnt_overregion_key_num]!=0):
			#print "TOTAL R: RATIO HERE"
			#print "NUM=",float(annMap[sum_repl_overregion_key_num]),"key=",sum_repl_overregion_key_num
			#print "DNM=",float(annMap[sum_slnt_overregion_key_num]),"key=",sum_slnt_overregion_key_num
			#print "RTO=",float(float(annMap[sum_repl_overregion_key_num])/float(annMap[sum_slnt_overregion_key_num])),"key=",total_rs_key
			annMap[total_rs_key]=float(float(annMap[sum_repl_overregion_key_num])/float(annMap[sum_slnt_overregion_key_num]))
		else:
			annMap[total_rs_key]=None
			

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
		v_end_frm=getTheFrameForThisReferenceAtThisPosition(vInfor['subject ids'],organism,imgtdb_obj,v_end_aln)
		sys.exit(0)
		
	







def appendAnnToFileWithMap(fHandl,m,rid,read_name,desiredKeys=None,defaultValue="None",logHandle=None):


	#apapend intro data
	keys_to_append=[

	#PART 1 : basic alignment/germline data
	"Read identifier",
	"V gene",
	"J gene",
	"D gene",
	"Sequence similarity",

	#PART 2 : post-alignment filter reports for standard use
	"Out-of-frame junction",
	"Missing CYS",
	"Missing TRP/PHE",
	"Stop Codon?",
	"Indels Found",
	"Only Frame-Preserving Indels Found"
	]


	#"Stop codons? (over V and J)",
	#"Productive CDR3 rearrangement (T/F)",
	#"Base substitution freq% (over V)",
	#"Base substitutions (over V)",
	#"Length (over V)",
	#"Replacement mutation freq% (over V)",
	#"Replacement mutations (codons) (over V)",
	#"Silent mutation freq% (over V)",
	#"Silent mutations (codons) (over V)",
	#"Indel frequency (over V and J)",
	#"Insertion count (over V and J)",
	#"Deletion count (over V and J)"]
	

	
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
			#keys_to_append.append(region+" insertion count ("+mode+")")
			#keys_to_append.append(region+" deletion count ("+mode+")")



	
	#keys_to_append.append("Homology% (over V and J)")
	#keys_to_append.append("Indel frequency (over V)")
	#keys_to_append.append("Insertion count (over V)")
	#keys_to_append.append("Length (over V and J)")
	#keys_to_append.append("Mutations (over V and J)")
	#keys_to_append.append("Mutations (over V)")
	#keys_to_append.append("Nonsynonymous base substitutions (over V and J)")
	#keys_to_append.append("Nonsynonymous base substitutions (over V)")
	#keys_to_append.append("AA (over V and J)")
	#keys_to_append.append("AA (over V)")
	#keys_to_append.append("Nucleotide read (over V and J)")
	#keys_to_append.append("Nucleotide read (over V)")
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
	#keys_to_append.append("Stop codons? (over V)")
	#keys_to_append.append("Synonymous base substitutions (over V and J)")
	#keys_to_append.append("Synonymous base substitutions (over V)")
	segments=get_segment_list()
	for segment in sorted(segments,reverse=True):
		alt_key="Alternate "+segment+" gene"
		keys_to_append.append(alt_key)

	#keys_to_append.append("AGS Q Note")
	#keys_to_append.append("AGS Codon Muts")
	#keys_to_append.append("AGS Amino Muts")
	#keys_to_append.append("AGS Silent Codon Muts")
	#keys_to_append.append("AGS Silent Amino Muts")
	keys_to_append.append("Release Version Tag")
	keys_to_append.append("Release Version Hash")






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
	m[keys[len(keys)-2]]=release_info_tag
	m[keys[len(keys)-1]]=release_info_hash




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
		
		


#add on the kinds of arguments this program accepts
def add_rep_char_args_to_parser(parser):
	version_info=getVersionInfo()
	parser.description+=' Generate read-level repertoire-characterization data. VERSION '+version_info
	parser.add_argument('-json_out',type=str,nargs=1,default="/dev/stdout",help="output file for the JSON segment count IMGT hierarchy")
	parser.add_argument('-combo_out',type=str,nargs=1,default="/dev/stdout",help="output file for JSON segment recombination count data")
	parser.add_argument('-cdr3_hist_out',type=str,nargs=1,default="/dev/stdout",help="output file for the CDR3 histogram of lengths (both kabat and imgt systems)")
	parser.add_argument('-skip_char',action='store_true', default=False,help="If this is set, then characterization besides germline segment assignment and CDR3 length is skipped")
	parser.add_argument('char_out',type=str,nargs=1,help="the path to the output TSV file of read-level repertoire characterization data")
	parser.add_argument('vdj_db_root',type=str,nargs=1,help="path to the VDJ directory root REQUIRED")
	parser.add_argument('qry_fasta',type=str,nargs=1,help="path to the input fasta file of query (Rep-Seq) data input to IgBLAST")
	#parser.add_argument('-region_out'
	parser.add_argument('db_organism',type=str,nargs=1,default="human",help="the organism IgBLASTed against;must exist under vdj_db_root",choices=["human","Mus_musculus"])
	return parser


		


if (__name__=="__main__"):
	parser=makeParserArgs()
	parser=add_rep_char_args_to_parser(parser)
	args = parser.parse_args()
	if(args):
		organism=extractAsItemOrFirstFromList(args.db_organism)
		mf=makeVDJMLDefaultMetaAndFactoryFromArgs(args)
		meta=mf[0]
		fact=mf[1]
		imgtdb_obj=imgt_db(extractAsItemOrFirstFromList(args.vdj_db_root))
		#print "a fact is ",fact
		#print "the file is ",args.igblast_in[0]
		query_fasta=extractAsItemOrFirstFromList(args.qry_fasta)
		if(os.path.exists(query_fasta) and os.path.isfile(query_fasta)):
			#tot_reads=countFastaReads(query_fasta)
			tot_reads=None
			#do we really wanna count the reads?????
		else:
			tot_reads=None
		fasta_reader=SeqIO.parse(open(query_fasta, "r"), "fasta")
		modes=['kabat','imgt']
		my_cdr3_map=histoMapClass(modes)
		cdr3_hist_out_file=args.cdr3_hist_out
		segment_counter=IncrementMapWrapper()
		combo_counter=recombFreqManager()
		read_num=1
		segments_json_out=extractAsItemOrFirstFromList(args.json_out)
		#print "to write vdjml to ",args.vdjml_out
		rrw = vdjml.Vdjml_writer(extractAsItemOrFirstFromList(args.vdjml_out), meta)
		rep_char_out=extractAsItemOrFirstFromList(args.char_out)
		fHandle=open(rep_char_out,'w')
		if(False):
			logFile=rep_char_out+".log"
			logHandle=open(logFile,'w')
		else:
			logHandle=None
		if(tot_reads==None):
			every_read=100
		else:
			every_read=math.log(tot_reads,10.0)
			if(every_read<=1):
				every_read=1
			else:
				every_read=int(every_read)
				every_read=int(10**every_read)
			#print "every_read at ",every_read

		#test release info
		loadReleaseInfoTagAndHash()
		if(release_info_tag==None or release_info_hash==None):
			print "ERROR IN RELEASE_INFO ?!"



		for read_result_obj in scanOutputToVDJML(extractAsItemOrFirstFromList(args.igblast_in),fact,query_fasta):


			#prepare for the iteration and give a possible status message...
			if(read_num>1 and read_num%every_read==0):
				print "Processing read",read_num,"..."
			elif(read_num==1):
				print "Processing reads..."
			query_record=fasta_reader.next()

			#analyze a read's results
			read_analysis_results=rep_char_read(read_result_obj,meta,organism,imgtdb_obj,query_record,args.skip_char)

			#handle cdr3 length/histogram
			if('cdr3_length_results' in read_analysis_results):
				cdr3_res=read_analysis_results['cdr3_length_results']
				for mode in modes:
					my_cdr3_map.inc(mode,cdr3_res[mode])

			#handle segment counting
			if('VDJ' in read_analysis_results):
				segments=read_analysis_results['VDJ']
				for s in segments:
					actual=segments[s]
					if(actual is not None):
						segment_counter.increment(actual)
				vseg=segments['V']
				dseg=segments['D']
				jseg=segments['J']
				combo_counter.addVDJRecombination(vseg,dseg,jseg,False)

			#process for writing
			rrw(read_result_obj)

			#write rep-char
			appendAnnToFileWithMap(fHandle,read_analysis_results['ann_map'],read_num,query_record.id,None,"None",logHandle)

			#increment the read number
			read_num+=1

		print "Processed total of",str(read_num-1),"reads!"

		#close rep_char out
		fHandle.close()

		#close log
		if(logHandle is not None):
			logHandle.close()

		#write the CDR3 hist when non-dev-null
		if(type(cdr3_hist_out_file)==list):
			cdr3_hist_out_file=extractAsItemOrFirstFromList(cdr3_hist_out_file)
		if(not(cdr3_hist_out_file=="/dev/null")):
			my_cdr3_map.writeToFile(cdr3_hist_out_file)
			print "Wrote CDR3 lengths histogram to ",cdr3_hist_out_file

		#write AGS/NMO information
		if(organism=="human"):
			print "AGS_SCORE\t"+str(myCodonCounter.computeAGS())+"\tAGS6_TOT_RM\t"+str(myCodonCounter.computeAGS6TotRM())+"\tTOT_RM\t"+str(myCodonCounter.computeSampTotRM())
			print "AGS5_SCRE\t"+str(myCodonCounter.computeAGS5())+"\tAGS5_TOT_RM\t"+str(myCodonCounter.computeAGS5TotRM())+"\tTOT_RM\t"+str(myCodonCounter.computeSampTotRM())
			print "NMO_SCORE\t"+str(myCodonCounter.computeNMO())+"\tNMO_QUERIES_RM\t"+str(myCodonCounter.queriesWithRM)+"\tNMO_SAMP_NUC_TOT\t"+str(myCodonCounter.computeSampNMOTot())
		else:
			print "No AGS/NMO data for non-human analysis!"


		#write the segment counts when non-dev-null
		if(type(segments_json_out)==list):
			segments_json_out=segments_json_out
		print "Writing JSON segment counts output to ",segments_json_out,"..."
		segment_counter.JSONIFYToFile(
			extractAsItemOrFirstFromList(args.vdj_db_root),
			organism,
			segments_json_out,
			False,
			imgtdb_obj.getPickleFullPath()
			)
		print "Writing JSON segment counts output complete!"
		recomb_out_file=extractAsItemOrFirstFromList(args.combo_out)
		print "Writing JSONS segment combination frequency data to ",recomb_out_file
		combo_counter.writeJSONToFile(recomb_out_file)
		




	else:
		#print "error in args!"
		parser.print_help()




