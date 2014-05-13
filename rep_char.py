#!/usr/bin/env python

import vdjml
from vdjml_igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs
#from utils import printMap,get_domain_modes,biopythonTranslate,makeAllMapValuesVal,getNumberBpInAlnStr,repStr
from utils import *
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
from vdjml_utils import getTopVDJItems,getRegionsObjsFromSegmentCombo,getHitInfo,getProductiveRearrangmentFlag,getVDJServerRegionAlignmentFromLargerVAlignmentPyObj,getAlignmentString,getTypeNameScoreMap,getVDJTieMap
from Bio import SeqIO
from segment_utils import IncrementMapWrapper,getVRegionsList,getRegionSpecifcCharacterization,getEmptyRegCharMap,getAdjustedCDR3StartFromRefDirSetAllele,getTheFrameForThisReferenceAtThisPosition,getCDR3RegionSpecificCharacterizationSubAln,getVRegionStartAndStopGivenRefData,getADJCDR3EndFromJAllele,alleleIsTR
from char_utils import getNumberBaseSubsFromBTOP,getNumberIndelsFromBTOP,getIndelMapFromBTOP
from alignment import alignment
from CharacterizationThread import CharacterizationThread
from vdjml_igblast_parse import rev_comp_dna,getRegPosFromInvertedPos
import re
import Queue
import threading
import time
import multiprocessing
import math


global_key_base="vdj_server_ann_"



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

	#printMap(read_ann_map,True)
	#sys.exit(0)

	return return_obj
	





#specialized annotation for CDR3 region!
def readAnnotate_cdr3(read_result_obj,meta,organism,imgtdb_obj,read_rec,read_ann_map,cdr3_length_results):
	global global_key_base
	for mode in get_domain_modes():
		f=cdr3_length_results[mode+'_from']
		t=cdr3_length_results[mode+'_to']
		len_key_aa="CDR3 AA length ("+mode+")"
		aa_key="CDR3 AA ("+mode+")"
		
		if(f!=(-1) and t!=(-1) and cdr3_length_results[mode]!=(-1) ):
			qw=str(read_rec.seq)
			if(cdr3_length_results['qry_rev']):
				qw=rev_comp_dna(qw)
			query_cdr3=qw[f-1:t]
			read_ann_map[global_key_base+mode+'_cdr3_na']=query_cdr3
			read_ann_map[aa_key]=biopythonTranslate(read_ann_map[global_key_base+mode+'_cdr3_na'])
			#add in lengths!
			read_ann_map[global_key_base+mode+'_cdr3_na_len']=len(query_cdr3)
			read_ann_map[len_key_aa]=len(read_ann_map[aa_key])
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




# get CDR3 spec char. but return q_pos of last_analyzed position and also limit to analysis greater than given char
def getCDR3SpecifcCharacterization(s_aln,q_aln,frame_mask,q_aln_start,limit):
	empty_map=getEmptyRegCharMap()
	#the only keys are zero!
	#print "s",s_aln
	#print "q",q_aln
	#print "f",frame_mask
	#print "x",q_aln_start
	#print "l",limit
	q_pos=q_aln_start
	temp=0
	#while(temp<len(s_aln)):
		
	
	
	

#add to the characterization map CDR3 region characterization
#each vmap,dmap,jmap are the alignment infos as returned from getHitInfo
def cdr3RegCharAnalysis(vMap,dMap,jMap,mode,cdr3_anal_map,organism,imgtdb_obj):
	empty_map=getEmptyRegCharMap()
	#first do input validation
	if(vMap==None or jMap==None):
		empty_map=makeAllMapValuesVal(empty_map,-1)
		return empty_map
	if(not(mode in cdr3_anal_map)):
		empty_map=makeAllMapValuesVal(empty_map,-1)
		return empty_map
	elif(cdr3_anal_map[mode]==(-1)):
		empty_map=makeAllMapValuesVal(empty_map,-1)
		return empty_map
	for m in [vMap,dMap,jMap]:
		if('subject seq'  not in m):
			empty_map=makeAllMapValuesVal(empty_map,-1)
			return empty_map
		elif('query seq' not in m):
			empty_map=makeAllMapValuesVal(empty_map,-1)
			return empty_map
	#accumulate counts for V
	max_v_anal=(-1)
	max_d_anal=(-1)
	q_cdr3_start=cdr3_anal_map[mode+'_from']
	q_cdr3_end=cdr3_anal_map[mode+'_to']
	v_s_aln=vMap['subject seq']
	v_q_aln=vMap['query seq']
	temp_v_pos=int(vMap['s. start'])
	temp_v=0
	temp_v_q_pos=int(vMap['q. start'])
	VrefName=vMap['subject ids']
	v_ref_cdr3_start=getAdjustedCDR3StartFromRefDirSetAllele(VrefName,imgtdb_obj,organism,mode)
	frame_mask=list()
	qry_cdr3_start_last_frame=(-1)
	cdr3_q_start_for_return=(-1)
	cdr3_s_aln=""
	cdr3_q_aln=""
	v_tot=0
	v_same=0
	#print "to enter while"
	while(temp_v<len(v_s_aln)):
		if(temp_v_pos>=v_ref_cdr3_start):
			cdr3_s_aln+=v_s_aln[temp_v]
			cdr3_q_aln+=v_q_aln[temp_v]
			qry_cdr3_start_last_frame=getTheFrameForThisReferenceAtThisPosition(VrefName,organism,imgtdb_obj,temp_v_pos)
			if(v_s_aln[temp_v]!="-" and v_q_aln[temp_v]!="-"):
				v_tot+=1
				if(v_s_aln[temp_v]==v_q_aln[temp_v]):
					v_same+=1
			frame_mask.append(qry_cdr3_start_last_frame)
		if(v_q_aln[temp_v]!="-" and cdr3_q_start_for_return!=(-1)):
			cdr3_q_start_for_return=temp_v_q_pos
		if(v_q_aln[temp_v]!="-"):
			temp_v_q_pos+=1
		if(v_s_aln[temp_v]!="-"):
			lastVPos=temp_v_pos
			temp_v_pos+=1
		temp_v+=1
	#cdr3_v_char_map=getRegionSpecifcCharacterization(cdr3_s_aln,cdr3_q_aln,"CDR3",frame_mask,dMode)
	#print "THE CDR3 ALN IN V is (q top, s bottom)"
	#print cdr3_q_aln
	#print frame_mask
	#print "stop here"
	#print "last q",temp_v_q_pos
	getCDR3SpecifcCharacterization(cdr3_s_aln,cdr3_q_aln,frame_mask,temp_v_q_pos,(-1))
	map_map=dict()
	super_s=makeEmptyArrayOfStringsOfLen(jMap['q. end'])
	#sys.exit(0)

	sys.exit(0)





def getValnWholeSeqStopFlag(vInfo,dInfo,jInfo,imgtdb_obj,organism,annMap,seq_rec):
	if(not(vInfo==None)):
		vAlnObj=alignment(vInfo['query seq'],vInfo['subject seq'],vInfo['q. start'],vInfo['q. end'],vInfo['s. start'],vInfo['s. end'])
		vAlnNoGapStart=vAlnObj.getGEQAlignmentFirstNonGap()
		s_start=vAlnNoGapStart.s_start
		s_start_frame=getTheFrameForThisReferenceAtThisPosition(vInfo['subject ids'],organism,imgtdb_obj,s_start)
		numNToAdd=s_start_frame%3
		Ns=repStr("N",numNToAdd)
		q_start=vAlnNoGapStart.q_start
		#print "Initial q_start is ",q_start
		#translate all the way to the end of the query here!!!!
		q_end=vAlnNoGapStart.q_end
		#if(vInfo['is_inverted']):
		#	q_start=getRegPosFromInvertedPos(q_start,len(seq_rec.seq))
		#	print "q_start after inversion : ",q_start
		#	q_end=getRegPosFromInvertedPos(q_end,len(seq_rec.seq))
		if(not(jInfo==None)):
			jAlnObj=alignment(jInfo['query seq'],jInfo['subject seq'],jInfo['q. start'],jInfo['q. end'],jInfo['s. start'],jInfo['s. end'])
			j_q_end=jAlnObj.q_end
			#print "Initial q_end from j",j_q_end			
			#if(jInfo['is_inverted']):
			#	j_q_end=getRegPosFromInvertedPos(j_q_end,len(seq_rec.seq))
			#print "Used q_end from j",j_q_end
			#the purpose of this "max" function is to let q_end stay where it is in the case the J aligns BEFORE V does!
			q_end=max(q_end,j_q_end)
		else:
			#no J alignment so try to get end of D alignment
			if(not(dInfo==None)):
				dAlnObj=alignment(dInfo['query seq'],dInfo['subject seq'],dInfo['q. start'],dInfo['q. end'],dInfo['s. start'],dInfo['s. end'])
				d_q_end=dAlnObj.q_end
				#print "Initial q_end from d",d_q_end
				q_end=max(q_end,d_q_end)
			else:
				#d is none
				pass
		if(vInfo['is_inverted']):
			query_to_exm=rev_comp_dna(str(seq_rec.seq))[q_start-1:q_end+1]
		else:
			query_to_exm=str(seq_rec.seq)[q_start-1:q_end+1]
		query_to_exm=Ns+query_to_exm
		trx_to_exm=biopythonTranslate(query_to_exm)
		#print "The query and translation to examine for ",seq_rec.id," (with q_start=",q_start," and q_end=",q_end,")  :"
		#print query_to_exm
		#print trx_to_exm
		if(trx_to_exm.find("*")!=(-1)):
			return True
		else:
			return False
	else:
		return None
		




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
			#if no info is found, return the map from just V
			#return preCDR3Aln.characterize()
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
def returnVAlnAllCharMap(vInfo,imgtdb_obj):
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
		new_map=dict()
		for k in vCharMap:
			newKey=k
			newKey=newKey[0].upper()+newKey[1:]
			new_map[newKey+ " (over V)"]=vCharMap[k]
		new_map['V length (nucleotides)']=abs(vInfo['q. start']-vInfo['q. end'])
		return new_map





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

	numN=postCDR3AlnObj.q_start-preCDR3Aln.q_end-1
	NSub=repStr("N",numN)
	NQry=repStr("N",numN)
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
		topkey=seg+" gene (highest score)"
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



	#characterization beside segments/CDR3
	global characterization_thread_set
	global characterization_queue
	global characterization_queue_results
	num_submitted_jobs=0
	if(not(skip_char)):
		#V seq characterization (OVER v) only
		V_map=returnVAlnAllCharMap(vInfo,imgtdb_obj)
		if(V_map is not None):
			for k in V_map:
				annMap[k]=V_map[k]

		#whole seq characterization (over V and J)
		whole_char_map=returnWholeSeqCharMap(vInfo,jInfo,imgtdb_obj,organism,annMap)
		for w in whole_char_map:
			#new_key=global_key_base+'whole_seq_'+w
			new_key=w+" (over V and J)"
			new_key=new_key[0].upper()+new_key[1:]
			annMap[new_key]=whole_char_map[w]
		whole_seq_stp_cdn_Tot_flag=getValnWholeSeqStopFlag(vInfo,dInfo,jInfo,imgtdb_obj,organism,annMap,read_rec)
		annMap['vdj_server_whole_vj_stp_cdn']=whole_seq_stp_cdn_Tot_flag


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
			#num_threads=10
			#print "USING NUM_THREADS=",num_threads
			characterization_thread_set=set()
			for c in range(num_threads):
				#print "ABOUT TO CALL CONST...."
				temp_thread=CharacterizationThread(characterization_queue,characterization_queue_results)
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
						#print "Now to analyze region ",region," in mode",mode
						#raInfo=getVDJServerRegionAlignmentFromLargerVAlignmentPyObj(read_result_obj,meta,organism,mode,region,imgtdb_obj,False,read_rec)
						#regionAlignment=getVDJServerRegionAlignmentFromLargerVAlignmentPyObj(read_result_obj,meta,organism,mode,region,imgtdb_obj,False,read_rec)
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
						characterization_queue.put(char_job_dict)
						num_submitted_jobs+=1
						#print "JUST SUBMITTED A CHAR_QUEUE JOB :",char_job_dict
						#if(regionAlignment!=None and not(noneSeg_flag)):
						#	reg_ann_show_msg=False
						#	#reg_ann_msg="characterization for region="+region+" mode="+mode+" for read="+read_rec.id
						#	reg_ann_msg=None
						#	reg_ann=regionAlignment.characterize(reg_ann_msg,reg_ann_show_msg)
						#	#if(mode=="imgt" and read_rec.id=="FR3_STOP" and region=="FR3"):
						#	#	print "got ann "
						#	#	printMap(reg_ann)
						#	#	print "The alignment nice : "
						#	#	print regionAlignment.getNiceString()
						#	#	sys.exit(0)
						#	for key in reg_ann:
						#		annMap[key_base+region+"_"+key]=reg_ann[key]
						#	annMap[key_base+region+'_qry_aln']=regionAlignment.q_aln
						#	annMap[key_base+region+'_qry_srt']=regionAlignment.q_start
						#	annMap[key_base+region+'_qry_end']=regionAlignment.q_end
						#	annMap[key_base+region+'_ref_aln']=regionAlignment.s_aln
						#	annMap[key_base+region+'_frm_msk']=regionAlignment.s_frame_mask
						#else:
						#	#print raInfo
						#	#reg_ann=getEmptyRegCharMap()
						#	#for key in reg_ann:
						#	#	annMap[key_base+region+"_"+key]=(-1)
						#	pass
				else:
					#if mode=kabat and type is TR, then skip the regions!
					pass
		else:
			#no annotation possible if V is empty!
			pass

	if(num_submitted_jobs>0):
		#print "CHAR_QUEUE pre-join...."
		characterization_queue.join()
		get_res=0
		while(get_res<num_submitted_jobs):
			temp_results=characterization_queue_results.get()
			if(temp_results is not None):
				for temp_key in temp_results:
					annMap[temp_key]=temp_results[temp_key]			
			get_res+=1
	
	annMap["Read identifier"]=read_rec.id
	#getAlignmentString(read_result_obj,meta,query_record,imgtdb_obj,organism)
	#analyze_combinations(read_result_obj,meta,organism,imgtdb_obj,read_rec,annMap)
	#printMap(annMap,True)
	t_map=getTypeNameScoreMap(read_result_obj,meta)
	tie_map=getVDJTieMap(t_map,topVDJ)
	segments=get_segment_list()
	for segment in segments:
		#key=global_key_base+segment+'_hitscorelist'
		key=segment+" gene hit/score list"
		if(segment in t_map):
			annMap[key]=t_map[segment]
		else:
			annMap[key]=None
		#key=global_key_base+segment+'_tie'
		key=segment+" gene tie"
		if(segment in tie_map):
			annMap[key]=tie_map[segment]
		else:
			annMap[key]="None"

	return annMap





def analyze_combinations(read_result_obj,meta,organism,imgtdb_obj,read_rec,read_ann_map):
	segment_combinations=read_result_obj.segment_combinations()
	for s in range(len(segment_combinations)):
		#print "LOOKING AT COMBINATION # ",str(int(s+1))," of ",str(len(segment_combinations))," FOR READ ID=",read_result_obj.id()
		segment_combination=segment_combinations[s]
		analyzeRegionObjsFromSegmentCombo(segment_combination,meta,organism,imgtdb_obj,read_rec,read_ann_map)




def analyzeRegionObjsFromSegmentCombo(segment_combo,meta,organism,imgtdb_obj,read_rec,read_ann_map):
	global global_key_base
	#print "got a segment combo : ",segment_combo
	#print "regions is ",segment_combo.regions
	#print "now calling....",
	regions=segment_combo.regions()
	#print "the type of regions is ",str(type(regions))
	for region in regions:
		#print "region???",region
		#print "range : ",region.range_
		#print region.range_.pos_0()
		rgn_start=region.range_.pos_1() 
		rgn_end=rgn_start+region.range_.length()-1
		#print "start/end : ",rgn_start," and ",rgn_end
		rgn_ns=int(str(region.num_system_))
		#print "rgn_ns",rgn_ns
		if(rgn_ns==1):
			rgn_ns="imgt"
		else:
			rgn_ns="kabat"
		#print "rgn_ns",rgn_ns
		rgn_nm=meta[region.region_]
		rgn_nm=re.sub(r'W','',rgn_nm)
		lookup_key_base=global_key_base+rgn_ns+"_"+rgn_nm+"_"
		lookup_qry_start=lookup_key_base+"qry_srt"
		lookup_qry_end=lookup_key_base+"qry_end"
		#print "using lookups ",lookup_qry_start," and ",lookup_qry_end
		if((lookup_qry_start in read_ann_map) and (lookup_qry_end in read_ann_map)):
			#print "found lookup!"
			igblast_start=int(rgn_start)
			igblast_end=int(rgn_end)
			vdj_start=read_ann_map[lookup_qry_start]
			vdj_end=read_ann_map[lookup_qry_end]
			strt_match=(vdj_start==igblast_start)
			end_match=(vdj_end==igblast_end)
			#print "start and end matches : ",strt_match," and ",end_match
			if(strt_match and end_match):
				#print "total match!"
				pass
			else:
				#print "at least one mismatch! igblast start/end : ",igblast_start," and ",igblast_end," vdj start/end     : ",vdj_start," and ",vdj_end," : query=",read_ann_map['read_name']," seg=",read_ann_map['top_V']," rgn=",rgn_nm," rrv : ",getVRegionStartAndStopGivenRefData(read_ann_map['top_V'],organism,imgtdb_obj,rgn_nm,rgn_ns)
				pass
		else:
			#print "no lookup found for ",lookup_qry_start," and ",lookup_qry_end
			pass



		









#using the PYVDJML, extract information from the IGBLAST-notated
#regions (FR1, CDR1, etc) and compare them with VDJ-server lookedup values
def vSegmentRegionVDJAnalyse(read_result_obj,meta,organism,imgtdb_obj,read_rec):
	topVDJ=getTopVDJItems(read_result_obj,meta)
	topV=topVDJ['V']
	segment_combinations=read_result_obj.segment_combinations()
	for s in range(len(segment_combinations)):
		segment_combination=segment_combinations[s]
		#print "Looking at combination #",str(s+1)
		getRegionsObjsFromSegmentCombo(segment_combination)
		#sc1 = vdjml.Segment_combination(scb)
		#scb.add_region(name=reg_name,read_range=reg_interval,metric=mm)
		#		scb.insert_region(vdjml.Num_system.kabat,vdjml.Gene_region_type.cdr2,reg_interval,mm)
		#sc1.insert_region(
		#                          vdjml.Num_system.imgt,
		#                          vdjml.Gene_region_type.fr1,
		#                          vdjml.Interval.first_last_1(1,54),
		#                          vdjml.Match_metrics(100, 54)
		#                          )
		#	return read_result_obj





def appendAnnToFileWithMap(fHandl,m,rid,desiredKeys=None,defaultValue="None",logHandle=None):
	keys=[

	"Read sequence #",
	"Read identifier",
	"V gene (highest score)",
	"J gene (highest score)",
	"D gene (highest score)",
	"CDR3 AA length (imgt)",
	"CDR3 AA (imgt)",
	"CDR3 AA length (kabat)",
	"Stop codons? (over V and J)",
	"Productive CDR3 rearrangement (T/F)",
	"Insertion count (over V and J)",
	"Deletion count (over V and J)",
	"Indel frequency (over V and J)",
	"Homology% (over V)",
	"Base substitutions (over V)",
	"Length (over V)",
	"Base substitution freq% (over V)",
	"Replacement mutations (codons) (over V)",
	"Silent mutations (codons) (over V)",
	"Replacement mutation freq% (over V)",
	"Silent mutation freq% (over V)",
	"FR1 AA (imgt)",
	"FR1 AA (kabat)",
	"FR1 R:S ratio (imgt)",
	"FR1 R:S ratio (kabat)",
	"FR1 Stop codons? (imgt)",
	"FR1 Stop codons? (kabat)",
	"FR1 base substitution freq% (imgt)",
	"FR1 base substitution freq% (kabat)",
	"FR1 base substitutions (imgt)",
	"FR1 base substitutions (kabat)",
	"FR1 deletion count (imgt)",
	"FR1 deletion count (kabat)",
	"FR1 homology% (imgt)",
	"FR1 homology% (kabat)",
	"FR1 indel frequency (imgt)",
	"FR1 indel frequency (kabat)",
	"FR1 insertion count (imgt)",
	"FR1 insertion count (kabat)",
	"FR1 length (imgt)",
	"FR1 length (kabat)",
	"FR1 nucleotide read (imgt)",
	"FR1 nucleotide read (kabat)",
	"FR1 replacement mutation freq% (imgt)",
	"FR1 replacement mutation freq% (kabat)",
	"FR1 replacement mutations (codons) (imgt)",
	"FR1 replacement mutations (codons) (kabat)",
	"FR1 silent mutation freq% (imgt)",
	"FR1 silent mutation freq% (kabat)",
	"FR1 silent mutations (codons) (imgt)",
	"FR1 silent mutations (codons) (kabat)",
	"CDR1 AA (imgt)",
	"CDR1 AA (kabat)",
	"CDR1 R:S ratio (imgt)",
	"CDR1 R:S ratio (kabat)",
	"CDR1 Stop codons? (imgt)",
	"CDR1 Stop codons? (kabat)",
	"CDR1 base substitution freq% (imgt)",
	"CDR1 base substitution freq% (kabat)",
	"CDR1 base substitutions (imgt)",
	"CDR1 base substitutions (kabat)",
	"CDR1 deletion count (imgt)",
	"CDR1 deletion count (kabat)",
	"CDR1 homology% (imgt)",
	"CDR1 homology% (kabat)",
	"CDR1 indel frequency (imgt)",
	"CDR1 indel frequency (kabat)",
	"CDR1 insertion count (imgt)",
	"CDR1 insertion count (kabat)",
	"CDR1 length (imgt)",
	"CDR1 length (kabat)",
	"CDR1 nucleotide read (imgt)",
	"CDR1 nucleotide read (kabat)",
	"CDR1 replacement mutation freq% (imgt)",
	"CDR1 replacement mutation freq% (kabat)",
	"CDR1 replacement mutations (codons) (imgt)",
	"CDR1 replacement mutations (codons) (kabat)",
	"CDR1 silent mutation freq% (imgt)",
	"CDR1 silent mutation freq% (kabat)",
	"CDR1 silent mutations (codons) (imgt)",
	"CDR1 silent mutations (codons) (kabat)",
	"FR2 AA (imgt)",
	"FR2 AA (kabat)",
	"FR2 R:S ratio (imgt)",
	"FR2 R:S ratio (kabat)",
	"FR2 Stop codons? (imgt)",
	"FR2 Stop codons? (kabat)",
	"FR2 base substitution freq% (imgt)",
	"FR2 base substitution freq% (kabat)",
	"FR2 base substitutions (imgt)",
	"FR2 base substitutions (kabat)",
	"FR2 deletion count (imgt)",
	"FR2 deletion count (kabat)",
	"FR2 homology% (imgt)",
	"FR2 homology% (kabat)",
	"FR2 indel frequency (imgt)",
	"FR2 indel frequency (kabat)",
	"FR2 insertion count (imgt)",
	"FR2 insertion count (kabat)",
	"FR2 length (imgt)",
	"FR2 length (kabat)",
	"FR2 nucleotide read (imgt)",
	"FR2 nucleotide read (kabat)",
	"FR2 replacement mutation freq% (imgt)",
	"FR2 replacement mutation freq% (kabat)",
	"FR2 replacement mutations (codons) (imgt)",
	"FR2 replacement mutations (codons) (kabat)",
	"FR2 silent mutation freq% (imgt)",
	"FR2 silent mutation freq% (kabat)",
	"FR2 silent mutations (codons) (imgt)",
	"FR2 silent mutations (codons) (kabat)",
	"CDR2 AA (imgt)",
	"CDR2 AA (kabat)",
	"CDR2 R:S ratio (imgt)",
	"CDR2 R:S ratio (kabat)",
	"CDR2 Stop codons? (imgt)",
	"CDR2 Stop codons? (kabat)",
	"CDR2 base substitution freq% (imgt)",
	"CDR2 base substitution freq% (kabat)",
	"CDR2 base substitutions (imgt)",
	"CDR2 base substitutions (kabat)",
	"CDR2 deletion count (imgt)",
	"CDR2 deletion count (kabat)",
	"CDR2 homology% (imgt)",
	"CDR2 homology% (kabat)",
	"CDR2 indel frequency (imgt)",
	"CDR2 indel frequency (kabat)",
	"CDR2 insertion count (imgt)",
	"CDR2 insertion count (kabat)",
	"CDR2 length (imgt)",
	"CDR2 length (kabat)",
	"CDR2 nucleotide read (imgt)",
	"CDR2 nucleotide read (kabat)",
	"CDR2 replacement mutation freq% (imgt)",
	"CDR2 replacement mutation freq% (kabat)",
	"CDR2 replacement mutations (codons) (imgt)",
	"CDR2 replacement mutations (codons) (kabat)",
	"CDR2 silent mutation freq% (imgt)",
	"CDR2 silent mutation freq% (kabat)",
	"CDR2 silent mutations (codons) (imgt)",
	"CDR2 silent mutations (codons) (kabat)",
	"FR3 AA (imgt)",
	"FR3 AA (kabat)",
	"FR3 R:S ratio (imgt)",
	"FR3 R:S ratio (kabat)",
	"FR3 Stop codons? (imgt)",
	"FR3 Stop codons? (kabat)",
	"FR3 base substitution freq% (imgt)",
	"FR3 base substitution freq% (kabat)",
	"FR3 base substitutions (imgt)",
	"FR3 base substitutions (kabat)",
	"FR3 deletion count (imgt)",
	"FR3 deletion count (kabat)",
	"FR3 homology% (imgt)",
	"FR3 homology% (kabat)",
	"FR3 indel frequency (imgt)",
	"FR3 indel frequency (kabat)",
	"FR3 insertion count (imgt)",
	"FR3 insertion count (kabat)",
	"FR3 length (imgt)",
	"FR3 length (kabat)",
	"FR3 nucleotide read (imgt)",
	"FR3 nucleotide read (kabat)",
	"FR3 replacement mutation freq% (imgt)",
	"FR3 replacement mutation freq% (kabat)",
	"FR3 replacement mutations (codons) (imgt)",
	"FR3 replacement mutations (codons) (kabat)",
	"FR3 silent mutation freq% (imgt)",
	"FR3 silent mutation freq% (kabat)",
	"FR3 silent mutations (codons) (imgt)",
	"FR3 silent mutations (codons) (kabat)",
	"Homology% (over V and J)",
	"Indel frequency (over V)",
	"Insertion count (over V)",
	"Length (over V and J)",
	"Mutations (over V and J)",
	"Mutations (over V)",
	"Nonsynonymous base substitutions (over V and J)",
	"Nonsynonymous base substitutions (over V)",
	"AA (over V and J)",
	"AA (over V)",
	"Nucleotide read (over V and J)",
	"Nucleotide read (over V)",
	"R:S ratio (over V and J)",
	"R:S ratio (over V)",
	"Replacement mutation freq% (over V and J)",
	"Replacement mutations (codons) (over V and J)",
	"Base substitutions (over V and J)",
	"CDR3 Stop codon? (over V and J)",
	"Deletion count (over V)",
	"CDR3 AA (kabat)",
	"Silent mutation freq% (over V and J)",
	"Silent mutations (codons) (over V and J)",
	"Base substitution freq% (over V and J)",
	"Stop codons? (over V)",
	"Synonymous base substitutions (over V and J)",
	"Synonymous base substitutions (over V)",
	"V gene hit/score list",
	"V gene tie",
	"D gene hit/score list",
	"D gene tie",
	"J gene hit/score list",
	"J gene tie"
	]

	m[keys[0]]=rid
	#keys=m.keys()
	#keys.sort()

	

	if(rid==1):
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
	parser.description+=' Generate read-level repertoire-characterization data.'
	parser.add_argument('-json_out',type=str,nargs=1,default="/dev/stdout",help="output file for the JSON segment count IMGT hierarchy")
	parser.add_argument('-cdr3_hist_out',type=str,nargs=1,default="/dev/stdout",help="output file for the CDR3 histogram of lengths (both kabat and imgt systems)")
	parser.add_argument('-skip_char',action='store_true', default=False,help="If this is set, then characterization besides asignment and CDR3 length is skipped")
	parser.add_argument('vdj_db',type=str,nargs=1,help="path to the VDJ root REQUIRED")
	parser.add_argument('qry_fasta',type=str,nargs=1,help="path to the input fasta file of query (Rep-Seq) data input to IgBLAST")
	#parser.add_argument('-region_out'
	parser.add_argument('db_organism',type=str,nargs=1,default="human",help="db organism",choices=["human","Mus_musculus"])
	return parser


		


if (__name__=="__main__"):
	parser=makeParserArgs()
	parser=add_rep_char_args_to_parser(parser)
	args = parser.parse_args()
	if(args):
		organism=extractAsItemOrFirstFromList(args.db_organism)
		args.db_species=args.db_organism
		mf=makeVDJMLDefaultMetaAndFactoryFromArgs(args)
		meta=mf[0]
		fact=mf[1]
		imgtdb_obj=imgt_db(extractAsItemOrFirstFromList(args.vdj_db))
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
		read_num=1
		segments_json_out=extractAsItemOrFirstFromList(args.json_out)
		#print "to write vdjml to ",args.vdjml_out
		rrw = vdjml.Vdjml_writer(extractAsItemOrFirstFromList(args.vdjml_out), meta)
		rep_char_out=extractAsItemOrFirstFromList(args.char_out)
		fHandle=open(rep_char_out,'w')
		logFile=rep_char_out+".log"
		logHandle=open(logFile,'w')
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

			#process for writing
			rrw(read_result_obj)

			#write rep-char
			appendAnnToFileWithMap(fHandle,read_analysis_results['ann_map'],read_num,None,"None",logHandle)

			#increment the read number
			read_num+=1

		print "Processed total of",str(read_num-1),"reads!"

		#close rep_char out
		fHandle.close()

		#close log
		logHandle.close()

		#write the CDR3 hist when non-dev-null
		if(type(cdr3_hist_out_file)==list):
			cdr3_hist_out_file=extractAsItemOrFirstFromList(cdr3_hist_out_file)
		if(not(cdr3_hist_out_file=="/dev/null")):
			my_cdr3_map.writeToFile(cdr3_hist_out_file)
			print "Wrote CDR3 lengths histogram to ",cdr3_hist_out_file


		#write the segment counts when non-dev-null
		if(type(segments_json_out)==list):
			segments_json_out=segments_json_out
		if(not(segments_json_out=="/dev/null")):
			print "Writing JSON segment counts output to ",segments_json_out,"..."
			segment_counter.JSONIFYToFile(extractAsItemOrFirstFromList(args.vdj_db),organism,segments_json_out)
			print "Writing JSON segment counts output complete!"

	else:
		#print "error in args!"
		parser.print_help()




