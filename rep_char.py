#!/usr/bin/env python

import vdjml
from vdjml_igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs
from utils import printMap,get_domain_modes,biopythonTranslate,makeAllMapValuesVal,getNumberBpInAlnStr,repStr
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
from vdjml_utils import getTopVDJItems,getRegionsObjsFromSegmentCombo,getHitInfo,getProductiveRearrangmentFlag,getVDJServerRegionAlignmentFromLargerVAlignmentPyObj,getAlignmentString
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
		if(f!=(-1) and t!=(-1) and cdr3_length_results[mode]!=(-1) ):
			qw=str(read_rec.seq)
			if(cdr3_length_results['qry_rev']):
				qw=rev_comp_dna(qw)
			query_cdr3=qw[f-1:t]
			read_ann_map[global_key_base+mode+'_cdr3_na']=query_cdr3
			read_ann_map[global_key_base+mode+'_cdr3_tr']=biopythonTranslate(read_ann_map[global_key_base+mode+'_cdr3_na'])
			#add in lengths!
			read_ann_map[global_key_base+mode+'_cdr3_na_len']=len(query_cdr3)
			read_ann_map[global_key_base+mode+'_cdr3_tr_len']=len(read_ann_map[global_key_base+mode+'_cdr3_tr'])
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
			read_ann_map[global_key_base+mode+'_cdr3_tr']=""
			#add in lengths!
			read_ann_map[global_key_base+mode+'_cdr3_na_len']=""
			read_ann_map[global_key_base+mode+'_cdr3_tr_len']=""
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
	print "THE CDR3 ALN IN V is (q top, s bottom)"
	print cdr3_q_aln
	print frame_mask
	print "stop here"
	print "last q",temp_v_q_pos
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
		print "Initial q_start is ",q_start
		#translate all the way to the end of the query here!!!!
		q_end=vAlnNoGapStart.q_end
		#if(vInfo['is_inverted']):
		#	q_start=getRegPosFromInvertedPos(q_start,len(seq_rec.seq))
		#	print "q_start after inversion : ",q_start
		#	q_end=getRegPosFromInvertedPos(q_end,len(seq_rec.seq))
		if(not(jInfo==None)):
			jAlnObj=alignment(jInfo['query seq'],jInfo['subject seq'],jInfo['q. start'],jInfo['q. end'],jInfo['s. start'],jInfo['s. end'])
			j_q_end=jAlnObj.q_end
			print "Initial q_end from j",j_q_end			
			#if(jInfo['is_inverted']):
			#	j_q_end=getRegPosFromInvertedPos(j_q_end,len(seq_rec.seq))
			print "Used q_end from j",j_q_end
			#the purpose of this "max" function is to let q_end stay where it is in the case the J aligns BEFORE V does!
			q_end=max(q_end,j_q_end)
		else:
			#no J alignment so try to get end of D alignment
			if(not(dInfo==None)):
				dAlnObj=alignment(dInfo['query seq'],dInfo['subject seq'],dInfo['q. start'],dInfo['q. end'],dInfo['s. start'],dInfo['s. end'])
				d_q_end=dAlnObj.q_end
				print "Initial q_end from d",d_q_end
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
		print "The query and translation to examine for ",seq_rec.id," (with q_start=",q_start," and q_end=",q_end,")  :"
		print query_to_exm
		print trx_to_exm
		if(trx_to_exm.find("*")!=(-1)):
			return True
		else:
			return False
	else:
		return None
		




def returnWholeSeqCharMap(vInfo,jInfo,imgtdb_obj,organism,annMap):
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
			return preCDR3Aln.characterize()
		else:
			#work on V and J
			jAlnObj=alignment(jInfo['query seq'],jInfo['subject seq'],jInfo['q. start'],jInfo['q. end'],jInfo['s. start'],jInfo['s. end'])
			J_cdr3_end=getADJCDR3EndFromJAllele(jInfo['subject ids'],imgtdb_obj,organism,"imgt")
			#print "jAlnObj"
			#print jAlnObj.getNiceString()
			if(J_cdr3_end==(-1)):
				#could not detect CDR3 end! so just return V map data
				return preCDR3Aln.characterize()
			elif(jAlnObj.q_end<=preCDR3Aln.q_start or jAlnObj.q_end<=preCDR3Aln.q_end):
				#j alignment within V so just return V ; avoid double countings
				return preCDR3Aln.characterize()
			else:
				postCDR3AlnStart=J_cdr3_end+1
				#okay, got J alignment post-CDR3 position
			if(postCDR3AlnStart>jAlnObj.s_end):
				#CDR3 start is after the s_aln end which means J aligns only in CDR3, ....so just return V map
				return preCDR3Aln.characterize()
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
				print "postCDR3AlnStart=",postCDR3AlnStart," and postCDR3AlnObj.s_start=",postCDR3AlnObj.s_start," so setting SFM as ",str(int(int(postCDR3AlnObj.s_start-postCDR3AlnStart)%3))
				postCDR3AlnObj.setSFM((postCDR3AlnObj.s_start-postCDR3AlnStart)%3)
			else:
				#postCDR3AlnStart>=postCDR3AlnObj.s_start ?????
				print "postCDR3AlnStart=",postCDR3AlnStart
				print "postCDR3AlnObj.s_start",postCDR3AlnObj.s_start
				print "vINFO",vInfo
				print "jINFO",jInfo
				print "J_cdr3_end=",J_cdr3_end
				print "VALNOBJ"
				print vAlnObj.getNiceString()
				print "postCDR3"
				print postCDR3AlnObj.getNiceString()
				raise Exception('postCDR3AlnStart>postCDR3AlnObj.s_start ?????  Error in alignment class???')
				#ALIGNMENT. Q_FROM=45 Q_TO=51 S_FROM=19 S_TO=25
				#QURY:GCAGGAG
				#MDLN:X||X|XX
				#SBCT:ACATGGA

				
			numN=postCDR3AlnObj.q_start-preCDR3Aln.q_end-1
			NSub=repStr("N",numN)
			NQry=repStr("N",numN)
			vCharMap=preCDR3Aln.characterize()
			print "the vCharMap annotation map and alignment (for whole seq named ",vInfo['query id'],") is "
			printMap(vCharMap)
			print preCDR3Aln.getNiceString()
			jCharMap=postCDR3AlnObj.characterize()
			#to get/compute the 'whole' map, add up the "non-special" values
			#otherwise, compute them specially
			toSkip=['bsb_freq','pct_id','ns_rto']
			totCharMap=dict()
			for k in vCharMap:
				if(not(k in toSkip)):
					totCharMap[k]=vCharMap[k]+jCharMap[k]
			if(vCharMap['stp_cdn'] or jCharMap['stp_cdn']):
				totCharMap['stp_cdn']=True
			
			if('vdj_server_ann_imgt_cdr3_tr' in annMap):
				starPos=annMap['vdj_server_ann_imgt_cdr3_tr'].find("*")
				if(starPos!=(-1)):
					totCharMap['cdr3_stp_cdn']=True
				else:
					totCharMap['cdr3_stp_cdn']=False
			else:
				totCharMap['cdr3_stp_cdn']=False
				
			#handle pct id and bsb_freq specially
			totQ=getNumberBpInAlnStr(postCDR3AlnObj.q_aln)+getNumberBpInAlnStr(preCDR3Aln.q_aln)
			totBS=totCharMap['base_sub']
			if(totQ!=0):
				totCharMap['pct_id']=float((float(totQ)-float(totBS))/float(totQ))
				totCharMap['bsb_freq']=float(totCharMap['base_sub'])/float(totQ)
			else:
				totCharMap['pct_id']=0.0
				totCharMap['bsb_freq']=0
			#handle ns_ratio
			if(totCharMap['synonymous_bsb']!=0):
				totCharMap['ns_rto']=float(totCharMap['nonsynonymous_bsb'])/float(totCharMap['synonymous_bsb'])
			else:
				totCharMap['ns_rto']=0
			return totCharMap
	emptyMap=getEmptyRegCharMap()
	return emptyMap
	


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
		if(topVDJ[seg] is not None):
			annMap['top_'+seg]=topVDJ[seg]
		else:
			annMap['top_'+seg]="None"			

	vInfo=None
	dInfo=None
	jInfo=None
	noneSeg_flag=False
	for seg in topVDJ:
		print "in second loop...."
		print "seg=",seg
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
		#whole seq characterization
		whole_char_map=returnWholeSeqCharMap(vInfo,jInfo,imgtdb_obj,organism,annMap)
		for w in whole_char_map:
			new_key=global_key_base+'whole_seq_'+w
			annMap[new_key]=whole_char_map[w]
		whole_seq_stp_cdn_Tot_flag=getValnWholeSeqStopFlag(vInfo,dInfo,jInfo,imgtdb_obj,organism,annMap,read_rec)
		annMap['vdj_server_whole_vj_stp_cdn']=whole_seq_stp_cdn_Tot_flag

		#productive rearrangement 
		#if IMGT CDR3 length is a multiple of 3 and it's not (-1), then consider it productive
		if(not(noneSeg_flag)):
			if(cdr3_map['imgt']!=(-1) and cdr3_map['imgt']>0):
				if((cdr3_map['imgt'])%3==0):
					annMap['vdj_server_ann_productive_rearrangement']=True
				else:
					annMap['vdj_server_ann_productive_rearrangement']=False
			else:
				#annMap['vdj_server_ann_productive_rearrangement']="N/A"
				pass

		#regions characterization FR1, FR2, FR3, CDR1, CDR2 (imgt and kabat)
		mode_list=get_domain_modes()

		if(characterization_thread_set==None):
			num_cpus=multiprocessing.cpu_count()
			num_threads=max(1,num_cpus)
			#num_threads=1
			#num_threads=10
			print "USING NUM_THREADS=",num_threads
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
						print "JUST SUBMITTED A CHAR_QUEUE JOB :",char_job_dict
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
		print "CHAR_QUEUE pre-join...."
		characterization_queue.join()
		get_res=0
		while(get_res<num_submitted_jobs):
			temp_results=characterization_queue_results.get()
			if(temp_results is not None):
				for temp_key in temp_results:
					annMap[temp_key]=temp_results[temp_key]			
			get_res+=1
	
	
	annMap['read_name']=read_rec.id
	#getAlignmentString(read_result_obj,meta,query_record,imgtdb_obj,organism)
	#analyze_combinations(read_result_obj,meta,organism,imgtdb_obj,read_rec,annMap)
	printMap(annMap,True)
	return annMap





def analyze_combinations(read_result_obj,meta,organism,imgtdb_obj,read_rec,read_ann_map):
	segment_combinations=read_result_obj.segment_combinations()
	for s in range(len(segment_combinations)):
		print "LOOKING AT COMBINATION # ",str(int(s+1))," of ",str(len(segment_combinations))," FOR READ ID=",read_result_obj.id()
		segment_combination=segment_combinations[s]
		analyzeRegionObjsFromSegmentCombo(segment_combination,meta,organism,imgtdb_obj,read_rec,read_ann_map)




def analyzeRegionObjsFromSegmentCombo(segment_combo,meta,organism,imgtdb_obj,read_rec,read_ann_map):
	global global_key_base
	print "got a segment combo : ",segment_combo
	print "regions is ",segment_combo.regions
	print "now calling....",
	regions=segment_combo.regions()
	print "the type of regions is ",str(type(regions))
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
				print "total match!"
			else:
				print "at least one mismatch! igblast start/end : ",igblast_start," and ",igblast_end," vdj start/end     : ",vdj_start," and ",vdj_end," : query=",read_ann_map['read_name']," seg=",read_ann_map['top_V']," rgn=",rgn_nm," rrv : ",getVRegionStartAndStopGivenRefData(read_ann_map['top_V'],organism,imgtdb_obj,rgn_nm,rgn_ns)
		else:
			print "no lookup found for ",lookup_qry_start," and ",lookup_qry_end
			pass



		









#using the PYVDJML, extract information from the IGBLAST-notated
#regions (FR1, CDR1, etc) and compare them with VDJ-server lookedup values
def vSegmentRegionVDJAnalyse(read_result_obj,meta,organism,imgtdb_obj,read_rec):
	topVDJ=getTopVDJItems(read_result_obj,meta)
	topV=topVDJ['V']
	segment_combinations=read_result_obj.segment_combinations()
	for s in range(len(segment_combinations)):
		segment_combination=segment_combinations[s]
		print "Looking at combination #",str(s+1)
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





def appendAnnToFileWithMap(fHandl,m,rid,desiredKeys=None,defaultValue="None"):
	keys=[
	"read_id#",
	"read_name",
	"top_D",
	"top_J",
	"top_V",
	"vdj_server_ann_imgt_CDR1_base_sub",
	"vdj_server_ann_imgt_CDR1_bsb_freq",
	"vdj_server_ann_imgt_CDR1_deletions",
	"vdj_server_ann_imgt_CDR1_frm_msk",
	"vdj_server_ann_imgt_CDR1_indel_freq",
	"vdj_server_ann_imgt_CDR1_insertions",
	"vdj_server_ann_imgt_CDR1_mutations",
	"vdj_server_ann_imgt_CDR1_nonsynonymous_bsb",
	"vdj_server_ann_imgt_CDR1_ns_rto",
	"vdj_server_ann_imgt_CDR1_pct_id",
	"vdj_server_ann_imgt_CDR1_qry_aln",
	"vdj_server_ann_imgt_CDR1_qry_end",
	"vdj_server_ann_imgt_CDR1_qry_srt",
	"vdj_server_ann_imgt_CDR1_ref_aln",
	"vdj_server_ann_imgt_CDR1_stp_cdn",
	"vdj_server_ann_imgt_CDR1_synonymous_bsb",
	"vdj_server_ann_imgt_CDR2_base_sub",
	"vdj_server_ann_imgt_CDR2_bsb_freq",
	"vdj_server_ann_imgt_CDR2_deletions",
	"vdj_server_ann_imgt_CDR2_frm_msk",
	"vdj_server_ann_imgt_CDR2_indel_freq",
	"vdj_server_ann_imgt_CDR2_insertions",
	"vdj_server_ann_imgt_CDR2_mutations",
	"vdj_server_ann_imgt_CDR2_nonsynonymous_bsb",
	"vdj_server_ann_imgt_CDR2_ns_rto",
	"vdj_server_ann_imgt_CDR2_pct_id",
	"vdj_server_ann_imgt_CDR2_qry_aln",
	"vdj_server_ann_imgt_CDR2_qry_end",
	"vdj_server_ann_imgt_CDR2_qry_srt",
	"vdj_server_ann_imgt_CDR2_ref_aln",
	"vdj_server_ann_imgt_CDR2_stp_cdn",
	"vdj_server_ann_imgt_CDR2_synonymous_bsb",
	"vdj_server_ann_imgt_FR1_base_sub",
	"vdj_server_ann_imgt_FR1_bsb_freq",
	"vdj_server_ann_imgt_FR1_deletions",
	"vdj_server_ann_imgt_FR1_frm_msk",
	"vdj_server_ann_imgt_FR1_indel_freq",
	"vdj_server_ann_imgt_FR1_insertions",
	"vdj_server_ann_imgt_FR1_mutations",
	"vdj_server_ann_imgt_FR1_nonsynonymous_bsb",
	"vdj_server_ann_imgt_FR1_ns_rto",
	"vdj_server_ann_imgt_FR1_pct_id",
	"vdj_server_ann_imgt_FR1_qry_aln",
	"vdj_server_ann_imgt_FR1_qry_end",
	"vdj_server_ann_imgt_FR1_qry_srt",
	"vdj_server_ann_imgt_FR1_ref_aln",
	"vdj_server_ann_imgt_FR1_stp_cdn",
	"vdj_server_ann_imgt_FR1_synonymous_bsb",
	"vdj_server_ann_imgt_FR2_base_sub",
	"vdj_server_ann_imgt_FR2_bsb_freq",
	"vdj_server_ann_imgt_FR2_deletions",
	"vdj_server_ann_imgt_FR2_frm_msk",
	"vdj_server_ann_imgt_FR2_indel_freq",
	"vdj_server_ann_imgt_FR2_insertions",
	"vdj_server_ann_imgt_FR2_mutations",
	"vdj_server_ann_imgt_FR2_nonsynonymous_bsb",
	"vdj_server_ann_imgt_FR2_ns_rto",
	"vdj_server_ann_imgt_FR2_pct_id",
	"vdj_server_ann_imgt_FR2_qry_aln",
	"vdj_server_ann_imgt_FR2_qry_end",
	"vdj_server_ann_imgt_FR2_qry_srt",
	"vdj_server_ann_imgt_FR2_ref_aln",
	"vdj_server_ann_imgt_FR2_stp_cdn",
	"vdj_server_ann_imgt_FR2_synonymous_bsb",
	"vdj_server_ann_imgt_FR3_base_sub",
	"vdj_server_ann_imgt_FR3_bsb_freq",
	"vdj_server_ann_imgt_FR3_deletions",
	"vdj_server_ann_imgt_FR3_frm_msk",
	"vdj_server_ann_imgt_FR3_indel_freq",
	"vdj_server_ann_imgt_FR3_insertions",
	"vdj_server_ann_imgt_FR3_mutations",
	"vdj_server_ann_imgt_FR3_nonsynonymous_bsb",
	"vdj_server_ann_imgt_FR3_ns_rto",
	"vdj_server_ann_imgt_FR3_pct_id",
	"vdj_server_ann_imgt_FR3_qry_aln",
	"vdj_server_ann_imgt_FR3_qry_end",
	"vdj_server_ann_imgt_FR3_qry_srt",
	"vdj_server_ann_imgt_FR3_ref_aln",
	"vdj_server_ann_imgt_FR3_stp_cdn",
	"vdj_server_ann_imgt_FR3_synonymous_bsb",
	"vdj_server_ann_imgt_cdr3_na",
	"vdj_server_ann_imgt_cdr3_na_len",
	"vdj_server_ann_imgt_cdr3_tr",
	"vdj_server_ann_imgt_cdr3_tr_len",
	"vdj_server_ann_kabat_CDR1_base_sub",
	"vdj_server_ann_kabat_CDR1_bsb_freq",
	"vdj_server_ann_kabat_CDR1_deletions",
	"vdj_server_ann_kabat_CDR1_frm_msk",
	"vdj_server_ann_kabat_CDR1_indel_freq",
	"vdj_server_ann_kabat_CDR1_insertions",
	"vdj_server_ann_kabat_CDR1_mutations",
	"vdj_server_ann_kabat_CDR1_nonsynonymous_bsb",
	"vdj_server_ann_kabat_CDR1_ns_rto",
	"vdj_server_ann_kabat_CDR1_pct_id",
	"vdj_server_ann_kabat_CDR1_qry_aln",
	"vdj_server_ann_kabat_CDR1_qry_end",
	"vdj_server_ann_kabat_CDR1_qry_srt",
	"vdj_server_ann_kabat_CDR1_ref_aln",
	"vdj_server_ann_kabat_CDR1_stp_cdn",
	"vdj_server_ann_kabat_CDR1_synonymous_bsb",
	"vdj_server_ann_kabat_CDR2_base_sub",
	"vdj_server_ann_kabat_CDR2_bsb_freq",
	"vdj_server_ann_kabat_CDR2_deletions",
	"vdj_server_ann_kabat_CDR2_frm_msk",
	"vdj_server_ann_kabat_CDR2_indel_freq",
	"vdj_server_ann_kabat_CDR2_insertions",
	"vdj_server_ann_kabat_CDR2_mutations",
	"vdj_server_ann_kabat_CDR2_nonsynonymous_bsb",
	"vdj_server_ann_kabat_CDR2_ns_rto",
	"vdj_server_ann_kabat_CDR2_pct_id",
	"vdj_server_ann_kabat_CDR2_qry_aln",
	"vdj_server_ann_kabat_CDR2_qry_end",
	"vdj_server_ann_kabat_CDR2_qry_srt",
	"vdj_server_ann_kabat_CDR2_ref_aln",
	"vdj_server_ann_kabat_CDR2_stp_cdn",
	"vdj_server_ann_kabat_CDR2_synonymous_bsb",
	"vdj_server_ann_kabat_FR1_base_sub",
	"vdj_server_ann_kabat_FR1_bsb_freq",
	"vdj_server_ann_kabat_FR1_deletions",
	"vdj_server_ann_kabat_FR1_frm_msk",
	"vdj_server_ann_kabat_FR1_indel_freq",
	"vdj_server_ann_kabat_FR1_insertions",
	"vdj_server_ann_kabat_FR1_mutations",
	"vdj_server_ann_kabat_FR1_nonsynonymous_bsb",
	"vdj_server_ann_kabat_FR1_ns_rto",
	"vdj_server_ann_kabat_FR1_pct_id",
	"vdj_server_ann_kabat_FR1_qry_aln",
	"vdj_server_ann_kabat_FR1_qry_end",
	"vdj_server_ann_kabat_FR1_qry_srt",
	"vdj_server_ann_kabat_FR1_ref_aln",
	"vdj_server_ann_kabat_FR1_stp_cdn",
	"vdj_server_ann_kabat_FR1_synonymous_bsb",
	"vdj_server_ann_kabat_FR2_base_sub",
	"vdj_server_ann_kabat_FR2_bsb_freq",
	"vdj_server_ann_kabat_FR2_deletions",
	"vdj_server_ann_kabat_FR2_frm_msk",
	"vdj_server_ann_kabat_FR2_indel_freq",
	"vdj_server_ann_kabat_FR2_insertions",
	"vdj_server_ann_kabat_FR2_mutations",
	"vdj_server_ann_kabat_FR2_nonsynonymous_bsb",
	"vdj_server_ann_kabat_FR2_ns_rto",
	"vdj_server_ann_kabat_FR2_pct_id",
	"vdj_server_ann_kabat_FR2_qry_aln",
	"vdj_server_ann_kabat_FR2_qry_end",
	"vdj_server_ann_kabat_FR2_qry_srt",
	"vdj_server_ann_kabat_FR2_ref_aln",
	"vdj_server_ann_kabat_FR2_stp_cdn",
	"vdj_server_ann_kabat_FR2_synonymous_bsb",
	"vdj_server_ann_kabat_FR3_base_sub",
	"vdj_server_ann_kabat_FR3_bsb_freq",
	"vdj_server_ann_kabat_FR3_deletions",
	"vdj_server_ann_kabat_FR3_frm_msk",
	"vdj_server_ann_kabat_FR3_indel_freq",
	"vdj_server_ann_kabat_FR3_insertions",
	"vdj_server_ann_kabat_FR3_mutations",
	"vdj_server_ann_kabat_FR3_nonsynonymous_bsb",
	"vdj_server_ann_kabat_FR3_ns_rto",
	"vdj_server_ann_kabat_FR3_pct_id",
	"vdj_server_ann_kabat_FR3_qry_aln",
	"vdj_server_ann_kabat_FR3_qry_end",
	"vdj_server_ann_kabat_FR3_qry_srt",
	"vdj_server_ann_kabat_FR3_ref_aln",
	"vdj_server_ann_kabat_FR3_stp_cdn",
	"vdj_server_ann_kabat_FR3_synonymous_bsb",
	"vdj_server_ann_kabat_cdr3_na",
	"vdj_server_ann_kabat_cdr3_na_len",
	"vdj_server_ann_kabat_cdr3_tr",
	"vdj_server_ann_kabat_cdr3_tr_len",
	"vdj_server_ann_productive_rearrangement",
	"vdj_server_ann_whole_seq_base_sub",
	"vdj_server_ann_whole_seq_bsb_freq",
	"vdj_server_ann_whole_seq_cdr3_stp_cdn",
	"vdj_server_ann_whole_seq_deletions",
	"vdj_server_ann_whole_seq_indel_freq",
	"vdj_server_ann_whole_seq_insertions",
	"vdj_server_ann_whole_seq_mutations",
	"vdj_server_ann_whole_seq_nonsynonymous_bsb",
	"vdj_server_ann_whole_seq_ns_rto",
	"vdj_server_ann_whole_seq_pct_id",
	"vdj_server_ann_whole_seq_stp_cdn",
	"vdj_server_ann_whole_seq_synonymous_bsb",
	"vdj_server_whole_vj_stp_cdn"
	]

	m[keys[0]]=rid
	#keys=m.keys()
	keys.sort()

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

		


#add on the kinds of arguments this program accepts
def add_rep_char_args_to_parser(parser):
	parser.description+=' Generate read-level repertoire-characterization data.'
	parser.add_argument('-json_out',type=str,nargs=1,default="/dev/stdout",help="output file for the JSON segment count IMGT hierarchy")
	parser.add_argument('-cdr3_hist_out',type=str,nargs=1,default="/dev/stdout",help="output file for the CDR3 histogram of lengths (both kabat and imgt systems)")
	parser.add_argument('-skip_char',action='store_true', default=False,help="If this is set, then characterization besides asignment and CDR3 length is skipped")
	parser.add_argument('vdj_db',type=str,nargs=1,help="path to the VDJ root REQUIRED")
	parser.add_argument('qry_fasta',type=str,nargs=1,help="path to the input fasta file of query (Rep-Seq) data input to IgBLAST")
	#parser.add_argument('-region_out'
	#parser.add_argument('-db_species',type=str,nargs=1,default="human",help="species of the db")
	return parser


		


if (__name__=="__main__"):
	parser=makeParserArgs()
	parser=add_rep_char_args_to_parser(parser)
	args = parser.parse_args()
	if(args): 
		mf=makeVDJMLDefaultMetaAndFactoryFromArgs(args)
		meta=mf[0]
		fact=mf[1]
		imgtdb_obj=imgt_db(args.vdj_db[0])
		#print "a fact is ",fact
		#print "the file is ",args.igblast_in[0]
		query_fasta=args.qry_fasta[0]
		fasta_reader=SeqIO.parse(open(query_fasta, "r"), "fasta")
		modes=['kabat','imgt']
		organism=args.db_species
		my_cdr3_map=histoMapClass(modes)
		cdr3_hist_out_file=args.cdr3_hist_out
		segment_counter=IncrementMapWrapper()
		read_num=1
		segments_json_out=args.json_out
		print "to write vdjml to ",args.vdjml_out
		rrw = vdjml.Vdjml_writer(str(args.vdjml_out[0]), meta)
		rep_char_out=args.char_out[0]
		fHandle=open(rep_char_out,'w')
		for read_result_obj in scanOutputToVDJML(args.igblast_in[0],fact,query_fasta):
			#prepare for the iteration and give a possible status message...
			if(read_num>1 and read_num%1000==0):
				print "Processed read",read_num,"..."
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

			appendAnnToFileWithMap(fHandle,read_analysis_results['ann_map'],read_num)

			#increment the read number
			read_num+=1

		#close rep_char out
		fHandle.close()

		#write the CDR3 hist when non-dev-null
		if(type(cdr3_hist_out_file)==list):
			cdr3_hist_out_file=cdr3_hist_out_file[0]
		if(not(cdr3_hist_out_file=="/dev/null")):
			my_cdr3_map.writeToFile(cdr3_hist_out_file)

		#write the segment counts when non-dev-null
		if(type(segments_json_out)==list):
			segments_json_out=segments_json_out[0]
		if(not(segments_json_out=="/dev/null")):
			segment_counter.JSONIFYToFile(args.vdj_db[0],organism,segments_json_out)

	else:
		#print "error in args!"
		parser.print_help()


