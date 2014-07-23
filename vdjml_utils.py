#!/usr/bin/env python

from alignment import alignment
from pprint import pprint
from utils import *
from segment_utils import getNAProductiveRearrangmentFlagFromVJHitLineData,getRegionAlignmentFromLargerVAlignment,getVRegionsList,getVRegionStartAndStopGivenRefData
from igblast_parse import extractJunctionRegionSeq,rev_comp_dna
import vdjml

def getAlignmentString(read_result_obj,meta,query_record=None,imgtdb_obj=None,organism=None):
	topVDJ=getTopVDJItems(read_result_obj,meta)
	topVAllele=topVDJ['V']
	if(topVAllele is not None):
		vMap=getHitInfo(read_result_obj,meta,topVDJ['V'],query_record,imgtdb_obj,organism)
		junction_names=['vd_junction','dj_junction','vj_junction']
		junction_regions=dict()
		junction_seqs=dict()
		junction_alns=dict()
		all_dict=dict()
		all_dict['V']=topVDJ['V']
		for j in range(len(junction_names)):
			junction_regions[junction_names[j]]=getJunctionRegionByname(read_result_obj,meta,junction_names[j])
			if(not(junction_regions[junction_names[j]]==None)):
				junction_seqs[junction_names[j]]=extractJunctionRegionSeq(junction_regions[junction_names[j]],query_record)
				junction_alns[junction_names[j]]=makeJuncAlignmentObjFromJuncRegion(junction_regions[junction_names[j]],junction_seqs[junction_names[j]])
			else:
				junction_seqs[junction_names[j]]=None
				junction_alns[junction_names[j]]=None
			all_dict[junction_names[j]]=junction_seqs[junction_names[j]]
		all_dict['D']=topVDJ['D']
		all_dict['J']=topVDJ['J']
		print "ALL_DICT",all_dict
	else:
		return ""






#given a read result, just print to the screen AA mut entires in the object
def showAAMuts(read_result_obj):
	segment_matches=read_result_obj.segment_matches()
	for segment_match in segment_matches:
		aa_subs=segment_match.aa_substitutions()
		num_subs=0
		for aa_sub in aa_subs:
			print "GOT A SUB!!!!!!"
			num_subs+=1
		print "FOUND ",num_subs," subs!"


#given a list of X123Y (where X and Y are amino acid single-letter codes)
#and given a segment match name, scan a VDJML object and where the segment
#match name is found, insert the AA
#NOTE that the positions are relative to the read
def addAAMuts(read_result_obj,aa_mut_list,target_seg_match_name,meta):
	segment_matches=read_result_obj.segment_matches()
	for segment_match in segment_matches:
		for gls_match in segment_match.gl_segments():
			hit_name=meta[gls_match.gl_segment()].name()
			print "Considering to add ",aa_mut_list," to ",hit_name
			if(hit_name==target_seg_match_name):
				print "found it..."
				#__init__(_object*, unsigned int read_pos0, std::string read_aa, std::string gl_aa)
				for AA_mut in aa_mut_list:
					aa_parts=getMutList(AA_mut)
					aa_from=aa_parts[0]
					aa_pos=aa_parts[1]
					aa_to=aa_parts[2]
					print "got ",aa_parts
					if(not(aa_to=="X")):
						myAASub=vdjml.Aa_substitution(int(aa_pos)-1,aa_to,aa_from)
						print "adding here...."
						segment_match.insert(myAASub)
						showAAMuts(read_result_obj)
			else:
				print "nope, this ain't it!"

#given a junction pyObject and the sequence return an alignment object
def makeJuncAlignmentObjFromJuncRegion(jregObj,jseq):
	rgn_start=jregObj.range_.pos_1() 
	rgn_end=rgn_start+jregObj.range_.length()-1
	q_aln=jseq
	s_aln=repStr("-",rgn_end-rgn_start+1)
	q_start=rgn_start
	q_end=rgn_end
	s_start=q_start
	s_end=q_end
	return alignment(q_aln,s_aln,q_start,q_end,s_start,s_end)


#given a read object and junction name extract the junction region from the identified combination ID
def getJunctionRegionByname(read_result_obj,meta,junc_name,combID=0):
	#print "Trying to get junction region by name ",junc_name
	segment_combinations=read_result_obj.segment_combinations()
	if(len(segment_combinations)==0):
		return None
		#print "no combo so returning none...."
	else:
		if(combID>=len(segment_combinations)):
			return None
		segment_combo=segment_combinations[combID]
		regions=segment_combo.regions()
		for region in regions:
			rgn_nm=meta[region.region_]
			if(rgn_nm==junc_name):
				return region
		return None






def buildAlignmentWholeSeqsVDJML(btop,s,isRead=False):
	vdjml_btop=vdjml.Btop(btop)
	if(not(isRead)):
		vdjml_aln=vdjml.sequence_match(vdjml_btop,gl_seq=s)
	else:
		vdjml_aln=vdjml.sequence_match(vdjml_btop,read_seq=s)		
	aln=list()
	aln.append(vdjml_aln.seq_[0])	#QUERY
	aln.append("")
	aln.append(vdjml_aln.seq_[1])	#SUBJECT
	return aln





#given a read object, meta object, and an allele name
#known to ALREADY HIT to the read, find the following info
#and put it in a dict 'q. start', 'q. end', 's. start', 's. end', 'query id', 'subject ids', & btop
def getHitInfo(read_result_obj,meta,alleleName,query_record=None,imgtdb_obj=None,organism=None):
	#print "Need to extract info for ",alleleName
	segment_matches=read_result_obj.segment_matches() 
	ret_map=dict()
	for segment_match in segment_matches:
			#segment_match = read_result_obj[i]
			#print "got a segment match, id=",i," from combination # ",str(int(s+1))
			for gls_match in segment_match.gl_segments():
				#print "in inner most loop...."
				#print "gls_match=",gls_match
				hit_name=meta[gls_match.gl_segment()].name()
				# print meta[gls_match.gl_segment()].name()
				#print "the hit name is ",hit_name
				if(hit_name==alleleName):
					#print "TARGET FOUND!"
					match_range=segment_match.gl_range()
					#print "got a range : ",match_range," for GL in ",alleleName," for query=",read_result_obj.id()
					pos1=match_range.pos1()
					pose=pos1+match_range.length()
					ret_map['s. start']=int(pos1)
					ret_map['s. end']=int(pose)-1
					query_range=segment_match.read_range()
					query_start=query_range.pos1()
					query_end=query_start+query_range.length()
					ret_map['q. start']=query_start
					ret_map['q. end']=query_end-1
					btop=str(segment_match.btop())
					ret_map['btop']=str(btop)
					ret_map['subject ids']=alleleName
					segment_match_metrics=segment_match.match_metrics()
					ret_map['is_inverted']=segment_match_metrics.inverted()
					ret_map['query id']=read_result_obj.id()
					if(ret_map['is_inverted']):
						ret_map['query id']="reversed|"+read_result_obj.id()
					if(query_record!=None):
						#attempt to use BTOP
						q_start_line=ret_map['q. start']
						q_end_line=ret_map['q. end']
						s_start_line=ret_map['s. start']
						s_end_line=ret_map['s. end']
						query_str=str(query_record.seq).upper()
						if(ret_map['is_inverted']):
							query_str=rev_comp_dna(query_str)
						query_for_btop=query_str[q_start_line-1:q_end_line]
						q_m_s_vdjml=buildAlignmentWholeSeqsVDJML(btop,query_for_btop,True)
						ret_map['query seq']=q_m_s_vdjml[0]   #query string/alignment
						ret_map['subject seq']=q_m_s_vdjml[2] #subject string/alignment
					return ret_map
				else:
					#print "skipping ",hit_name
					pass
	#print "For ",alleleName," the map is :"
	#printMap(ret_map)
	return ret_map



#is 'b' falling between 'a' and 'c'???
def between(a,b,c,exclu=True):
	if(exclu):
		if(a<b and b<c):
			return True
		else:
			return False
	else:
		if(a<=b and b<=c):
			return True
		else:
			return False
	


def getRegionAlignmentFromVDJML(read_result_obj,org,mode,region_name,imgtdb_obj,v_hit_info):
	valid_regions=getVRegionsList()
	if(not(region_name in valid_regions)):
		print region_name," is an invalid region!!!!"
		return None
	sub_name=v_hit_info['subject ids']
	ref_region_interval=getVRegionStartAndStopGivenRefData(sub_name,org,imgtdb_obj,region_name,mode)
	reg_region_len=abs(ref_region_interval[0]-ref_region_interval[1])
	ref_region_transcript_start=getVRegionStartAndStopGivenRefData(sub_name,org,imgtdb_obj,region_name,"imgt")[0]
	if(ref_region_interval[0]==(-1) or ref_region_interval[1]==(-1)):
		#print "Can't get region info, start or stop of ref region is neg 1...."
		return None	
	sub_aln_strt=v_hit_info['s. start']
	sub_aln_ends=v_hit_info['s. end']
	sub_aln_len=abs(sub_aln_ends-sub_aln_strt)
	#first condition  : region start is inside the sub-alignment
	#second condition : region end is within the sub-alignment
	#third  condition : this last condition is for the sub-alignment being wholly-contained within the region.  In this case, the sub-alignment is SMALL!
	if(
		between(sub_aln_strt,ref_region_interval[0],sub_aln_ends) 
		or 
		between(sub_aln_strt,ref_region_interval[1],sub_aln_ends) 
		or
		(between(ref_region_interval[0],sub_aln_strt,ref_region_interval[1]) and between(ref_region_interval[0],sub_aln_ends,ref_region_interval[1]))
		):
		#good!  some part of the region falls within the subject alignment!
		rel_reg_strt=(ref_region_interval[0]-sub_aln_strt)-1 #do the subtraction of 1, because vdjml uses 0-based offsets
		rel_reg_strt=max(0,rel_reg_strt) 
		#take the max here because the rel_reg_strt might be negative.  this is in case the region is partially covered in the alignment
		#at this point rel_reg_strt is in the space of the sub-alignment.  it's index is 0-based with respect to the sub-alignment in the reference.
		rel_reg_ends=rel_reg_strt+reg_region_len+1
		#now, get the reg end to also be in the sub_alignment space
		rel_reg_ends=min(rel_reg_ends,sub_aln_len+1)
		#use min here to make sure it doesn't go off the end of the sub-alignment!
		VDJML_aln=buildAlignmentWholeSeqsVDJML(v_hit_info['btop'],v_hit_info['subject seq'])
		vdjmlBTOP=vdjml.Btop(v_hit_info['btop'])
		print "REGION : ",region_name
		print "BTOP sub start and end : ",v_hit_info['s. start']," and ",v_hit_info['s. end']
		print "Rel reg start=",rel_reg_strt
		print "Rel reg ends =",rel_reg_ends
		s = vdjml.sequence_match(
			vdjmlBTOP,
			gl_start=rel_reg_strt,
			gl_end=rel_reg_ends,
			gl_seq=v_hit_info['subject seq']
			)
		subaln_reg_qry_strt=s.start_[0]	
		subaln_reg_ref_strt=s.start_[1]
		subaln_reg_qry_end=s.end_[0]
		subaln_reg_ref_end=s.end_[1]
		subaln_reg_query=s.seq_[0]
		subaln_reg_ref=s.seq_[1]
		reg_aln_obj=alignment(
			subaln_reg_query,
			subaln_reg_ref,
			subaln_reg_qry_strt,
			subaln_reg_qry_end,
			subaln_reg_ref_strt,
			subaln_reg_ref_end
			)
		#aquire frame information
		#region_frame_start=(ref_region_interval[0]-ref_region_transcript_start)%3
		region_frame_start=0		
		reg_aln_obj.setSFM(region_frame_start)
		print "A VDJML REGION : ",
		print reg_aln_obj.getNiceString()
		return reg_aln_obj
		
	else:
		#no part of the region falls within the subject alignment
		return None



#given a V segment alignment extract from it the sub-portion for a given
#region.  Return a 4-length array with subject and query alignment data  [region_alignment,frame_mask,r_q_start,r_q_end]  (region_alignment has query first, then sub)
#return None if alignment too short or doesn't cover region
def getVDJServerRegionAlignmentFromLargerVAlignmentPyObj(read_result_obj,meta,org,mode,region_name,imgtdb_obj,wholeOnly,read_rec):
	topVDJ=getTopVDJItems(read_result_obj,meta)
	if('V' in topVDJ):
		topV=topVDJ['V']
		if(topV!=None):
			v_info=getHitInfo(read_result_obj,meta,topV,read_rec,imgtdb_obj,org)
			raInfo=getRegionAlignmentFromLargerVAlignment(v_info,org,mode,region_name,imgtdb_obj,wholeOnly)
			#oldAlnString=raInfo.getNiceString()
			#print ">"+v_info['subject ids']+"\n"+imgtdb_obj.getRefDirSetFNASeqGivenOrgAndAllele(v_info['subject ids'],org)+"\n\n"
			#print "\n",oldAlnString
			#raInfoVDJML=getRegionAlignmentFromVDJML(read_result_obj,org,mode,region_name,imgtdb_obj,v_info)
			return raInfo
		else:
			return None
	else:
		return None




	
#map VDJ to True or False indicating a tie or not in the scores
def getVDJTieMap(t_map,topVDJ):
	#print "t_map is ",t_map
	tie_map=dict()
	for seg in t_map:
		#print "current seg : ",seg
		seg_map=t_map[seg]
		#print "seg_map ",seg_map
		score_list=list()
		for name in seg_map:
			score=seg_map[name]
			score_list.append(score)
		tie_map[seg]=False
		if(topVDJ!=None):
			if(seg in topVDJ):
				top_score=seg_map[topVDJ[seg]]
				score_count=score_list.count(top_score)
				if(score_count!=1 and score_count>1):
					tie_map[seg]=True
	return tie_map




#return a map assigning BLAST scores to germline segment hits
def getTypeNameScoreMap(read_result_obj,meta):
	t_map=dict()
	#print "in getHitsAndScores"
	segment_matches=read_result_obj.segment_matches()
	#print "got segment_matches ",segment_matches
	for segment_match in segment_matches:
		#print "Looking at a segment_match ",segment_match
		for gls_match in segment_match.gl_segments():
			vdjml_look_seg_type=vdjml.segment_type(gls_match,meta)
			#print "seg_type=",vdjml_look_seg_type
			name=meta[gls_match.gl_segment()].name()
			#print "The name is ",name
			segment_match_metrics=segment_match.match_metrics()
			score=segment_match_metrics.score()
			#print "score=",score
			if(vdjml_look_seg_type in t_map):
				pass
			else:
				t_map[vdjml_look_seg_type]=dict()
			t_map[vdjml_look_seg_type][name]=score
	#print "GRAND RESULT\n"
	#print t_map
	return t_map
		




#using the read object and meta data
#find the names of the segments in the top
#VDJ combintion and add these it
def getTopVDJItems(read_result_obj,meta):
	#print "I want top hits from a read object!"
	names=dict()
	segTypes=["V","D","J"]
	for st in segTypes:
		names[st]=None
	segment_combinations=read_result_obj.segment_combinations()
	for s in range(len(segment_combinations)):
		#print "LOOKING AT COMBINATION # ",str(int(s+1))," of ",str(len(segment_combinations))," FOR READ ID=",read_result_obj.id()
		segment_combination=segment_combinations[s]
		#print "about to enter into region extraction...."
		#getRegionsObjsFromSegmentCombo(segment_combination,meta)
		#print "The number of segments in this combination is ",len(segment_combination.segments())
		#print "The ids of the segments are ",segment_combination.segments()
		for i in segment_combination.segments():
			segment_match = read_result_obj[i]
			#print "got a segment match, id=",i," from combination # ",str(int(s+1))
			#print "the seg type is ",seg_type
			for gls_match in segment_match.gl_segments():
				#seg_type=gls_match.segment_type_
				#print "in inner most loop...."
				#print "gls_match=",gls_match
				#print meta[gls_match.gl_segment_].name_
				#type_from_meta=meta[gls_match.gl_segment_].gst_
				#print "META_TYPE",
				vdjml_look_seg_type=vdjml.segment_type(gls_match,meta)
				if(names[vdjml_look_seg_type]==None):
					# print meta[gls_match.gl_segment()].name()
					names[vdjml_look_seg_type]=meta[gls_match.gl_segment()].name()
				else:
					#print "SKIPPING ",meta[gls_match.gl_segment_].name_," fail on emptiness"
					pass
		#print "\n\n\n\n\n"
	#print "Returning from TOP : "
	#printMap(names)
	return names
	



