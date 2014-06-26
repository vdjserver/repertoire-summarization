#!/usr/bin/env python

from pprint import pprint
from utils import *
from segment_utils import getNAProductiveRearrangmentFlagFromVJHitLineData,getRegionAlignmentFromLargerVAlignment
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






def buildAlignmentWholeSeqsVDJML(btop,s):
	vdjml_btop=vdjml.Btop(btop)
	vdjml_aln=vdjml.sequence_match(vdjml_btop,gl_seq=s)
	aln=list()
	aln.append(vdjml_aln.seq_[0])
	aln.append("")
	aln.append(vdjml_aln.seq_[1])
	return aln
#        
#           0123456789
#           .....AA... read
#           .....C-... germline
#           012345 678
#        
#  	 b = vdjml.Btop('5ACA-3')
#        s = vdjml.sequence_match(
#                                b, 
#                                gl_seq='gggggCggg'
#                                )
#        self.assertEqual(s.start_[0], 0)
#        self.assertEqual(s.start_[1], 0)
#        self.assertEqual(s.end_[0], 10)
#        self.assertEqual(s.end_[1], 9)
#        self.assertEqual(s.seq_[0], 'gggggAAggg')
#        self.assertEqual(s.seq_[1], 'gggggC-ggg')





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
						sbjct=imgtdb_obj.getRefDirSetFNASeqGivenOrgAndAllele(ret_map['subject ids'],organism)
						#print "current got subject ",sbjct
						sbjct_for_btop=sbjct[s_start_line-1:s_end_line]
						#print "from map :"
						#printMap(ret_map)
						#print "from full query=",query_str
						#print "from full sbjct=",sbjct
						#print "got btop ready query=",query_for_btop," of len ",len(query_for_btop)
						#print "got btop ready sbjct=",sbjct_for_btop," of len ",len(sbjct_for_btop)
						#q_m_s=buildAlignmentWholeSeqs(btop,query_for_btop,sbjct_for_btop)
						q_m_s_vdjml=buildAlignmentWholeSeqsVDJML(btop,sbjct_for_btop)
						#print q_s[0]
						#print q_s[1]
						#print q_s[2]
						#sys.exit(0)
						#print "VDJML : \n",q_m_s_vdjml,"\n"
						#print "CURRENT : \n",q_m_s,"\n"
						#assert q_m_s_vdjml[0]==q_m_s[0]	#assert VDJML-generated is same as currently is
						#assert q_m_s_vdjml[2]==q_m_s[2]

						
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
	if(between(sub_aln_strt,ref_region_interval[0],sub_aln_end) or between(sub_aln_strt,ref_region_interval[1],sub_aln_end)):
		#good!  some part of the region falls within the subject alignment!
		rel_reg_strt=ref_region_interval[0]-sub_aln_strt
		rel_reg_strt=max(0,rel_reg_strt)
		rel_reg_ends=min(rel_reg_strt+reg_region_len,sub_aln_ends-sub_aln_strt)
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




#return a map assiigning BLAST scores to germline segment hits
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
	



