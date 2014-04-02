#!/usr/bin/env python

from pprint import pprint
from utils import printMap,repStr
from segment_utils import getNAProductiveRearrangmentFlagFromVJHitLineData,getRegionAlignmentFromLargerVAlignment
from vdjml_igblast_parse import extractJunctionRegionSeq,rev_comp_dna
from igblast_utils import buildAlignmentWholeSeqs
from alignment import alignment
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
			for gls_match in segment_match.germline_segments():
				#print "in inner most loop...."
				#print "gls_match=",gls_match
				hit_name=meta[gls_match.gl_segment_].name_
				#print "the hit name is ",hit_name
				if(hit_name==alleleName):
					#print "TARGET FOUND!"
					match_range=gls_match.range_
					#print "got a range : ",match_range
					pos1=match_range.pos_1()
					pose=pos1+match_range.length()
					#print "retrieved ",pos1
					#pprint(match_range)
					#match_pos1=match_range.first_last_1
					#print "p1=",match_pos1("pos1")
					#match_len=match_range.length()
					ret_map['s. start']=int(pos1)
					ret_map['s. end']=int(pose)-1
					#	1       59      90      148
					query_range=segment_match.range()
					query_start=query_range.pos_1()
					query_end=query_start+query_range.length()
					ret_map['q. start']=query_start
					ret_map['q. end']=query_end-1
					btop=str(segment_match.btop())
					ret_map['btop']=str(btop)
					ret_map['subject ids']=alleleName
					ret_map['is_inverted']=gls_match.inverted()
					ret_map['query id']=read_result_obj.id()
					if(ret_map['is_inverted']):
						ret_map['query id']="reversed|"+read_result_obj.id()
					if(query_record!=None):
						#attempt to use BTOP
						q_start_line=ret_map['q. start']
						q_end_line=ret_map['q. end']
						s_start_line=ret_map['s. start']
						s_end_line=ret_map['s. end']
						query_str=str(query_record.seq)
						if(ret_map['is_inverted']):
							query_str=rev_comp_dna(query_str)
						query_for_btop=query_str[q_start_line-1:q_end_line]
						sbjct=imgtdb_obj.getRefDirSetFNASeqGivenOrgAndAllele(ret_map['subject ids'],organism)
						sbjct_for_btop=sbjct[s_start_line-1:s_end_line]
						#print "from map :"
						#printMap(ret_map)
						#print "from full query=",query_str
						#print "from full sbjct=",sbjct
						#print "got btop ready query=",query_for_btop," of len ",len(query_for_btop)
						#print "got btop ready sbjct=",sbjct_for_btop," of len ",len(sbjct_for_btop)
						q_m_s=buildAlignmentWholeSeqs(btop,query_for_btop,sbjct_for_btop)
						#print q_s[0]
						#print q_s[1]
						#print q_s[2]
						#sys.exit(0)
						ret_map['query seq']=q_m_s[0]
						ret_map['subject seq']=q_m_s[2]
					return ret_map
				else:
					#print "skipping ",hit_name
					pass
	#print "For ",alleleName," the map is :"
	#printMap(ret_map)
	return ret_map



#return either true or false, no N/A yet
def getProductiveRearrangmentFlag(read_result_obj,meta,organism,imgtdb_obj):
	topVDJ=getTopVDJItems(read_result_obj,meta)
	productive_flag=False
	top_V=topVDJ['V']
	top_J=topVDJ['J']
	if(top_V!=None and top_J!=None):
		vHitInfo=getHitInfo(read_result_obj,meta,top_V)
		jHitInfo=getHitInfo(read_result_obj,meta,top_J)
		productive_flag=getNAProductiveRearrangmentFlagFromVJHitLineData(vHitInfo,jHitInfo,organism,imgtdb_obj)
		if(productive_flag is None):
			#there is no "N/A" yet!
			return False
		else:
			return productive_flag
	else:
		#note, there is no "N/A" or "inconclusive" ....
		return productive_flag




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



#from a segment combination, extract the region
def getRegionsObjsFromSegmentCombo(segment_combo,meta,ann_map):
	print "got a segment combo : ",segment_combo
	print "regions is ",segment_combo.regions
	print "now calling....",
	regions=segment_combo.regions()
	print "the type of regions is ",str(type(regions))
	for region in regions:
		#print "region???",region
		print "range : ",region.range_
		#print region.range_.pos_0()
		rgn_start=region.range_.pos_1() 
		rgn_end=rgn_start+region.range_.length()-1
		print "start/end : ",rgn_start," and ",rgn_end
		ann_map_key='vdj_server_ann_'
		





#using the read object and meta data
#find the names of the segments in the top
#VDJ combintion and add these it
def getTopVDJItems(read_result_obj,meta):
	#print "I want top hits from a read object!"
	top_segs=dict()
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
		seg_id=0
		for i in segment_combination.segments():
			segment_match = read_result_obj[i]
			#print "got a segment match, id=",i," from combination # ",str(int(s+1))
			#seg_type=segTypes[seg_id]
			#seg_type=segment_match.vdj
			#seg_type=segment_match.vdj_
			#print "the seg type is ",seg_type
			for gls_match in segment_match.germline_segments():
				#seg_type=gls_match.segment_type_
				#print "in inner most loop...."
				#print "gls_match=",gls_match
				#print meta[gls_match.gl_segment_].name_
				#type_from_meta=meta[gls_match.gl_segment_].gst_
				#print "META_TYPE",
				vdjml_look_seg_type=vdjml.segment_type(gls_match,meta)
				if(names[vdjml_look_seg_type]==None):
					#print "USING ",meta[gls_match.gl_segment_].name_
					names[vdjml_look_seg_type]=meta[gls_match.gl_segment_].name_
				else:
					#print "SKIPPING ",meta[gls_match.gl_segment_].name_," fail on emptiness"
					pass
				seg_id+=1
		#print "\n\n\n\n\n"
	#print "Returning from TOP : "
	#printMap(names)
	return names
	



