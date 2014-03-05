#!/usr/bin/env python

from pprint import pprint
from utils import printMap
from segment_utils import getNAProductiveRearrangmentFlagFromVJHitLineData

#given a read object, meta object, and an allele name
#known to already HIT to the read, find the following info
#and put it in a dict 'q. start', 'q. end', 's. start', 's. end', & btop
def getHitInfo(read_result_obj,meta,alleleName):
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
					quern_end=query_start+query_range.length()
					ret_map['q. start']=query_start
					ret_map['q. end']=quern_end-1
					btop=segment_match.btop()
					ret_map['btop']=str(btop)
					ret_map['subject ids']=alleleName
					ret_map['is_inverted']=gls_match.is_inverted()
					ret_map['query id']=read_result_obj.id()
					if(ret_map['is_inverted']):
						ret_map['query id']="reversed|"+read_result_obj.id()
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





def getRegionsObjsFromSegmentCombo(segment_combo):
	print "got a segment combo : ",segment_combo
	print "regions is ",segment_combo.regions
	print "now calling....",
	regions=segment_combo.regions()
	print "the type of regions is ",str(type(regions))
	#regions.new()
	#pprint(regions)
	for region in regions:
		print "region???"
	#for r in range(len(region_list)):
	#	print "Looking at region "+str(r+1)




#given a read result object and meta data
#get the names of the top V,D, and J alleles
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
		#print "The number of segments in this combination is ",len(segment_combination.segments())
		#print "The ids of the segments are ",segment_combination.segments()
		seg_id=0
		for i in segment_combination.segments():
			segment_match = read_result_obj[i]
			#print "got a segment match, id=",i," from combination # ",str(int(s+1))
			seg_type=segTypes[seg_id]
			#print "the seg type is ",seg_type
			for gls_match in segment_match.germline_segments():
				#print "in inner most loop...."
				#print "gls_match=",gls_match
				#print meta[gls_match.gl_segment_].name_
				if(names[seg_type]==None):
					#print "USING ",meta[gls_match.gl_segment_].name_
					names[seg_type]=meta[gls_match.gl_segment_].name_
				else:
					#print "SKIPPING ",meta[gls_match.gl_segment_].name_
					pass
			seg_id+=1
		#print "\n\n\n\n\n"
	#printMap(names)
	return names
	



