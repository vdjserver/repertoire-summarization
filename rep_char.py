#!/usr/bin/env python

import vdjml
from vdjml_igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs
from utils import printMap
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
from vdjml_utils import getTopVDJItems
from Bio import SeqIO
from segment_utils import IncrementMapWrapper

# use a read_result_ojbect return several things:
# 1) the segment counts (a list of 1, 2, or 3 items) , 
# 2) the cdr3_lengths (kabat and imgt) ( a dict() with two keys)
# 3) a clone of the object, but with additional (with additional read information)
def rep_char_read(read_result_obj,meta,organism,imgtdb_obj,read_rec):
	print "In loop analyzing ",read_result_obj.id()


	#prepare an object to return
	return_obj=dict()

	#retrieve the top VDJ
	#topVDJ=getTopVDJItems(read_result_obj,meta)
	#printMap(topVDJ)
	#return_obj['VDJ']=topVDJ

	#retrieve the CDR3 lengths
	cdr3_length_results=CDR3LengthAnalysisVDMLOBJ(read_result_obj,meta,organism,imgtdb_obj,read_rec)
	#print "THE RETURNED RESULTS ARE :"
	#printMap(cdr3_length_results)
	return_obj['cdr3_length_results']=cdr3_length_results


	return return_obj
	








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
	




if (__name__=="__main__"):
	parser=makeParserArgs()
	args = parser.parse_args()
	#args.
	if(args): 
		mf=makeVDJMLDefaultMetaAndFactoryFromArgs(args)
		meta=mf[0]
		fact=mf[1]
		imgtdb_obj=imgt_db("/home/data/DATABASE/01_22_2014/")
		print "a fact is ",fact
		print "the file is ",args.igblast_in[0]
		fasta_reader=SeqIO.parse(open("/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna", "r"), "fasta")
		modes=['kabat','imgt']
		my_cdr3_map=histoMapClass(modes)
		for read_result_obj in scanOutputToVDJML(args.igblast_in[0],fact):
			query_record=fasta_reader.next()
			#print "got a result from ",args.igblast_in[0]
			read_analysis_results=rep_char_read(read_result_obj,meta,"human",imgtdb_obj,query_record)
			cdr3_res=read_analysis_results['cdr3_length_results']
			for mode in modes:
				print "mode and val are ",mode," , ",cdr3_res[mode]
				my_cdr3_map.inc(mode,cdr3_res[mode])			
		#my_cdr3_map.printMaps()
		my_cdr3_map.writeToFile("/dev/stdout")
	else:
		#print "error in args!"
		parser.print_help()


