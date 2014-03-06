#!/usr/bin/env python

import vdjml
from vdjml_igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs
from utils import printMap,get_domain_modes,biopythonTranslate,rev_comp_dna
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
from vdjml_utils import getTopVDJItems,getRegionsObjsFromSegmentCombo,getHitInfo,getProductiveRearrangmentFlag,getVDJServerRegionAlignmentFromLargerVAlignmentPyObj
from Bio import SeqIO
from segment_utils import IncrementMapWrapper,getVRegionsList,getRegionSpecifcCharacterization
from char_utils import getNumberBaseSubsFromBTOP,getNumberIndelsFromBTOP,getIndelMapFromBTOP


global_key_base="vdj_server_ann_"

# use a read_result_ojbect return several things:
# 1) the segment counts (a list of 1, 2, or 3 items) , 
# 2) the cdr3_lengths (kabat and imgt) ( a dict() with two keys)
# 3) a clone of the object, but with additional (with additional read information)
def rep_char_read(read_result_obj,meta,organism,imgtdb_obj,read_rec):
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
	#read_result_obj=vSegmentRegionVDJAnalyse(read_result_obj,meta,organism,imgtdb_obj,read_rec)
	#vSegmentRegionVDJAnalyse(read_result_obj,meta,organism,imgtdb_obj,read_rec)
	read_ann_map=readAnnotate(read_result_obj,meta,organism,imgtdb_obj,read_rec)

	#add some CDR3 information to the ann_amp
	global global_key_base
	for mode in get_domain_modes():
		f=cdr3_length_results[mode+'_from']
		t=cdr3_length_results[mode+'_to']
		if(f!=(-1) or t!=(-1) or cdr3_length_results[mode]==(-1) ):
			qw=str(read_rec.seq)
			if(cdr3_length_results['qry_rev']):
				qw=rev_comp_dna(qw)
			query_cdr3=qw[f-1:t]
			read_ann_map[global_key_base+mode+'_cdr3_na']=query_cdr3
			read_ann_map[global_key_base+mode+'_cdr3_tr']=biopythonTranslate(read_ann_map[global_key_base+mode+'_cdr3_na'])
		else:
			read_ann_map[global_key_base+mode+'_cdr3_na']=""
			read_ann_map[global_key_base+mode+'_cdr3_tr']=read_ann_map[global_key_base+mode+'_cdr3_na']
			

	printMap(read_ann_map,True)
	sys.exit(0)
	return return_obj
	


def readAnnotate(read_result_obj,meta,organism,imgtdb_obj,read_rec):
	print "To annotate a read...."
	topVDJ=getTopVDJItems(read_result_obj,meta)
	annMap=dict()
	for seg in topVDJ:
		annMap['top_'+seg]=topVDJ[seg]
	#print "ann map from segs :"
	#printMap(annMap)
	#getHitInfo(read_result_obj,meta,alleleName):
	hit_infos=dict()
	whole_seq_number_base_subs=0
	whole_seq_number_indels=0
	whole_seq_number_insertions=0
	whole_seq_number_deletions=0
	for seg in topVDJ:
		#print "in second loop...."
		if(topVDJ[seg] is not None):
			#print "to get info...."
			hit_info=getHitInfo(read_result_obj,meta,topVDJ[seg])
			btop=hit_info['btop']
			#print "got btop ",seg," (",topVDJ[seg],") :",btop
			whole_seq_number_base_subs+=getNumberBaseSubsFromBTOP(btop)
			whole_seq_number_indels+=getNumberIndelsFromBTOP(btop)
			indel_map=getIndelMapFromBTOP(btop)
			whole_seq_number_insertions+=indel_map['insertions']
			whole_seq_number_deletions+=indel_map['deletions']
			#print "from ",btop," got ",muts," substitutions"
			#printMap(hit_info)
			#printMap(hit_infos[seg])
		else:
			print "got a none seg!"
	annMap['whole_seq_number_base_subs']=whole_seq_number_base_subs
	annMap['whole_seq_number_indels']=whole_seq_number_indels
	annMap['whole_seq_number_insertions']=whole_seq_number_insertions
	annMap['whole_seq_number_deletions']=whole_seq_number_deletions
	annMap['productive_rearrangement']=getProductiveRearrangmentFlag(read_result_obj,meta,organism,imgtdb_obj)
	#printMap(annMap)
	#VDJSERVER ANNOTATIONS
	mode_list=get_domain_modes()
	for mode in mode_list:
		for region in getVRegionsList():
			#print "Now to analyze region ",region," in mode",mode
			raInfo=getVDJServerRegionAlignmentFromLargerVAlignmentPyObj(read_result_obj,meta,organism,mode,region,imgtdb_obj,False,read_rec)
			if(raInfo!=None):
				#print "got no info!"
				aln=raInfo[0]
				q_aln=aln[0]
				s_aln=aln[1]
				mask=raInfo[1]
				q_start=raInfo[2]
				q_end=raInfo[3]
				reg_ann=getRegionSpecifcCharacterization(s_aln,q_aln,region,mask,mode)
				#print "got ann "
				#printMap(reg_ann)
				key_base=global_key_base+mode+"_"
				for key in reg_ann:
					annMap[key_base+region+"_"+key]=reg_ann[key]
			else:
				#print raInfo
				pass
	#printMap(annMap,True)
	#sys.exit(0)
	return annMap








#using the PYVDJML, add VDJserver specific tags
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
	


def add_rep_char_args_to_parser(parser):
	parser.add_argument('-json_out',type=str,nargs=1,default="/dev/stdout",help="output file for the JSON segment count IMGT hierarchy")
	parser.add_argument('-cdr3_hist_out',type=str,nargs=1,default="/dev/stdout",help="output file for the JSON segment count IMGT hierarchy")
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
		rrw = vdjml.Result_writer(str(args.vdjml_out[0]), meta)
		for read_result_obj in scanOutputToVDJML(args.igblast_in[0],fact):
			#prepare for the iteration and give a possible status message...
			if(read_num>1 and read_num%1000==0):
				print "Processed read",read_num,"..."
			query_record=fasta_reader.next()
			
			#analyze a read's results
			read_analysis_results=rep_char_read(read_result_obj,meta,organism,imgtdb_obj,query_record)

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

			#increment the read number
			read_num+=1

		#write the CDR3 hist	
		my_cdr3_map.writeToFile(cdr3_hist_out_file)

		#write the segment counts
		segment_counter.JSONIFYToFile(args.vdj_db[0],organism,segments_json_out)

	else:
		#print "error in args!"
		parser.print_help()


