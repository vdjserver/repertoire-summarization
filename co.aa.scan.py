#!/usr/bin/env python

from utils import glob_walk,extractAsItemOrFirstFromList
import os
from sample_cdr3 import getBatchIDAndSampleIDFromPath
import argparse
import re
from co_scan import isControlBC,mapToControl,obtainBatchSampleBSMapFromLookupFile,get_bccdr3map_cdr3countmap_paircountmap,create_cdr3_bcset_mapping,generate_whittle_down,getCDR3sBlackListed




def scanTSVAppendingTripletBSAndCDR3AAAndCountToFile(tsv_base,file_to_append_to,bs_lookup):
	b_s_bs_map=obtainBatchSampleBSMapFromLookupFile(bs_lookup)
	print "Loaded lookup from ",bs_lookup
	tsv_glob=tsv_base+"/*/*out.tsv"
	writer=open(file_to_append_to,'w')
	for tsv in glob_walk(tsv_glob):
		print "Now to scan TSV ",tsv," to append SB-CDR3-COUNT to file ",file_to_append_to
		batch_and_sample=getBatchIDAndSampleIDFromPath(tsv)
		print "Its batch and sample : ",batch_and_sample
		batch=batch_and_sample[0]
		sample=batch_and_sample[1]
		subject_barcode=b_s_bs_map[batch][sample]
		control_mapped_subject_barcode=mapToControl(subject_barcode)
		reader=open(tsv,'r')
		line_num=0
		for line in reader:
			if(line_num>=1):
				temp_line=line.strip()
				pieces=temp_line.split('\t')
				cdr3_aa=pieces[18] #use AA instead of NA
				status=pieces[183]
				count=int(pieces[len(pieces)-1])
				#print "in file ",tsv,"on line=",line_num," got cdr3=",cdr3_na," status=",status," and count=",str(count)
				if(status=="OK"):
					#got OK status!
					if(re.match(r'^[GPAVLIMCFYWHKRQNEDST]+$',cdr3_aa)):
						#got a good CDR3!
						line_to_write=control_mapped_subject_barcode+"\t"+cdr3_aa+"\t"+str(count)
						writer.write(line_to_write+"\n")
						#print line_to_write
					else:
						#not incorporating because had non-AA! (perhaps X)?
						pass
			line_num+=1
		reader.close()
	writer.close()	






#generate the whittle down blacklist file of blacklisted CDR3s with sample-barcodes and percent given:
#1) a barcode-sample CDR3 -> count mapping (one-one) and
#2) a set of blacklisted CDR3s that are in two or more sample-barcodes
#nothing returned, but the output file is written!
def generate_AA_whittle_down(bc_cdr3_pair_count,cdr3_black_list_set,whittled_down_output,cdr3_count):
	whittled_down=whittled_down_output
	whittled_writer=open(whittled_down,'w')
	#write whittled down from the pair counts writing percentages instead of counts
	#but only for CDR3s that are on the blacklist
	for pair in bc_cdr3_pair_count:
		pair_pieces=pair.split('\t')
		#print "Pair_pieces are ",pair_pieces
		pair_bc=pair_pieces[0]
		pair_c=pair_pieces[1]
		total=cdr3_count[pair_c]
		pair_count=bc_cdr3_pair_count[pair]
		if(pair_c in cdr3_black_list_set):
			ratio=float(float(pair_count)/float(total))
			line_to_write=pair_bc+"\t"+pair_c+"\t"+str(ratio)
			whittled_writer.write(line_to_write+"\n")
	whittled_writer.close()
		





if (__name__=="__main__"):
	print "in main"
	parser=argparse.ArgumentParser("AA-based CO")
	parser.add_argument('tsv_base',type=str,nargs=1,help="base directory for TSVs where glog base/*/*out.tsv is used to scan for CDR3s,TSVs, and counts.  In this directory the file cdr3_black_list_whittle_pct.tsv (which is)  will be created by the script as well as cdr3_sample_barcode_id.OK.CONTROL.tsv (which is created by writing instances of triples of SUBJECT_BARCODE[TAB]CDR3_AA(not 'None')[TAB]weighted count from scanning TSVs)")
	parser.add_argument('b_s_lookup_path',type=str,nargs=1,help="path to the BATCH_SAMPLE_BARCODESAMPLE lookup file")
	args=parser.parse_args()
	if(args):
		import os
		tsv_base=extractAsItemOrFirstFromList(args.tsv_base)
		bs_lookup=extractAsItemOrFirstFromList(args.b_s_lookup_path)
		err_state=False
		if(os.path.exists(tsv_base) and os.path.isdir(tsv_base)):
			#good
			if(os.path.isfile(bs_lookup) and os.path.exists(bs_lookup)):
				#good
				file_to_append_to=tsv_base+"/cdr3_sample_barcode_id.OK.CONTROL.tsv"
				print "Using "+file_to_append_to+" for gathering CDR3 weighted count data (with 1) control-mapped subject-barcode,2) CDR3 AA, and 3) weighted count per line)"
				if(not(os.path.exists(file_to_append_to))):
					whittled_down_output=tsv_base+"/cdr3_black_list_whittle_pct.tsv"
					if(not(os.path.exists(whittled_down_output))):
						scanTSVAppendingTripletBSAndCDR3AAAndCountToFile(tsv_base,file_to_append_to,bs_lookup)
						print "File  ",file_to_append_to," created...now scanning it ..."
						ontained_counting_info=get_bccdr3map_cdr3countmap_paircountmap(file_to_append_to)
						bc_cdr3_dict=ontained_counting_info[0] #map, given a BS, gives a set of CDR3s found in it (not weighted, just an indicator)
						cdr3_count=ontained_counting_info[1] # map, given a CDR3 NA seq, get how many times it appeared (WEIGHTED)
						bc_cdr3_pair_count=ontained_counting_info[2] #WEIGHTED count of BC/CDR3 pairings
						print "Now acquired a BS(barcode-sample)->CDR3  set map	(telling which CDR3s are associated with a given BS)"
						print "Now acquired a CDR3 count map (given a CDR3 seq, how many times it appeared (weighted)"
						print "Now acquired a BS-CDR3 pair map (given a BS-CDR3 (as a pair), how many times that pair appeared (weighted))"
						#given  weighted counts, obtain a mapping where you give it a CDR3, 
						#it tells/gives you a SET of sample-barcodes that the CDR3 falls in
						cdr3_bc_dict=create_cdr3_bcset_mapping(bc_cdr3_dict,cdr3_count)
						print "Now obtained a CDR3->BS dict (given a CDR3, obtain a set of BS that it appears in)"
						#given the CDR3->BS mapping obtain a set of 'blacklisted CDR3s (appearing in more than one sample-barcode)
						cdr3_black_list_set=getCDR3sBlackListed(cdr3_bc_dict)
						print "From that, now obtained a black-listed set of CDR3s...\n"
						print "To write to file ",whittled_down_output," CDR3s and their distribution among BS!"
						generate_AA_whittle_down(bc_cdr3_pair_count,cdr3_black_list_set,whittled_down_output,cdr3_count)

					else:
						print "File ",whittled_down_output," found!  must delete it first to proceed!"
						parser.print_help()
				else:
					print "File ",file_to_append_to," found!  must delete it first to proceed!"
					parser.print_help()
			else:
				print "BATCH/SAMPLE barcode lookup file ",bs_lookup," not found!"
				parser.print_help()
		else:
			print "Directory ",tsv_base," not found!"
			parser.print_help()
	else:
		parser.print_help()





