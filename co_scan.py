#!/usr/bin/env python

from utils import glob_walk,extractAsItemOrFirstFromList
import os
from sample_cdr3 import getBatchIDAndSampleIDFromPath
import argparse
import re


#for BS in `find */*.fasta` ; do echo "BS IS $BS" ; BID=`echo $BS|tr "/" "\t"|cut -f1` ; echo "BATCH_ID  IS $BID" ; SID=`echo $BS|tr "/" "\t"|cut -f2|grep -Po '^.+\.'|tr -d "."`; echo "SID IS $SID" ; CODE=`awk -v BATCH=$BID -v SAMP=$SID   '{if($1==BATCH && $2==SAMP ) {   print $1 "\t" $2 "\t" $3} }' ../B.S.SBC.Lookup.tsv|cut -f3` ; echo "the code is $CODE" ; TSV="$BS.igblast.out.rc_out.tsv" ; awk -v code=$CODE   -F"\t" '{if($184=="OK") {print code "\t" $21}}' $TSV >> cdr3_sample_barcode_id.OK.tsv  ;        done ;
#add filter so that "None" and blank fall off
#we need to take redundancy into account here..... (we need 3 columsn and the 3rd column is the copy-count)




#is a given sample-barcode a controL?  return true or false
#M105,M106,N105,N106 are classified as CONTROLS all others not.
def isControlBC(bc):
	if(
		bc=="M105" or
		bc=="M106" or
		bc=="N105" or
		bc=="N106"):
		return True
	else:
		return False



#given a sample barcode, map it to itself (if it's not control) and to CONTROL othertise
def mapToControl(bc):
	if(isControlBC(bc)):
		return "CONTROL"
	else:
		return bc





#utility to read the BATCH_SAMPLE_BARCODE file and return a lookup datastructure from it
def obtainBatchSampleBSMapFromLookupFile(lookup_file_path):
	lu=dict()
	reader=open(lookup_file_path,'r')
	for line in reader:
		temp_line=line.strip()
		pieces=temp_line.split('\t')
		if(not(len(pieces)==3)):
			err_msg="Error, got "+str(len(pieces))+" tab-separated values in file "+lookup_file_path+" but expected 3 : BATCH[TAB]SAMPLE[TAB]SUBJECT_BARCODE !"
			raise Exception(err_msg)
			sys.exit(0)
		batch=pieces[0]
		sample=pieces[1]
		bs=pieces[2]
		if(not(batch in lu)):
			lu[batch]=dict()
		lu[batch][sample]=bs
	reader.close()
	return lu
	#/home/esalina2/diogenix_2014_contract/B.S.SBC.Lookup.tsv




def scanTSVAppendingTripletBSAndCDR3NAAndCountToFile(tsv_file_base,file_to_append_to,batch_sample_barcode_lookup_file):
	b_s_bs_map=obtainBatchSampleBSMapFromLookupFile(batch_sample_barcode_lookup_file)
	print "Loaded lookup from ",batch_sample_barcode_lookup_file
	tsv_glob=tsv_file_base+"/*/*out.tsv"
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
				cdr3_na=pieces[20]
				status=pieces[183]
				count=int(pieces[len(pieces)-1])
				#print "in file ",tsv,"on line=",line_num," got cdr3=",cdr3_na," status=",status," and count=",str(count)
				if(status=="OK"):
					#got OK status!
					if(re.match(r'^[ACGT]+$',cdr3_na)):
						#got a good CDR3!
						line_to_write=control_mapped_subject_barcode+"\t"+cdr3_na+"\t"+str(count)
						writer.write(line_to_write+"\n")
			line_num+=1
		reader.close()
	writer.close()
	
		




#Given the input file (BS_CDR3NA_COUNT) , return 3 things:
#1) return a sample-barcode -> CDR3 mapping (one -> many) telling, for a given sample-barcode, the set of CDR3s found in it (this is an indicator-style type structure)
#2) a CDR3->count dict, telling, given a CDR3, the number of times it was found in the input file (this keeps a weighted count ; weighted based on input counts from 3rd column)
#3) a pair->count file telling, given a [barcode-sample,CDR3] pair, the number of times it was encountered in the input file (this keeps a weighted count ; weighted based on input counts from third column)
def get_bccdr3map_cdr3countmap_paircountmap(bc_cdr3_count_file):
	#given a barcode, get the set of CDR3s present in it
	bc_cdr3_dict=dict()
	reader=open(bc_cdr3_count_file,'r')
	#dict to count number of CDR3 appearances
	cdr3_count=dict()
	print "Now mapping codes to a set of CDRs...so that given a code, a set of CDR3s may be obtained."
	print "Also counting the CDR3s....."
	print "Also counting sample-barcode CDR3 paired counts! (weighted by redundancy counts)"
	bc_cdr3_pair_count=dict()
	for line in reader:
		temp_line=line.strip()
		#print "Read line '"+temp_line+"' for mapping a code to a set of CDR3s...."
		pieces=temp_line.split('\t')
		code=pieces[0]
		c=pieces[1]
		this_INPUT_count=int(pieces[2])
		#will look for 3 columns soon
		pair=code+"\t"+c
		if(not(pair in bc_cdr3_pair_count)):
			bc_cdr3_pair_count[pair]=0
		bc_cdr3_pair_count[pair]+=this_INPUT_count
		if(not(c in cdr3_count)):
			cdr3_count[c]=0
		cdr3_count[c]+=this_INPUT_count
		if(code in bc_cdr3_dict):
			bc_cdr3_dict[mapToControl(code)].add(c)
		else:
			bc_cdr3_dict[mapToControl(code)]=set()
			bc_cdr3_dict[mapToControl(code)].add(c)
		#print "After this line, code "+code+" maps to ",bc_cdr3_dict[code]
	reader.close()
	return [bc_cdr3_dict,cdr3_count,bc_cdr3_pair_count]







#given
#1) a barcode->CDR3 (one->many) mapping and
#2) a CDR3 count map (one-one)
#obtain/return a CDR3->barcode (one-many) mapping so that, given a CDR3, a set of barcode-samples it falls in may be known/obrtained
def create_cdr3_bcset_mapping(bc_cdr3_dict,cdr3_count):
	#given a CDR3 get the set of barcodes it's present in 
	cdr3_bc_dict=dict()
	print "Now mapping CDR3s to codes so that, given a CDR3, a set of codes it falls in may be obtained...."
	for cdr3 in cdr3_count:
		temp_cdr3=cdr3.strip()
		#print "Read a cdr3 ",temp_cdr3
		if(not(temp_cdr3 in cdr3_bc_dict)):
			cdr3_bc_dict[temp_cdr3]=set()
		for barcode in bc_cdr3_dict:
			if(temp_cdr3 in bc_cdr3_dict[barcode]):
				cdr3_bc_dict[temp_cdr3].add(barcode)			
		if(not(temp_cdr3 in cdr3_bc_dict)):
			raise Exception(str(temp_cdr3)+"never got assigned to a barcode!")
		#print "After the line, the set of codes it is associated with is ",cdr3_bc_dict[temp_cdr3]
	return cdr3_bc_dict




#given a CDR3->sample-barcode (one-many) mapping, obtain a set of 'blacklisted' CDR3s
#which, by definition, belong in more than one sample-barcode
def getCDR3sBlackListed(cdr3_bc_dict):
	#NOW, scan each CDR3 and print a message if it is in more than one barcode!
	#simultaneously create a 'black-listed' set of CDR3s.....
	cdr3_black_list_set=set()
	num_cdr3_in_multiple_code=0
	print "now in final scan...."
	for cdr3 in cdr3_bc_dict:
		bc_set=cdr3_bc_dict[cdr3]
		if(len(bc_set)>1):
			cdr3_black_list_set.add(cdr3)
			print "CDR3 ",cdr3," belongs in more than one barcode : ",bc_set
			num_cdr3_in_multiple_code+=1
	print "The number of CDR3s being found associated with more than one barcode is ",num_cdr3_in_multiple_code
	return cdr3_black_list_set





#generate the whittle down blacklist file of blacklisted CDR3s with sample-barcodes and percent given:
#1) a barcode-sample CDR3 -> count mapping (one-one) and
#2) a set of blacklisted CDR3s that are in two or more sample-barcodes
#nothing returned, but the output file is written!
def generate_whittle_down(bc_cdr3_pair_count,cdr3_black_list_set,whittled_down_output):
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
	parser=argparse.ArgumentParser(description='Perfrom subject-barcode CDR3 cross-over initial scan writing output of sample-barcode, CDR3 and cross-over percent.  First, scan for CDR3s being contained in multiple barcode-samples.  If a CDR3 is contained in exactly 1 barcode-sample it it NOT blacklisted and does NOT appear in the output.  If a CDR3 DOES appear in two or more barcode-samples then it it \'blacklisted\' and DOES appear in the output.  For each CDR3 blacklisted, for each barcode-sample it appears in, there will be a corresponding line in the output.  Moreover, each of those lines in the output is accompanied by a percentage indicating the relative percent that the CDR3 appears in the sample-barcode as a percentage of all the sample-barcodes it appears in.  NOTE : any one of the 4 : M105,M106,N105,N106 subject-barcodes will be interpeted as CONTROL!')
	#parser.add_argument('cdr3_sample_barcode_OK_tsv',type=str,nargs=1,help="path to a file (which is created by this script) of occurences of sample-barcode resolved CDR3s conditioned on OK AGS status ; each line has [BARCODE-SAMPLE] [tab] [CDR3 NUCL. ACID SEQ] [tab] [COUNT]; these are extracted from TSVs and a barcode-sample/batch-sample lookup table conditioned on the OK AGS status ; a line may appear multiple times.  these are not counts ; they are instances")
	#parser.add_argument('cdr3_whittle_down',type=str,nargs=1,help="the output file of ordered triples of [sample barcode][tab][CDR3 nucl. acid][tab][percentage]") 
	parser.add_argument('tsv_base',type=str,nargs=1,help="base directory for TSVs where glog base/*/*out.tsv is used to scan for CDR3s,TSVs, and counts.  In this directory the file cdr3_black_list_whittle_pct.tsv (which is)  will be created by the script as well as cdr3_sample_barcode_id.OK.CONTROL.tsv (which is created by writing instances of triples of SUBJECT_BARCODE[TAB]CDR3_NA(not 'None')[TAB]Count from scanning TSVs)")
	parser.add_argument('b_s_lookup_path',type=str,nargs=1,help="path to the BATCH_SAMPLE_BARCODESAMPLE lookup file")
	args=parser.parse_args()
	if(args):
		#success!
		import os
		tsv_base=extractAsItemOrFirstFromList(args.tsv_base)
		bs_lookup=extractAsItemOrFirstFromList(args.b_s_lookup_path)
		err_state=False
		if(os.path.exists(tsv_base) and os.path.isdir(tsv_base)):
			#good
			if(os.path.isfile(bs_lookup) and os.path.exists(bs_lookup)):
				#good
				file_to_append_to=tsv_base+"/cdr3_sample_barcode_id.OK.CONTROL.tsv"
				if(not(os.path.exists(file_to_append_to))):
					whittled_down_output=tsv_base+"/cdr3_black_list_whittle_pct.tsv"
					if(not(os.path.exists(whittled_down_output))):
						#good
						#good, can't overwrite a non-existing file!! now proceed to generate it!
						scanTSVAppendingTripletBSAndCDR3NAAndCountToFile(tsv_base,file_to_append_to,bs_lookup)
						ontained_counting_info=get_bccdr3map_cdr3countmap_paircountmap(file_to_append_to)
						bc_cdr3_dict=ontained_counting_info[0] #map, given a BS, gives a set of CDR3s found in it (not weighted, just an indicator)
						cdr3_count=ontained_counting_info[1] # map, given a CDR3 NA seq, get how many times it appeared (WEIGHTED)
						bc_cdr3_pair_count=ontained_counting_info[2] #WEIGHTED count of BC/CDR3 pairings
						#given  weighted counts, obtain a mapping where you give it a CDR3, 
						#it tells/gives you a SET of sample-barcodes that the CDR3 falls in
						cdr3_bc_dict=create_cdr3_bcset_mapping(bc_cdr3_dict,cdr3_count)
						#given the CDR3->BS mapping obtain a set of 'blacklisted CDR3s (appearing in more than one sample-barcode)
						cdr3_black_list_set=getCDR3sBlackListed(cdr3_bc_dict)					
						generate_whittle_down(bc_cdr3_pair_count,cdr3_black_list_set,whittled_down_output)
					else:
						print "File ",whittled_down_output," found! Must first delete it to proceed!"
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
		#print "failure"
		parser.print_help()




