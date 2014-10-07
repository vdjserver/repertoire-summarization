#!/usr/bin/env python

from utils import glob_walk,extractAsItemOrFirstFromList
import os
from sample_cdr3 import getBatchIDAndSampleIDFromPath



#for BS in `find */*.fasta` ; do echo "BS IS $BS" ; BID=`echo $BS|tr "/" "\t"|cut -f1` ; echo "BATCH_ID  IS $BID" ; SID=`echo $BS|tr "/" "\t"|cut -f2|grep -Po '^.+\.'|tr -d "."`; echo "SID IS $SID" ; CODE=`awk -v BATCH=$BID -v SAMP=$SID   '{if($1==BATCH && $2==SAMP ) {   print $1 "\t" $2 "\t" $3} }' ../B.S.SBC.Lookup.tsv|cut -f3` ; echo "the code is $CODE" ; TSV="$BS.igblast.out.rc_out.tsv" ; awk -v code=$CODE   -F"\t" '{if($184=="OK") {print code "\t" $21}}' $TSV >> cdr3_sample_barcode_id.OK.tsv  ;        done ;
#add filter so that "None" and blank fall off
#we need to take redundancy into account here..... (we need 3 columsn and the 3rd column is the copy-count)





def isControlBC(bc):
	if(
		bc=="M105" or
		bc=="M106" or
		bc=="N105" or
		bc=="N106"):
		return True
	else:
		return False




def mapToControl(bc):
	if(isControlBC(bc)):
		return "CONTROL"
	else:
		return bc



#Given the input file, return 3 things:
#1) return a sample-barcode -> CDR3 mapping (one -> many) telling, for a given sample-barcode, the set of CDR3s found in it
#2) a CDR3->count dict, telling, given a CDR3, the number of times it was found in the input file
#3) a pair->count file telling, given a [barcode-sample,CDR3] pair, the number of times it was encountered in the input file
def get_bccdr3map_cdr3countmap_paircountmap(bc_cdr3_file):
	#given a barcode, get the set of CDR3s present in it
	bc_cdr3_dict=dict()
	reader=open(bc_cdr3_file,'r')
	#dict to count number of CDR3 appearances
	cdr3_count=dict()
	print "Now mapping codes to a set of CDRs...so that given a code, a set of CDR3s may be obtained."
	print "Also counting the CDR3s....."
	print "Also counting sample-barcode CDR3 paired counts!"
	bc_cdr3_pair_count=dict()
	for line in reader:
		temp_line=line.strip()
		#print "Read line '"+temp_line+"' for mapping a code to a set of CDR3s...."
		pieces=temp_line.split('\t')
		code=pieces[0]
		c=pieces[1]
		#will look for 3 columns soon
		pair=str(temp_line)
		if(not(pair in bc_cdr3_pair_count)):
			bc_cdr3_pair_count[pair]=0
		bc_cdr3_pair_count[pair]+=1
		if(not(c in cdr3_count)):
			cdr3_count[c]=0
		cdr3_count[c]+=1
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
	parser.add_argument('cdr3_sample_barcode_OK_tsv',type=str,nargs=1,help="path to a file of occurences of sample-barcode resolved CDR3s conditioned on OK AGS status ; each line has [BARCODE-SAMPLE] [tab] [CDR3 NUCL. ACID SEQ] ; these are extracted from TSVs and a barcode-sample/batch-sample lookup table conditioned on the OK AGS status ; a line may appear multiple times.  these are not counts ; they are instances'
	parser.add_argument('cdr3_whittle_down',type=str,nargs=1,help="the output file of ordered triples of [sample barcode][tab][CDR3 nucl. acid][tab][percentage]") 
	args=parser.parse_args()
	if(args):
		print "success"
	else:
		print "failure"









