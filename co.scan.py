#!/usr/bin/env python

from utils import glob_walk
import os
from sample_cdr3 import getBatchIDAndSampleIDFromPath
from Bio import SeqIO



#for BS in `find */*.fasta` ; do echo "BS IS $BS" ; BID=`echo $BS|tr "/" "\t"|cut -f1` ; echo "BATCH_ID  IS $BID" ; SID=`echo $BS|tr "/" "\t"|cut -f2|grep -Po '^.+\.'|tr -d "."`; echo "SID IS $SID" ; CODE=`awk -v BATCH=$BID -v SAMP=$SID   '{if($1==BATCH && $2==SAMP ) {   print $1 "\t" $2 "\t" $3} }' ../B.S.SBC.Lookup.tsv|cut -f3` ; echo "the code is $CODE" ; TSV="$BS.igblast.out.rc_out.tsv" ; awk -v code=$CODE   -F"\t" '{if($184=="OK") {print code "\t" $21}}' $TSV >> cdr3_sample_barcode_id.OK.tsv  ;        done ;

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


cdr3_uniq_file="/home/esalina2/diogenix_2014_contract/set2014_r0_su/uniq.cdr3s.txt"
#barcode_cdr3_file="/home/esalina2/diogenix_2014_contract/set2014_r0_su/cdr3_sample_barcode_id.OK.tsv"
barcode_cdr3_file="/home/esalina2/diogenix_2014_contract/set2014_r0_su/cdr3_sample_barcode_id.OK.CONTROL_REN.tsv"

#given a barcode, get the set of CDR3s present in it
bc_cdr3_dict=dict()
reader=open(barcode_cdr3_file,'r')
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



#given a CDR3 get the set of barcodes it's present in 
cdr3_bc_dict=dict()
cdr3_reader=open(cdr3_uniq_file,'r')
print "Now mapping CDR3s to codes so that, given a CDR3, a set of codes it falls in may be obtained...."
for cdr3 in cdr3_reader:
	temp_cdr3=cdr3.strip()
	#print "Read a cdr3 ",temp_cdr3
	if(not(temp_cdr3 in cdr3_bc_dict)):
		cdr3_bc_dict[temp_cdr3]=set()
	for barcode in bc_cdr3_dict:
		if(temp_cdr3 in bc_cdr3_dict[barcode]):
			cdr3_bc_dict[temp_cdr3].add(barcode)			
	if(not(temp_cdr3 in cdr3_bc_dict)):
		print temp_cdr3,"never got assigned to a barcode!"
	#print "After the line, the set of codes it is associated with is ",cdr3_bc_dict[temp_cdr3]
cdr3_reader.close()


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


whittled_down="/home/esalina2/diogenix_2014_contract/set2014_r0_su/cdr3_black_list_whittle_pct.tsv"
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
		





#Now rewrite the data, but skipping over reads who generated the black_listed CDR3s.....
#co_filtered_base="/home/esalina2/diogenix_2014_contract/set2014_r0_su_post_CO_scan"
#if(os.path.exists(co_filtered_base)):
#	print "output base ",co_filtered_base," exists, so not writing to it and aborting!"
#else:
#	os.makedirs(co_filtered_base)
#	black_listed_readID_set=set()
#	white_listed_readID_set=set()
#	tsvs=glob_walk("/home/esalina2/diogenix_2014_contract/set2014_r0_su/*/*out.tsv")
#	#find black/white READ IDs depending on black_listed CDR3 status
#	for tsv in tsvs:
#		#print "In TSV collecting reads to blacklist...."
#		tsv_reader=open(tsv,'r')
#		for line in tsv_reader:
#			temp_line=line
#			temp_line=temp_line.strip()
#			pieces=temp_line.split('\t')
#			this_cdr3=pieces[20]
#			this_readid=pieces[1]
#			this_status=pieces[183]
#			#print "\n\n\nfile ",tsv
#			#print "cdr3=",this_cdr3
#			#print "this readid=",this_readid
#			#print "this status=",this_status
#			if((this_cdr3 in cdr3_black_list_set) or not( (this_status=="OK"))):
#				black_listed_readID_set.add(this_readid)
#			else:
#				white_listed_readID_set.add(this_readid)
#	fastas=glob_walk("/home/esalina2/diogenix_2014_contract/set2014_r0_su/*/*.fasta")
#	for fasta in fastas:
#		#print "Looking at fasta ",fasta
#		bs=getBatchIDAndSampleIDFromPath(fasta)
#		batch=bs[0]
#		sample=bs[1]
#		print "bs from ",fasta," is ",bs
#		target_fasta_dir=co_filtered_base+"/"+batch
#		if(not(os.path.exists(target_fasta_dir))):
#			os.makedirs(target_fasta_dir)
#		target_fasta=target_fasta_dir+"/"+sample+".fasta"
#		print "The filtered target is ",target_fasta
#		to_write_list=list()
#		for fasta_record in  SeqIO.parse(open(fasta, "r"), "fasta"):
#			fasta_rec_id=fasta_record.id
#			fasta_rec_seq=fasta_record.seq
#			if(fasta_rec_id in white_listed_readID_set):
#				to_write_list.append(fasta_record)
#		SeqIO.write(to_write_list,target_fasta,"fasta")
#
#




