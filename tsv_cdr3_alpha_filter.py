#!/usr/bin/env python

from utils import glob_walk
import os
from sample_cdr3 import getBatchIDAndSampleIDFromPath
from Bio import SeqIO
from ags_mgr import ags_manager


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




def niceAGS(arr):
	if(arr==None):
		return "None\t0.0\t0.0"
	ags_str=str(arr[0])+"\t"+str(arr[1])+"\t"+str(arr[2])
	return ags_str


#define inputs
base_dir="/home/esalina2/diogenix_2014_contract/set2014_r0_su"
cdr3_black_list_whittle="/home/esalina2/diogenix_2014_contract/set2014_r0_su/cdr3_black_list_whittle_pct.tsv"
alpha=0.95
BS_SBC_lookup="/home/esalina2/diogenix_2014_contract/B.S.SBC.Lookup.tsv"


#create batch->sample->barcode lookup
b_s_b_lookup=dict()
#esalina2@eddiecomp:~/diogenix_2014_contract$ head B.S.SBC.Lookup.tsv 
#S1401136_Run1	Sample00001	3SB27
#S1401136_Run1	Sample00002	N106
#S1401136_Run1	Sample00003	4UKMS
#S1401136_Run1	Sample00005	GIFTL
#S1401136_Run1	Sample00006	XURET
treader=open(BS_SBC_lookup,'r')
for line in treader:
	temp_line=line.strip()
	pieces=temp_line.split('\t')
	batch=pieces[0]
	sample=pieces[1]
	bc=pieces[2]
	if(not(batch in b_s_b_lookup)):
		b_s_b_lookup[batch]=dict()
	if(not(sample in b_s_b_lookup[batch])):
		b_s_b_lookup[batch][sample]=mapToControl(bc)


#create the new black list
new_black_list_cdr3_bc=dict()
w_reader=open(cdr3_black_list_whittle,'r')



for line in w_reader:
	temp_line=line.strip()
	pieces=temp_line.split('\t')
	sbc=pieces[0]
	cdr3=pieces[1]
	pct=float(pieces[2])
	if(pct>=alpha):
		#don't add to blacklist!
		pass
	elif(not(cdr3=="None")):
		if(not(cdr3 in new_black_list_cdr3_bc)):
			new_black_list_cdr3_bc[cdr3]=set()
		else:
			pass
		new_black_list_cdr3_bc[cdr3].add(sbc)




#define the output base
output_filtered_base="/home/esalina2/diogenix_2014_contract/cdr3_alpha_"+str(alpha)+"_CO_filter/"
if(not(os.path.exists(output_filtered_base))):
	os.makedirs(output_filtered_base)



#now write out the post-filtered TSVs
tsvs=glob_walk(base_dir+"/*/*out.tsv")
for tsv in tsvs:
	print "looking at TSV ",tsv
	#print "BS : ",
	b_s=getBatchIDAndSampleIDFromPath(tsv)
	batch=b_s[0]
	sample=b_s[1]
	subject_barcode=b_s_b_lookup[batch][sample]
	output_dir=output_filtered_base+batch+"/"
	if(not(os.path.exists(output_dir))):
		os.makedirs(output_dir)
	output_path=output_dir+os.path.basename(tsv)
	print "output path is ",output_path
	tsv_writer=open(output_path,'w')
	num_reads_OK=0
	num_reads_OK_None_CDR3=0
	num_reads_OK_Have_CDR3_ON_black_list=0
	num_written=0
	line_num=0
	tsv_reader=open(tsv,'r')
	myAGSMgr=ags_manager("test")
	for line in tsv_reader:
		temp_line=line.strip()
		if(line_num==0):
			#just write out the header!
			tsv_writer.write(line)
			pass
		else:
			pieces=temp_line.split('\t')
			cdr3=pieces[20]
			status=pieces[183]
			rm_list_str=pieces[185]
			rm_list=eval(rm_list_str)
			#print "the rm_list is ",rm_list
			if(status=="OK"):
				num_reads_OK+=1
				if(cdr3=="None"):
					num_reads_OK_None_CDR3+=1
				else:
					#the read is OK and the CDR3 isn't none.
					#see if it's blacklisted for this barcode!
					if(not(cdr3 in new_black_list_cdr3_bc)):
						#write out the line, not blacklisted for any barcode
						num_written+=1
						tsv_writer.write(line)
						for rm in rm_list:
							myAGSMgr.receive_numbered_mut(rm)
						pass
					else:
						if(subject_barcode in new_black_list_cdr3_bc[cdr3]):
							#it's blacklisted for this barcode!
							num_reads_OK_Have_CDR3_ON_black_list+=1
						else:
							#write out the line cause the subject-barcode isn't blacklisted for this CDR3
							tsv_writer.write(line)
							num_written+=1
							for rm in rm_list:
								myAGSMgr.receive_numbered_mut(rm)
							pass
			else:
				pass
		line_num+=1
	ags_6=myAGSMgr.compute_ags("AGS6")
	ags_5=myAGSMgr.compute_ags("AGS5")
	print "REPORT\t"+batch+"\t"+sample+"\t"+str(num_reads_OK)+"\t"+str(num_reads_OK_None_CDR3)+"\t"+str(num_reads_OK_Have_CDR3_ON_black_list)+"\t"+str(num_written)+"\t"+niceAGS(ags_6)+"\t"+niceAGS(ags_5)







