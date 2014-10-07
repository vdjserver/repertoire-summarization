#!/usr/bin/env python



#make a map from a base CSV file
def create_read_id_count_map(csv_file_path):
	reader=open(csv_file_path)
	effective_line_num=0
	read_count_map=dict()
	for line in reader:
		temp_line=line.strip()
		if(temp_line[0]=="#"):
			#ignore lines starting with #
			pass
		else:
			if(effective_line_num==0):
				#ignore header 
				pass
			else:
				pieces=temp_line.split('\t')
				#print "Pieces read in are ",pieces
				read_id=pieces[0]
				count_str=pieces[1]
				count_val=int(count_str)
				read_count_map[read_id]=count_val
			effective_line_num+=1
	reader.close()
	return read_count_map
		


#make list of ordered READ IDs from the CSV
def getOrderedReadIDListFromCSV(csv_file_path):
	reader=open(csv_file_path)
	effective_line_num=0
	eff_list=list()
	for line in reader:
		temp_line=line.strip()
		if(temp_line[0]=="#"):
			#ignore lines starting with #
			pass
		else:
			if(effective_line_num==0):
				#ignore header 
				pass
			else:
				pieces=temp_line.split('\t')
				read_id=pieces[0]
				count_str=pieces[1]
				eff_list.append(read_id)
			effective_line_num+=1
	reader.close()
	return eff_list




#given a fasta path, generate the dup.csv path
def makeCSVPathFromFastaPath(fp):
	if(fp.endswith(".fasta")):
		pattern="\.fasta$"
		#Sample00001dup.csv
		#Sample00001.fasta
		replacement="dup.csv"
		import re
		fp=re.sub(pattern,replacement,fp)
		return fp
	else:
		#not a .fasta?
		return None

#from a TSV, generate the FASTA
def makeFASTAPathFromTSV(tsv_path):
	import re
	fasta=re.sub(".fasta[^/]*$",".fasta",tsv_path)
	return fasta


#from TSV generate CSV PATH
def fromTSVGenerateCSVPath(tsv_path):
	return makeCSVPathFromFastaPath(makeFASTAPathFromTSV(tsv_path))



#add the column to the TSV
def addExtraCountColToTSV(tsv_path):
	csv_path=fromTSVGenerateCSVPath(tsv_path)
	csv_read_ids=getOrderedReadIDListFromCSV(csv_path)
	#print "The length of csv_read_ids is ",len(csv_read_ids)
	read_id_to_count_map=create_read_id_count_map(csv_path)
	new_path=tsv_path+".with_counts.tsv"
	#first read in the file and write out the new one with countds
	reader=open(tsv_path,'r')
	writer=open(new_path,'w')
	ln=0
	for line in reader:
		temp_line=line.strip()
		if(ln==0):
			#header
			new_col_name="Count"
			pieces=temp_line.split('\t')
			pieces.append(new_col_name)
			new_line="\t".join(pieces)
			writer.write(new_line+"\n")
			#print "Just wrote header"
		else:
			#data
			pieces=temp_line.split('\t')
			read_id=pieces[1]
			#print "retreive TSV id is ",read_id
			ordered_id=csv_read_ids[ln-1]
			#print "retrived CSV id is ",ordered_id
			err_msg="Error in file CSV="+csv_path+" TSV="+tsv_path+" at read ID(TSV)="+read_id+" and read ID(CSV)="+ordered_id+
			if(not(ordered_id==read_id)):
				err_msg+=" Mismatch (line="+str(ln)+")"
				raise Exception(err_msg)
				return None
			if(not(read_id in read_id_to_count_map)):
				err_msg+=" READ_ID="+read_id+", not found in map to obtain count value!
			
			count_val=read_id_to_count_map[read_id]
			pieces.append(str(count_val))
			writer.write("\t".join(pieces)+"\n")
			#print "just wrote at ln=",ln
		ln+=1
	reader.close()
	#print "read from ",tsv_path," and wrote to ",new_path
	writer.close()
	import os
	import shutil
	#move the newly created file to write over 
	#the old file, thus adding the new column
	os.remove(tsv_path)
	shutil.move(new_path,tsv_path)
	


	



TSV_PATH="/home/esalina2/diogenix_2014_contract/set2014_su_counts/S1401136_Run1/Sample00001.fasta.igblast.out.rc_out.tsv"
addExtraCountColToTSV(TSV_PATH)


