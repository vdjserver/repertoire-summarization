#!/usr/bin/env python

import re
import os
from utils import *
from Bio.Blast import NCBIXML


#ROUTINES FOR MAKING LOOKUP TABLES OF KABAT MODE




#scan blastx output
#use motifs from http://www.bioinf.org.uk/abs/#cdrid
#to identify CDR3 end for KABAT
def writeKabatJCDR3End(k,o):
	reader=open(k,'r')
	blast_records = NCBIXML.parse(reader)
	#regex for CDR3 for heavy chains
	heavy_re=re.compile(r'WG.G')
	#regex for CDR3 for light chains
	light_re=re.compile(r'FG.G')
	cdr_map=dict()
	writer=open(o,'w')
	for record in blast_records:
		print "\n\n\n******************"
		print "query : ",record.query
		query_name=record.query
		if(len(record.alignments)<1):
			print "TOO FEW ALIGNMENTS",query_name
			cdr_map[query_name]=(-1)
			continue
		alignment=record.alignments[0]
		print "An alignment found for ",str(query_name)," has length ",alignment.length
		if(len(alignment.hsps)<1):
			cdr_map[query_name]=(-1)
			print "TOO FEW HSPs",query_name
			continue
		hsps=alignment.hsps
		hsp=alignment.hsps[0]
		query_trans=hsp.query
		query_start=hsp.query_start
		print "HSP",str(hsp),"\n"
		#use the regex for the chain type
		if(query_name.startswith("IGH")):
			myre=heavy_re
		else:
			myre=light_re
		sr=re.search(myre,query_trans)
		if(sr):
			search_position=sr.start()
			print "A search position was found using regex ",str(myre.pattern)
			print "The searched_position (AA space) is ",search_position
			#multipy by 3 to get into NA-space
			search_position*=3
			print "in na space no offset : ",search_position
			#add the offset for query start
			search_position+=(query_start-1)
			print "in na space WITH offset : ",search_position
			cdr_map[query_name]=int(search_position)
			print query_name+" ---> "+str(cdr_map[query_name])
			writer.write(str(query_name)+"\t"+str(cdr_map[query_name])+"\n")
		else:
			print "WARNING , regex "+str(myre.pattern)+" NOT FOUND at query = "+str(query_name)+" in file ",k," recording (-1) as CDR3 end!"
			writer.write(str(query_name)+"\t-1\n")
			
		#print "\n\n\n"
	writer.close()



def writeRegionsFromIGBLASTResult(k,o):
	infile=k
	reg_map=dict()
	regions=["FWR1","CDR1","FWR2","CDR2","FWR3","CDR3"]
	currentQuery=None
	currentMode=None
	reader=open(k,'r')
	for line in reader:
		line=line.strip()
		if(line.startswith("#")):
			if(line.startswith("# Alignment summary")):
				currentMode="summary"
			elif(line.startswith("# Query:")):
				## Query: IGHV5-a*04
				qre=re.compile(r'^#\s+Query:\s*(.*)$')
				rr=re.search(qre,line)
				if(rr):
					currentQuery=rr.group(1)
					currentQuery=currentQuery.strip()
			else:
				currentMode=None
		if( (not(currentQuery in reg_map)) and (not(currentQuery==None))     ):
			reg_map[currentQuery]=dict()
		if(currentMode=="summary"):
			for region in regions:
				if(not(region in reg_map[currentQuery])):
					reg_map[currentQuery][region]=dict()
					reg_map[currentQuery][region]["start"]=(-1)
					reg_map[currentQuery][region]["stop"]=(-1)
				if(line.startswith(region) or line.startswith(re.sub(r'W','',region))):
					pieces=line.split("\t")
					start=pieces[1]
					stop=pieces[2]
					if(isIntegral(start)):
						reg_map[currentQuery][region]["start"]=int(start)
					else:
						reg_map[currentQuery][region]["start"]=(-1)
					print "got start ",start," for q=",currentQuery," region=",region
					if(not(region=="CDR3")):
						if(isIntegral(stop)):
							reg_map[currentQuery][region]["stop"]=int(stop)
						else:
							reg_map[currentQuery][region]["stop"]=(-1)
						print "got stop ",stop," for q=",currentQuery," region=",region
					else:
						reg_map[currentQuery][region]["stop"]=(-1)
	key_list=reg_map.keys()
	key_list.sort()
	writer=open(o,'w')
	for i in range(len(key_list)):
		k=key_list[i]
		print "Now looking at ",str(k)
		printMap(reg_map[k])
		writer.write(str(k))
		for region in regions:
			writer.write("\t"+str(reg_map[k][region]["start"])+"\t"+str(reg_map[k][region]["stop"]))
		#if(i<len(key_list)-1):
		if(True):
			writer.write("\n")
	writer.close()

if (__name__=="__main__"):
	#organisms=["Mus_musculus","human"]
	#organisms=["Mus_musculus"]
	organisms=["human"]

	for org in organisms:
		igblast="/home/data/DATABASE/01_22_2014/"+org+"/ReferenceDirectorySet/KABAT_TR/igblastn.kabat.out"
		out="/home/data/DATABASE/01_22_2014/"+org+"/ReferenceDirectorySet/KABAT_TR/Vlookup.tsv"
		if(os.path.exists(igblast) and not os.path.exists(out)):
			print "proceed on kabat V"
			writeRegionsFromIGBLASTKabatResult(igblast,out)
		xml="/home/data/DATABASE/01_22_2014/"+org+"/ReferenceDirectorySet/KABAT_TR/blastx.out.xml"
		out="/home/data/DATABASE/01_22_2014/"+org+"/ReferenceDirectorySet/KABAT_TR/Jlookup.tsv"
		if(os.path.exists(xml) and not os.path.exists(out)):
			print "proceed on kabat J"
			writeKabatJCDR3End(xml,out)





#cmd for running igblastn to extract kabat info
#igblastn -query ../human_IG_V.fna  -germline_db_V /usr/local/igblast_from_lonestar/database/human_gl_V   -germline_db_D /usr/local/igblast_from_lonestar/database/human_gl_D  -germline_db_J /usr/local/igblast_from_lonestar/database/human_gl_J -auxiliary_data /usr/local/igblast_from_lonestar/optional_file/human_gl.aux    -domain_system kabat   -organism human  -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'

#blastx -query ../human_TR_J.fna -db /home/esalina2/Downloads/IMGT_02.13.2014/www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP  -outfmt 5  -max_target_seqs 1   > blastx.out.xml



