#!/usr/bin/env python


from os import listdir
from os.path import isfile, join
import os.path
import io
import subprocess
import re
from collections import defaultdict
import pprint
import json
import sqlite3
import random
import glob
import ntpath
import sys
import time










def makeRegionLookupTableFile(infile,outfile):
	INPUT=open(infile,'r')
	OUTPUT=open(outfile,'w')
	query_num=(-1)
	#CDR1	3	10	8	7	1	0	87.5
	#FWR2	11	61	51	50	1	0	98
	#CDR2	62	82	21	21	0	0	100
	#FWR3	83	196	114	113	1	0	99.1
	#CDR3 (V region only)	197	204	8	6	2	0	75
	regions=["FWR1","CDR1","FWR2","CDR2","FWR3","CDR3"]
	q_name=""
	region_dict=dict()
	for line in INPUT:
		line=line.strip()
		q_search=re.search('^#\s+Query:\s(.*)',line)
		if(q_search):
			query_num+=1
			if(query_num==0):
				OUTPUT.write("Query\t")
				for idx, val in enumerate(regions):
					OUTPUT.write(val+"_S\t"+val+"_E")
					if(idx<(len(regions)-1)):
						OUTPUT.write("\t")
						#add \t between items if not at the end				
					else:
						OUTPUT.write("\n")
					#write the data
				#write the headers only for first query!
			if(len(region_dict)>0):
				for idx, val in enumerate(regions):
					#OUTPUT.write("\t'"+str(val)+"'=>'"+str(region_dict[val])+"'\t\n")
					data=str(region_dict[val])
					pieces=data.split('\t')
					if(idx==0):
						OUTPUT.write(q_name+"\t")
					OUTPUT.write(pieces[0]+"\t"+pieces[1])
					if(idx<(len(regions)-1)):
						OUTPUT.write("\t")
					else:
						OUTPUT.write("\n")					
				#once the regions are written, flush them out to get new data
				for idx, val in enumerate(regions):
					region_dict[val]="(-1)\t(-1)"
			q_name=q_search.group(1)
		region_sep=str("|")
		region_regex=region_sep.join(regions)
		region_regex="^("+region_regex+")\s"
		region_search=re.search(region_regex,line)
		if(region_search):
			#print "got region_match line ",line
			#OUTPUT.write(line+"\n")
			#(from, to, length, matches, mismatches, gaps, percent identity)
			pieces=line.split('\t')
			pieces=trimList(pieces)
			if(line.startswith("CDR3")):
				#CDR3 logic
				region_dict["CDR3"]=str(pieces[1])+"\t"+str(pieces[2])
			else:
				region_dict[pieces[0]]=str(pieces[1])+"\t"+str(pieces[2])
	if(len(region_dict)>0):
		for idx, val in enumerate(regions):
			#OUTPUT.write("\t'"+str(val)+"'=>'"+str(region_dict[val])+"'\t\n")
			data=str(region_dict[val])
			pieces=data.split('\t')
			if(idx==0):
				OUTPUT.write(q_name+"\t")
			OUTPUT.write(pieces[0]+"\t"+pieces[1])
			if(idx<(len(regions)-1)):
				OUTPUT.write("\t")
			else:
				OUTPUT.write("\n")
	OUTPUT.close()
	INPUT.close()




def addMetadataTableToSQLITEDatabase(igblast_version,databases,domain_class,query_path):
	creation_sql="CREATE TABLE IF NOT EXISTS metadata (igblast_version,databases,domain_class,query_path);"
	insertion_sql="INSERT INTO metadata (igblast_version,databases,domain_class,query_path) VALUES ('"+igblast_version+"','"+databases+"','"+domain_class+"','"+query_path+"');"
	return creation_sql+"\n"+insertion_sql


def scanOutputToSQLite(input_file,output_file,fastaPath,debug=False):
	if(debug):
		print "debug is true"
	else:
		print "debug is false"
	#########################
	# variables
	igblast_version=""
	igblast_databases=""
	dom_class=""
	created_vdj_rearrangement=0
	created_junction_table=0
	created_summary_table=0
	created_hit_table=0
	current_query_id=(-1)
	current_hit_id=0
	mode="prologue"
	INPUT=open(input_file,'r')
	sqlite_path=os.environ['VDJ_SERVER_SQLITE3']
	if(sqlite_path==None):
		print "ERROR, FAILURE TO FIND VDJ_SERVER_SQLITE3_PATH!"
		return
	#p = subprocess.Popen([sqlite_path, output_file], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
	p = subprocess.Popen([sqlite_path, output_file], stdin=subprocess.PIPE)
	sqlite_stdin=p.stdin
	sqlite_stdout=p.stdout
	sqlite_stdin.write("begin;\n")
	for igblast_line in INPUT:
		igblast_line=igblast_line.strip()
		#print igblast_line
		if(current_query_id>0 and (current_query_id%10000)==0):
			sqlite_stdin.write("end;\n")
			sqlite_stdin.write("begin;\n")
		if(re.compile('^#').match(igblast_line)):
			rem=re.search('^#\s+IGBLASTN\s([^\s]+)\s*$',igblast_line)
			if rem:
				igblast_version=rem.group(1)
				mode="IGB_VERSION"
			rem=re.search('^#\s+Query:\s(.*)',igblast_line)
			if rem:
				current_query=rem.group(1)
				mode="query_retrieve"
				current_query_id+=1
				if(current_query_id<=0):
					sqlite_stdin.write("CREATE TABLE query_table(input_query_id INT,input_query_name);\n")
					sqlite_stdin.write("CREATE UNIQUE INDEX IF NOT EXISTS seq_id_index ON query_table (input_query_id);\n")
				query_insert=str("INSERT INTO query_table (input_query_id,input_query_name) VALUES ('"+str(current_query_id)+"',\""+str(current_query)+"\");\n")
				sqlite_stdin.write(query_insert)
			rem=re.search('^#\s+Database:\s(.*)$',igblast_line)
			if rem:
				igblast_databases=rem.group(1)
				mode="dbparse"
			rem=re.search('^#\sDomain\sclassification\srequested:\s(.*)$',igblast_line)
			if rem:
				dom_class=rem.group(1)
				mode="dom_class"
			rem=re.search('^#\s(Note.*)$',igblast_line)
			if rem:
				note_text=rem.group(1)
			rem=re.search('^#\sV\(D\)J\sjunction',igblast_line)
			if rem:
				mode="vdj_junction"
				srem=re.search('\((\s*V\send,[^\)]+)\)',igblast_line)    
				if srem:
					headers=srem.group(1)
					headers=re.compile('\s*,\s*').sub(',',headers,0)
					headers=re.compile('\ ').sub('_',headers,0)
					headers=re.compile('\-').sub('_',headers,0)
					junction_fields=headers.split(',')
					junction_fields=trimList(junction_fields)
					junction_fields.insert(0,"current_query_id")
			rem=re.search('^#\sBLAST\sprocessed\s(\d+)\squer/',igblast_line)
			if rem:
				reported_queries_processed=rem.group(1)
			rem=re.search('^#\sAlignment\ssummary',igblast_line)
			if rem:
				mode="vdj_alignment_summary"
				srem=re.search('\(([^\)]+)\)',igblast_line)
				if srem:
					headers=srem.group(1)
					headers=re.compile('\s*,\s*').sub(',',headers,0)
					headers=re.compile('\ ').sub('_',headers,0)
					headers=re.compile('\-').sub('_',headers,0)
					summary_fields=headers.split(',')
					summary_fields.insert(0,"region")
					summary_fields.insert(0,"input_query_id")
			rem=re.search('^#\sHit\stable',igblast_line)
			if rem:
				mode="hit_table_start"
			rem=re.search('^#\sV\(D\)J\srearrangement',igblast_line)
			if rem:
				mode="vdj_rearrangement"
				srem=re.search('\(([^\)]+)\):\s*$',igblast_line)
				if srem:
					headers=srem.group(1)
					#print "headers is ",headers
					#print "Now to replace commas with surrounding whitespace...."
					headers=re.sub(r'\s*,\s*',',',headers)
					#print "headers is ",headers
					headers=re.sub(' ','_',headers,0)
					#print "headers is ",headers
					headers=re.sub('\-','_',headers,0)
					#print "headers is ",headers
					rearrangment_summary_fields=headers.split(',')
					#print "rearrangment_summary_fields is ",rearrangment_summary_fields
					rearrangement_summary_fields=trimList(rearrangment_summary_fields)
					#print "rearrangment_summary_fields is ",rearrangment_summary_fields
					rearrangment_summary_fields.insert(0,"input_query_id")
					#print "rearrangment_summary_fields is ",rearrangment_summary_fields
			rem=re.search('#\sFields:',igblast_line)
			if rem:
				mode="hit_fields"
				temp_line=igblast_line
				headers=temp_line
				headers=re.compile('^#\s*Fields:\s*').sub('',headers,0)
				hit_fields=headers.split(',')
				hit_fields=trimList(hit_fields)
				hit_fields.insert(0,"segment")
				hit_fields.insert(0,"hit_id")
				hit_fields.insert(0,"input_query_id")
				for idx, val in enumerate(hit_fields):
					hit_fields[idx]=re.compile('[\s\.]+').sub('_',hit_fields[idx],0)
					hit_fields[idx]=re.compile('%').sub('pct',hit_fields[idx],0)
					hit_fields[idx]=re.compile('[^a-zA-Z]+$').sub('',hit_fields[idx],0)
					hit_fields[idx]=re.compile('^[^a-zA-Z]+').sub('',hit_fields[idx],0)
					hit_fields[idx]=re.compile('\/').sub('_',hit_fields[idx],0)
			rem=re.search('^#\s(\d+)\shits\sfound\s*',igblast_line)
			if rem:
				mode="hit_table"
				current_num_hits=rem.group(1)
		elif (not (re.compile('^\s*$').match(igblast_line))):
			#print "mode=",mode
			if(mode=="vdj_rearrangement"):
				#sdfasdf
				vdjr_vals=igblast_line.split('\t')
				vdjr_vals=trimList(vdjr_vals)
				vdjr_vals.insert(0,str(current_query_id))
				field_sep="','"
				joined_fields=field_sep.join(rearrangment_summary_fields)
				if(created_vdj_rearrangement==0):
					vdjr_sql="CREATE TABLE vdj_rearrangement"
					#print "fields is ",rearrangment_summary_fields
					#print "joined fields : ",joined_fields
					rearrangment_summary_fields_sql="('"+joined_fields+"');"
					vdjr_sql=vdjr_sql+rearrangment_summary_fields_sql
					sqlite_stdin.write(vdjr_sql)
					created_vdj_rearrangement=1
				sql_insert_rearrangment="INSERT INTO vdj_rearrangement "
				val_sep="','"
				joined_vals=val_sep.join(vdjr_vals)
				sql_insert_rearrangment=sql_insert_rearrangment+"('"+joined_fields+"') VALUES ('"+joined_vals+"');"
				#print "The sql for insertion is ",sql_insert_rearrangment
				sqlite_stdin.write(sql_insert_rearrangment)
			elif(mode=="vdj_junction"):
				junction_vals=igblast_line.split('\t')
				junction_vals=trimList(junction_vals)
				junction_vals.insert(0,str(current_query_id))
				vdjj_fields_sep="','"
				vdjj_fields=vdjj_fields_sep.join(junction_fields)
				if(created_junction_table==0):
					vdjj_sql="CREATE TABLE vdj_junction("
					vdjj_fields_sep="','"
					vdjj_fields=vdjj_fields_sep.join(junction_fields)
					vdjj_sql+="'"+vdjj_fields+"','V_J_junction');"
					#print "JUNC TABLE CREATION : ",vdjj_sql
					sqlite_stdin.write(vdjj_sql)
					print "THE CREATE JUNCTION WAS ",vdjj_sql
					created_junction_table=1
				j_v_map=dict()
				j_v_map['V_J_junction']='N/A'
				#print "field len :",len(junction_fields),junction_fields
				#print "val len: ",len(junction_vals),junction_vals
				for V in range(len(junction_fields)):
					#print "V is ",V
					#print "KEY is ",junction_fields[V]
					#print "VAL is ",junction_vals[V]
					j_v_map[junction_fields[V]]=junction_vals[V]
				junction_fields=list()
				junction_vals=list()
				for key in j_v_map:
					junction_fields.insert(0,key)
					junction_vals.insert(0,j_v_map[key])
				junction_vals_sql="INSERT INTO vdj_junction ('"+vdjj_fields_sep.join(junction_fields)+"') VALUES ('"+vdjj_fields_sep.join(junction_vals)+"');"
				#print "JVS IS ",junction_vals_sql
				sqlite_stdin.write(junction_vals_sql)
			elif(mode=="vdj_alignment_summary"):
				#print "doing alignment summary"
				#print "summary fields " , summary_fields
				summary_vals=igblast_line.split('\t')
				summary_vals=trimList(summary_vals)
				fields_sep="','"
				fields_sql="('"+fields_sep.join(summary_fields)+"')"
				if(created_summary_table==0):
					create_summary_sql="CREATE TABLE alignment_summary "+fields_sql+";"
					#print "create sum sql = ",create_summary_sql
					sqlite_stdin.write(create_summary_sql)
					created_summary_table=1					
				summary_vals.insert(0,str(current_query_id))
				summary_insert_sql="INSERT INTO alignment_summary "+fields_sql+" VALUES ('"+fields_sep.join(summary_vals)+"');"
				#print "insert summary sql "+summary_insert_sql
				sqlite_stdin.write(summary_insert_sql)
			elif(mode=="hit_table"):
				#print "in mode=hit_table",
				#print "hit_fields is ",hit_fields," with length=",len(hit_fields)
				fields_sep="','"
				#print "hit_vals before split",igblast_line
				hit_vals=igblast_line.split('\t')
				#print "hit_vals after split",hit_vals
				hit_vals.insert(0,str(current_hit_id))
				hit_vals.insert(0,str(current_query_id))
				hit_vals=trimList(hit_vals)
				if(created_hit_table==0):
					#print "hit_fields length : ",len(hit_fields)
					hit_fields_sql=fields_sep.join(hit_fields)
					hit_table_sql_create="CREATE TABLE hit_table ('"+hit_fields_sql+"');"
					#print "hit_table creation : ",hit_table_sql_create
					sqlite_stdin.write(hit_table_sql_create)
					created_hit_table=1
				#print "need insert data ",igblast_line
				insert_hit_sql="INSERT INTO hit_table ('"+fields_sep.join(hit_fields)+"') VALUES ('"+fields_sep.join(hit_vals)+"');"
				#print insert_hit_sql
				#print "hit_fields len : "+str(len(hit_fields))+", hit_vals len : "+str(len(hit_vals))
				#print "hit_fields=",hit_fields
				#print "hit_vals=",hit_vals
				sqlite_stdin.write(insert_hit_sql)
				current_hit_id+=1
			else:
				print "unhandled mode=",mode
	metadata_sql=addMetadataTableToSQLITEDatabase(igblast_version,igblast_databases,dom_class,fastaPath)
	sqlite_stdin.write(metadata_sql)
	sqlite_stdin.write("end;\n")
	p.communicate('')






def getNiceAlignment(aln):
	a=str("")
	for idx, val in enumerate(aln):
		a+=val+"\n"
	return a

def repeatString(s,n):
	if(n<=0):
		return ""
	else:
		r=str("")
		for i in range(n):
			r+=s
		return r
		


def getWholeChainAlignment(qvseq,vseq,qdseq,dseq,qjseq,jseq,btop_map):
	valn=buildAlignmentWholeSeqs(btop_map['V'],qvseq,vseq)
	daln=None
	if(dseq is not None):
		daln=buildAlignmentWholeSeqs(btop_map['D'],qdseq,dseq)
	jaln=buildAlignmentWholeSeqs(btop_map['J'],qjseq,jseq)
	print "V aln :\n",getNiceAlignment(valn),"\n\n"
	if(daln is not None):
		print "D aln :\n",getNiceAlignment(daln),"\n\n"
	print "J aln :\n",getNiceAlignment(jaln),"\n\n"
	




def buildAlignmentWholeSeqsDirect(q,s):
	aln=["","",""]
	aln[0]=q
	aln[2]=s
	for b in range(len(aln[0])):
		qbase=aln[0][b]
		sbase=aln[0][b]
		if(qbase=="-" or sbase=="-"):
			aln[1]+=" "
		elif(not(qbase==sbase)):
			aln[1]+="X"
		else:
			aln[1]+="|"
	return aln		



def printNiceAlignment(q,s):
	joiner="\n"
	return joiner.join(buildAlignmentWholeSeqsDirect(q,s))




def buildAlignmentWholeSeqs(btop,q,s,debug=False,level=0):
	if(debug):	
		print repeatString("\t",level)+"At begging of call, btop=",btop,"q=",q,"s=",s
	aln=["","",""]
	aln[0]=str("")
	aln[1]=str("")
	aln[2]=str("")
	if(debug):
		print "now performing digital tests..."
	dm=re.search('^(\d+)[^0-9]',btop)
	adm=re.search('^(\d+)$',btop)
	if adm:
		if(debug):		
			print repeatString("\t",level)+"matched all digital..."
		btopv=int(adm.group(1))
		aln[0]=q
		aln[1]=repeatString("|",btopv)
		aln[2]=s
		if(debug):
			print repeatString("\t",level)+"from all digital btop=",btop," returning : "
			for x in range(len(aln)):
				print repeatString("\t",level)+aln[x]
		return aln
	elif dm:
		if(debug):
			print repeatString("\t",level)+"matched start digital..."
		digitString=dm.group(1)
		size=int(digitString)
		aln[0]=q[:size]
		aln[2]=s[:size]
		aln[1]=repeatString("|",size)
		if(debug):
			print repeatString("\t",level)+"From little btop :",digitString," got "
			for x in range(len(aln)):
				print repeatString("\t",level)+aln[x]
		rec=buildAlignmentWholeSeqs(btop[len(digitString):],q[(size-0):],s[(size-0):],debug,level+1)
		aln[0]+=rec[0]
		aln[1]+=rec[1]
		aln[2]+=rec[2]
		return aln
	if(debug):
		print "no digital tests passed...doing letter tests...."
	firstTwoLetters=re.search('^([a-z\\-])([a-z\\-])',btop,re.IGNORECASE)
	if firstTwoLetters:
		if(debug):
			print "aligns first two letters..."
		firstLetter=firstTwoLetters.group(1)
		secondLetter=firstTwoLetters.group(2)
		aln[0]=firstLetter
		aln[1]=str("X")
		aln[2]=secondLetter
		if(debug):		
			print "First-two letters alignment = \n"+getNiceAlignment(aln)
		if(len(btop)>2):
			rec=["","",""]
			if(firstLetter=="-" or secondLetter=="-"):
				if(firstLetter=="-"):
					rec=buildAlignmentWholeSeqs(btop[2:],q,s[1:],debug,level+1)
				elif(secondLetter=="-"):
					rec=buildAlignmentWholeSeqs(btop[2:],q[1:],s,debug,level+1)
			else:
				 rec=buildAlignmentWholeSeqs(btop[2:],q[1:],s[1:],debug,level+1)
			aln[0]+=rec[0]
			aln[1]+=rec[1]
			aln[2]+=rec[2]
		else:
			if(debug):
				print "returning next-to-last aln"
			return aln
	if(debug):
		print "returning last aln"
	return aln





def get_rearrangement_segment_counts_from_db(dbfile):
	conn=sqlite3.connect(dbfile)
	c = conn.cursor()
	counts_map=dict()
	get_rearrangements_sql="SELECT Top_V_gene_match,Top_D_gene_match,Top_J_gene_match FROM vdj_rearrangement"
	c.execute(get_rearrangements_sql)
	for row in c:
		#print row
		#print "top v: ",str(row[0])
		for segment in row:
			if segment in counts_map:
				counts_map[segment]=counts_map[segment]+1
			else:
				counts_map[segment]=1
	conn.close()
	return counts_map



def getIGBlastExecutablePath():
	return "/usr/local/bin/igblastn"


def getDomainClasses():
	dcs=["imgt","kabat"]
	return dcs



def printVDJBestAssignments(i):
	reader=open(i,'r')
	currentQuery=None
	qre=re.compile(r'^#\s+Query:\s+([^\s]+)\s*$')
	captureFlag=False
	for line in reader:
		line=line.strip()
		qsr=re.search(qre,line)
		if(qsr):
			currentQuery=qsr.group(1)
		elif(line.startswith("# V(D)J rearrangement summary")):
			captureFlag=True
		elif(captureFlag):
			pieces=line.split('\t')
			for p in range(len(pieces)):
				if ',' in pieces[p]:
					#the piece is a string of comma-separated values
					comma_pieces=pieces[p].split(',')
					pieces[p]=comma_pieces[0]
			joiner=','
			line=joiner.join(pieces)			
			#line=re.sub(r'\t',',',line)
			print currentQuery+","+line
			captureFlag=False
		
		





def test():
	pass
	#scanOutputToSQLite("/home/mlevin/project/stanford_davis/run3/out4","/home/esalina2/FLORIAN/merged.igblast.db","/home/mlevin/project/stanford_davis/run3/out4/merged.fasta")
	printVDJBestAssignments("/home/mlevin/project/stanford_davis/run3/out4/merged.iglbast.out")
	#print "Running test...."
	#counts_map=
	#in_hier_dir="/home/data/vdj_server/pipeline/vdj_ann/17_way/"
	#out_hier_dir="/tmp/tmp_json"
	#count_map=dict()
	#count_map['TRBD2*01']=10
	#count_map['TRBD2*02']=30	
	#count_map['TRBD1*01']=25
	#count_map['TRBD']=1
	#name="TRBD"
	#tree_dat_file="/home/data/vdj_server/pipeline/vdj_ann/17_way/TRBD.dat"
	#tree=getHierarchyTreeFromFile(tree_dat_file)
	#print "The total is ",get_total_tree(tree,name,count_map)
	#hierarchy_jsonify_batch(in_hier_dir,out_hier_dir,count_map)
	#basic_tree_test()
	#root=getRootFromFile("/home/data/vdj_server/pipeline/vdj_ann/17_way/IGHM.dat")
	#print "the root is ",root
	#dbfile="/home/esalina2/round1/all_data.processed.r0.small.fna.imbg.igblastn.7.out.db"	
	#counts_map=get_rearrangement_segment_counts_from_db(dbfile)
	#print counts_map
	#test_tree()
	#igblast_file="/home/esalina2/round1/all_data.processed.r0.small.fna.imbg.igblastn.7.out"
	#output_db_file="/home/esalina2/round1/all_data.processed.r0.small.fna.imbg.igblastn.7.out.db"
	#fasta_file="/home/esalina2/round1/all_data.processed.r0.small.fna"
	#scanOutputToSQLite(igblast_file,output_db_file,fasta_file)
	#infile="/home/esalina2/region_annotation/human_gl_V.fna.imgt.7"
	#outfile="/dev/stdout"
	#makeRegionLookupTableFile(infile,outfile)
	#testDBFile="/tmp/delme.db"
	#if os.path.isfile(testDBFile):
	#	os.remove(testDBFile)
	#scanOutputToSQLite("/home/esalina/repserver/igblast_routines/igblast_2.out",testDBFile,"/tmp/fasta/fna")
	#	q="ATTACA"
	#	s="GATTACA"
	#btop="-G6"
	#aln=buildAlignmentWholeSeqs(btop,q,s)
	#print getNiceAlignment(aln)
	#v_name="IGHV4-34*01"
	#d_name="IGHD3-3*01"
	#j_name="IGHJ4*02"
	#q_name="AY671579.1.partial"
	#q_seq="GTCAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGATCCGCCAGCCCCCAGGGCAAGGGGCTGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGGCACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCTGTGTATTACTGTGCGAGAGGCAGTACCGGCCGATTTTTGGAGTGGTTATTATACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGGGAGTCGATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTGGCCGTTGGCTGCCTCGCACAGGACTTCCTTCCCGACTCCATCACTTTCTCCTGGAAATACAAGAACAACTCTGACATCAGCAGCACCCGGGGCTTCCCATCAGTCCTGAGA"
	#v_seq="CAGGTGCAGCTACAGCAGTGGGGCGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCGCTGTCTATGGTGGGTCCTTCAGTGGTTACTACTGGAGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCAATCATAGTGGAAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCTGTGTATTACTGTGCGAGAGG"
	#v_start=int(1)
	#v_stop=int(293)
	#q_v_start=int(3)
	#q_v_end=int(295)
	#q_v_btop=str("126CA1AG3GCCTTG79GA78")
	#q_btop=q_seq[(q_v_start-1):(q_v_end)]
	#v_btop=v_seq[(v_start-1):(v_stop)]
	#print "q_btop is ",q_btop
	#print "v_btop is ",v_btop
	#q_v_aln=buildAlignmentWholeSeqs(q_v_btop,q_btop,v_btop)
	#print "the qv nice is "
	#print "\n"+getNiceAlignment(q_v_aln)
	#d_seq="gtattacgatttttggagtggttattatacc"
	#j_seq="actactttgactactggggccagggaaccctggtcaccgtctcctcag"
	#q_d_start=int(306)
	#q_d_end=int(329)
	#d_start=int(7)
	#d_end=int(30)
	#q_d_btop="24"
	#q_btop=q_seq[(q_d_start-1):(q_d_end)]
	#d_btop=d_seq[(d_start-1):(d_end)]
	#q_d_aln=buildAlignmentWholeSeqs(q_d_btop,q_btop,d_btop)
	#print "len q_btop is ",len(q_btop)
	#print "len d_btop is ",len(d_btop)
	#print "the qd nice is "
	#print "\n"+getNiceAlignment(q_d_aln)
	#q_j_btop="46"
	#j_start=int(3)
	#j_end=int(48)
	#q_j_start=int(327)
	#q_j_end=int(372)
	#q_btop=q_seq[(q_j_start-1):(q_j_end)]
	#j_btop=j_seq[(j_start-1):j_end]
	#q_j_aln=buildAlignmentWholeSeqs(q_j_btop,q_btop,j_btop)
	#print "\n"+getNiceAlignment(q_j_aln)
	#v_d_junc=q_seq[(q_v_end):(q_d_start-1)]
	#print "The extracted VD junction is ",v_d_junc

if (__name__=="__main__"):
	import sys
	#print "call test"
	test()
