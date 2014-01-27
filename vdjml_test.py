#!/usr/bin/python

import vdjml
import re
from utils import trimList
from segment_utils import getFastaListOfDescs,getSegmentName




def printMap(m):
	i=0
	for k in m:
		print "#"+str(i)+"\t",k,"->",m[k]
		i+=1


def getListOfDescsFromDBList(db_fasta_list):
	l=list()
	for f in range(len(db_fasta_list)):
		descs=getFastaListOfDescs(db_fasta_list[f])
		l.append(descs)
	return l



def makeMap(col_list,val_tab_str):
	val_list=val_tab_str.split("\t")
	kvmap=dict()
	#if(len(val_list)!=len(col_list)):
		#print "ERROR, COL_LIST LENGTH NOT SAME AS VAL LENGTH : "
		#for x in range(max(len(val_list),len(col_list))):
		#	if(x<len(col_list)):
		#		print col_list[x]+"\t->\t",
		#	else:
		#		print "EMPTY\t->\t"
		#	if(x<len(val_list)):
		#		print val_list[x]
		#	else:
		#		print "EMPTY"
		#print "col list:"
		#print col_list
		#print "val list:"
		#print val_list
	for c in range(len(col_list)):
		if(c<len(val_list)):
			kvmap[col_list[c]]=val_list[c]
		else:
			kvmap[col_list[c]]=None
	return kvmap


def scanOutputToVDJML(input_file,output_file,fasta_paths,db_fasta_list,debug=False):
	if(debug):
		print "debug is true"
	else:
		print "debug is false"
	INPUT=open(input_file,'r')
	line_num=0
	current_query=None
	current_query_id=(-1)
	vdjr_vals_list=list()
	junction_vals_list=list()
	summary_vals_list=list()
	hit_vals_list=list()
	ddl=getListOfDescsFromDBList(db_fasta_list)
	vlist=ddl[0]
	dlist=ddl[1]
	jlist=ddl[2]
	#print "VLIST :",vlist
	#print "DLIST :",dlist
	#print "JLIST :",jlist
	meta=None
	fact=None
	rb1=None
	meta = vdjml.Results_meta()
	rrw1 = vdjml.Result_writer('/tmp/py_out_2.vdjml', meta)
	for igblast_line in INPUT:
		line_num+=1
		igblast_line=igblast_line.strip()
		print igblast_line
		if(re.compile('^#').match(igblast_line)):
			rem=re.search('^#\s+IGBLASTN\s([^\s]+)\s*$',igblast_line)
			if rem:
				igblast_version=rem.group(1)
				mode="IGB_VERSION"
			rem=re.search('^#\s+Query:\s(.*)',igblast_line)
			if rem:
				if(not(current_query==None)):
					print "NEED TO SERIALIZE FOR QUERY=",current_query
					if current_query_id==(0):
						#make new factories						
						fact = vdjml.Result_factory(meta)
						fact.set_default_aligner('IGBLAST', '123', '45', 67)
						fact.set_default_gl_database(
									     'human_gl_S', 
									     '123-5', 
									     'blah', 
									     'http://www.blah.org'
									     )
						fact.set_default_num_system(vdjml.Num_system.imgt)
						pass
					rb1 = fact.new_result(current_query)
					rrw1(rb1.get())
					#rrw1.close()
					print "\n\n\n"
					print "rearrangment summary "
					print "\tFields : "
					print rearrangement_summary_fields
					for f in range(len(rearrangement_summary_fields)):
						print "\t\tfield "+str(f)+" : "+rearrangement_summary_fields[f]
					print "\tValues:"
					print vdjr_vals_list
					for v in range(len(vdjr_vals_list)):
						print "\t\tval : "+str(v)+" : "+vdjr_vals_list[v]
					kvmap=makeMap(rearrangement_summary_fields,vdjr_vals_list[0])
					for k in kvmap:
						print "\t"+k+"\t->\t"+kvmap[k]
					print "\n\nJUNCTION summary:"
					print "\tJunction Fields list : ",junction_fields
					print "\tJunction vals list",junction_vals_list
					jmap=makeMap(junction_fields,junction_vals_list[0])
					print "\tJunction map : "
					printMap(jmap)
					print "\n\n\tAlignment summary Info : "
					print "Alignment summary fields : ",summary_fields
					print "Alignment summary vals : "
					for a in range(len(summary_vals_list)):
						print "\tSummary item "+str(a)+" : "
						print "\t",summary_vals_list[a]
						asMap=makeMap(summary_fields,summary_vals_list[a])
						print "\tThis is a particular alignment summary map : "
						printMap(asMap)
						print "\n"
					for h in range(len(hit_vals_list)):
						kvmap=makeMap(hit_fields,hit_vals_list[h])
						print "\n\n\nSHOWING A ",kvmap['segment']," segment! (hit #"+str(h+1)+" of "+str(len(hit_vals_list))+")"
						printMap(kvmap)
						if(kvmap['segment']=="D"):
							lookup=dlist
						elif(kvmap['segment']=="V"):
							lookup=vlist
						else:
							lookup=jlist
						print "The nice name is ",getSegmentName(kvmap['subject acc.ver'],lookup)
					print "\n\n\n"
					break;
				current_query=rem.group(1)
				mode="query_retrieve"
				current_query_id+=1
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
					summary_fields.insert(0,'region')
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
			rem=re.search('#\sFields:',igblast_line)
			if rem:
				mode="hit_fields"
				temp_line=igblast_line
				headers=temp_line
				headers=re.compile('^#\s*Fields:\s*').sub('',headers,0)
				hit_fields=headers.split(',')
				hit_fields=trimList(hit_fields)
				hit_fields.insert(0,"segment")
				#for idx, val in enumerate(hit_fields):
				#	hit_fields[idx]=re.compile('[\s\.]+').sub('_',hit_fields[idx],0)
				#	hit_fields[idx]=re.compile('%').sub('pct',hit_fields[idx],0)
				#	hit_fields[idx]=re.compile('[^a-zA-Z]+$').sub('',hit_fields[idx],0)
				#	hit_fields[idx]=re.compile('^[^a-zA-Z]+').sub('',hit_fields[idx],0)
				#	hit_fields[idx]=re.compile('\/').sub('_',hit_fields[idx],0)
			rem=re.search('^#\s(\d+)\shits\sfound\s*',igblast_line)
			if rem:
				mode="hit_table"
				current_num_hits=rem.group(1)
		elif (not (re.compile('^\s*$').match(igblast_line))):
			#print "mode=",mode
			if(mode=="vdj_rearrangement"):
				vdjr_vals_list.append(igblast_line)
			elif(mode=="vdj_junction"):
				junction_vals_list.append(igblast_line)
			elif(mode=="vdj_alignment_summary"):
				summary_vals_list.append(igblast_line)
			elif(mode=="hit_table"):
				#print "in mode=hit_table",
				#print "hit_fields is ",hit_fields," with length=",len(hit_fields)
				#print "hit_vals after split",hit_vals
				hit_vals_list.append(igblast_line)
			else:
				print "unhandled mode=",mode

    

scanOutputToVDJML("/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna.igblastn.imgt.out","/dev/null","/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna",
	["/home/esalina2/round1_imgt/human_IG_V.fna","/home/esalina2/round1_imgt/human_IG_D.fna","/home/esalina2/round1_imgt/human_IG_J.fna"])

