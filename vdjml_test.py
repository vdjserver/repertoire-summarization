#!/usr/bin/python

import vdjml
import re
from utils import trimList
from segment_utils import getFastaListOfDescs,getSegmentName
from igblast_utils import *



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



def makeRegionAlnMetricsFromValn(valn,qstart,regionmap):
	score=1
	identity=0
	insertions=0
	deletions=0
	substitutions=0
	aln_start_index=(int(regionmap['from']))-int(qstart)
	region_len=(int(regionmap['to']))-(int(regionmap['from']))
	counted=0
	temp_i=0
	print "Trying to get aln_start_index...."
	while(counted<int(regionmap['from'])-int(qstart)):
		if(not(valn[0][temp_i]=="-")):
			counted+=1
		temp_i+=1
	aln_start_index=temp_i
	print "qstart is ",qstart
	print "The whole aln : \n"
	print getNiceAlignment(valn)
	counted=0
	while(counted<region_len):
		if(not(valn[0][temp_i]=="-")):
			counted+=1
		temp_i+=1
	aln_end_index=temp_i+1
	counted=0
	print "The valn len =",len(valn[0])
	print "aln_start_index=",aln_start_index
	print "aln_end_index=",aln_end_index
	aln_region=list()
	for l in range(len(valn)):
		aln_region_temp=valn[l][aln_start_index:aln_end_index]
		aln_region.append(aln_region_temp)
	print "From a region "+regionmap['region']+" we have :\n"
	print getNiceAlignment(aln_region)
	print "\n\n"
	for i in range(len(aln_region[0])):
		qbase=aln_region[0][i]
		sbase=aln_region[2][i]
		#print "qbase=",qbase,"sbase=",sbase
		if(qbase==sbase):
			identity+=1
		elif(qbase=="-"):
			deletions+=1
		elif(sbase=="-"):
			insertions+=1
		else:
			substitutions+=1
	#aln_metric=vdjml.Match_metrics(score,100.0*float(identity)/float(region_len+1),insertions,deletions,substitutions)
	aln_metric =vdjml.Match_metrics(100.0*float(identity)/float(region_len+1),substitutions=substitutions,insertions=insertions,deletions=deletions)
	return aln_metric






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
	fact=vdjml.Result_factory(meta)
	fact.set_default_aligner('IGBLAST', '123', '45', 67)
	fact.set_default_gl_database(
				     'human_gl_S', 
				     '123-5', 
				     'blah', 
				     'http://www.blah.org'
				     )
	fact.set_default_num_system(vdjml.Num_system.imgt)
	rrw1 = vdjml.Result_writer('py_out_2.vdjml', meta)
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
			ref=re.search('^#\s+BLAST\s+processed\s+\d+\s+queries',igblast_line)
			if(rem or ref):
				if(not(current_query==None)):
					print "*******************************************\n"
					print "NEED TO SERIALIZE FOR QUERY=",current_query
					rb1 = fact.new_result(current_query)
					firstV=None
					firstD=None
					firstJ=None
					vaseq=None
					qvaseq=None
					daseq=None
					qdaseq=None
					jaseq=None
					qjaseq=None
					top_btop=dict()
					valn=None
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
						btop_to_add=kvmap['BTOP']
						#print "making interval with values : ",kvmap['q. start']," AND ",kvmap['q. end']
						query_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['q. start']),int(kvmap['q. end']))
						segment_type_to_add=kvmap['segment']
						hit_name_to_add=kvmap['subject ids']
						#hit_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['s. start']),int(kvmap['s. end']))
						if(int(kvmap['s. start'])<int(kvmap['s. end'])):
							hit_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['s. start']),int(kvmap['s. end']))
						else:
							#strand!
							hit_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['s. end']),int(kvmap['s. start']))							
						#metrics_to_add=vdjml.Match_metrics(int(kvmap['score']),float(kvmap['% identity']),int(9999),int(9999),int(kvmap['mismatches']))
						metrics_to_add=vdjml.Match_metrics(float(kvmap['% identity']),int(kvmap['score']),int(kvmap['mismatches']))
						smb1 = rb1.add_segment_match(
							btop_to_add,
							query_interval_to_add,
							segment_type_to_add,
							hit_name_to_add,
							hit_interval_to_add,
							metrics_to_add)
						if(firstV==None and segment_type_to_add=='V'):
							firstV=smb1
							vaseq=kvmap['subject seq']
							qvaseq=kvmap['query seq']
							top_btop['V']=kvmap['BTOP']
							alnDebugFlag=False
							if(current_query=="HZ8R54Q01C3N51"):
								alnDebugFlag=True
							#valn=buildAlignmentWholeSeqs(top_btop['V'],qvaseq,vaseq,alnDebugFlag)
							valn=buildAlignmentWholeSeqsDirect(qvaseq,vaseq)
							if(alnDebugFlag):
								print "The debug is :",
								print getNiceAlignment(valn)
								#sys.exit(0)
							valn_qstart=kvmap['q. start']
						if(firstD==None and segment_type_to_add=='D'):
							firstD=smb1
							daseq=kvmap['subject seq']
							qdaseq=kvmap['query seq']
							top_btop['D']=kvmap['BTOP']
						if(firstJ==None and segment_type_to_add=='J'):
							firstJ=smb1
							jaseq=kvmap['subject seq']
							qjaseq=kvmap['query seq']
							top_btop['J']=kvmap['BTOP']
					#getWholeChainAlignment(qvaseq,vaseq,qdaseq,daseq,qjaseq,jaseq,top_btop)
					scb=None
					#note, because igBLAST "anchors" on the V alignment, if V aligns, then D and J are both optional, but if V doesn't align, then no D alignment exists and no J alignment exists!
					if(not(firstV==None)):
						if(firstD==None):
							if(not(firstJ==None)):
								scb=rb1.add_segment_combination(firstV.get().id(),firstJ.get().id())
							else:
								#firstJ is none
								scb=rb1.add_segment_combination(firstV.get().id())
						else:
							#firstD not none
							if(not(firstJ==None)):
								scb=rb1.add_segment_combination(firstV.get().id(),firstD.get().id(),firstJ.get().id())						
							else:
								scb=rb1.add_segment_combination(firstV.get().id(),firstD.get().id())
					else:
						pass
					print "\n\n\tAlignment summary Info : "
					print "Alignment summary fields : ",summary_fields
					print "Alignment summary vals : "
					for a in range(len(summary_vals_list)):
						asMap=makeMap(summary_fields,summary_vals_list[a])
						if(not(asMap['region'].startswith("Total") or asMap['region'].startswith("CDR3"))):
							print "\tSummary item "+str(a)+" : "
							print "\t",summary_vals_list[a]
							print "\tThis is a particular alignment summary map : "
							printMap(asMap)	
							region_to_add=asMap['region']
							interval_to_add=vdjml.Interval.first_last_1(int(asMap['from']),int(asMap['to']))
							metrics_to_add=makeRegionAlnMetricsFromValn(valn,valn_qstart,asMap)
							scb.add_region(
								region_to_add,
								interval_to_add,
								metrics_to_add)
					rrw1(rb1.get())
					#rrw1.close()
					print "\n\n\n"
					if(len(vdjr_vals_list)>0):
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
					if(len(junction_vals_list)>0):
						print "\n\nJUNCTION summary:"
						print "\tJunction Fields list : ",junction_fields
						print "\tJunction vals list",junction_vals_list
						jmap=makeMap(junction_fields,junction_vals_list[0])
						print "\tJunction map : "
						printMap(jmap)
				if(rem):
					current_query=rem.group(1)
				else:
					current_query=None
				vdjr_vals_list=list()
				junction_vals_list=list()
				summary_vals_list=list()
				hit_vals_list=list()
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

    

scanOutputToVDJML("/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna.igblastn.imgt.out.first5.out","/dev/null","/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna",
	["/home/esalina2/round1_imgt/human_IG_V.fna","/home/esalina2/round1_imgt/human_IG_D.fna","/home/esalina2/round1_imgt/human_IG_J.fna"])

