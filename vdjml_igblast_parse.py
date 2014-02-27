#!/usr/bin/env python

import vdjml
import re
#from segment_utils import getFastaListOfDescs,getSegmentName,IncrementMapWrapper,looksLikeAlleleStr
import argparse
from Bio import SeqIO
#from cdr3_hist import CDR3LengthAnalysisVDMLOBJ
#from igblast_utils import getDomainClasses
#from segment_utils import getRegionAlignmentFromLargerVAlignment,getRegionSpecifcCharacterization,getCDR3RegionSpecificCharacterization,getVRegionStartAndStopGivenRefData,getADJCDR3EndFromJAllele,getTheFrameForThisReferenceAtThisPosition

#parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
#                   help='an integer for the accumulator')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                   help='sum the integers (default: find the max)')
#
#args = parser.parse_args()
#print args.accumulate(args.integers)




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
	#aln_metric=vdjml.Match_metrics(
	#		score,100.0*float(identity)/float(region_len+1),
	#		insertions,
	#		deletions,
	#		substitutions
	#)
	aln_metric =vdjml.Match_metrics(100.0*float(identity)/float(region_len+1),substitutions=substitutions,insertions=insertions,deletions=deletions)
	return aln_metric



#given a characterization map for a region, make a metrics object
def makeMetricFromCharMap(charMap):
	pct_id=float(charMap['pct_id'])
	substitutions=charMap['base_sub']
	insertions=charMap['insertions']
	deletions=charMap['deletions']
	aln_metric=vdjml.Match_metrics(
		pct_id,
		substitutions=substitutions,
		insertions=insertions,
		deletions=deletions
		)
	return aln_metric


def vdjml_read_serialize(
		current_query,					#current query name
		fact, 						#handle to VDJML factory
		hit_vals_list,					#hit values (list of lists)
		hit_fields,					#list of hit fields
		summary_fields,					#summary fields
		summary_vals_list,				#summary values list (list of lists)
		rearrangement_summary_fields,vdjr_vals_list,	#rearrangment summary (fields (list) and values (list of lists))
		junction_fields,junction_vals_list		#junction (fields (list) and values (list of lists))
		):
	#print "*******************************************\n"
	#print "NEED TO SERIALIZE FOR QUERY=",current_query
	rb1 = fact.new_result(current_query)
	#print "the base is ",imgtdb_obj.getDirBase()
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
	firstVMap=None
	firstDMap=None
	firstJMap=None
	firstVLine=None
	firstDLine=None
	firstJLine=None
	##################################################################################
	#if there are no hits, return an empty object..... that contains just the read name! :(
	if(len(hit_vals_list)==0):
		return rb1
	################################################################
	#this is where the hits are picked up and put in VDJML objects
	for h in range(len(hit_vals_list)):
		kvmap=makeMap(hit_fields,hit_vals_list[h])
		btop_to_add=kvmap['BTOP']
		query_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['q. start']),int(kvmap['q. end']))
		if(int(kvmap['q. start'])>int(kvmap['q. end'])):
			query_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['q. end']),int(kvmap['q. start']))
		segment_type_to_add=kvmap['segment']
		hit_name_to_add=kvmap['subject ids']
		read_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['q. start']),int(kvmap['q. end']))
		if(int(kvmap['s. start'])<int(kvmap['s. end'])):
			hit_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['s. start']),int(kvmap['s. end']))
		else:
			#strand!
			hit_interval_to_add=vdjml.Interval.first_last_1(int(kvmap['s. end']),int(kvmap['s. start']))							
		#printMap(kvmap)
		metrics_to_add=vdjml.Match_metrics(
			identity=float(kvmap['% identity']),
			score=int(kvmap['score']),
			)
		#this is where the hit is added to the general list
		smb1 = rb1.add_segment_match(
			btop=btop_to_add,
			read_range=read_interval_to_add,
			seg_name=hit_name_to_add,
			vdj=segment_type_to_add,
			gl_range=hit_interval_to_add,
			metric=metrics_to_add
			)
		#################################################################
		#this is where the top hits are retrieved for the combination/rearrangement
		if(firstV==None and segment_type_to_add=='V'):
			firstV=smb1
			firstVMap=kvmap
			vaseq=kvmap['subject seq']
			qvaseq=kvmap['query seq']
			top_btop['V']=kvmap['BTOP']
			alnDebugFlag=False
		elif(firstD==None and segment_type_to_add=='D'):
			firstD=smb1
			firstDMap=kvmap
			daseq=kvmap['subject seq']
			qdaseq=kvmap['query seq']
			top_btop['D']=kvmap['BTOP']
		elif(firstJ==None and segment_type_to_add=='J'):
			firstJMap=kvmap
			firstJ=smb1
			jaseq=kvmap['subject seq']
			qjaseq=kvmap['query seq']
			top_btop['J']=kvmap['BTOP']
	scb=None

	####################################################################################
	#this is where the combination is added based on the firstV,firstD,firstJ from above
	if(not(firstV==None)):
		if(firstD==None):
			#no D, so maybe a light chain?
			if(not(firstJ==None)):
				scb=rb1.add_segment_combination(firstV.get().id(),firstJ.get().id())#,productive=productive_flag)
			else:
				#firstJ is none
				scb=rb1.add_segment_combination(firstV.get().id())#,productive=False)
		else:
			#firstD not none
			if(not(firstJ==None)):
				scb=rb1.add_segment_combination(firstV.get().id(),firstD.get().id(),firstJ.get().id())#,productive=productive_flag)						
			else:
				scb=rb1.add_segment_combination(firstV.get().id(),firstD.get().id())#,productive=False)
	else:
		#firstV is none
		pass

	

	##############################################################
	# rearrangment summary
	# TOPV, TOPD, TOPJ
	#for f in range(len(rearrangement_summary_fields)):
	#	#print "\t\tfield "+str(f)+" : "+rearrangement_summary_fields[f]
	#	pass

	###############################################################
	#alignment summary (regions)
	for r in range(len(summary_vals_list)):
		kvmap=makeMap(summary_fields,summary_vals_list[r])
		region=kvmap['region'].strip()
		if(not(region.startswith("Total")) and not(region.startswith("CDR3"))):
			reg_from=kvmap['from']
			reg_to=kvmap['to']
			gaps=kvmap['gaps']
			length=kvmap['length']
			matches=kvmap['matches']
			pct_id=kvmap['percent_identity']
			try:
				reg_name=region
				reg_interval=vdjml.Interval.first_last_1(int(reg_from),int(reg_to))
				mm=vdjml.Match_metrics(identity=float(pct_id))
				#print "adding a region"
				scb.add_region(name=reg_name,read_range=reg_interval,metric=mm)
			except:
				#print "got an exception!"
				pass
	if(len(vdjr_vals_list)>0):
		#print "rearrangment summary "
		#print "\tFields : "
		#print rearrangement_summary_fields
		for f in range(len(rearrangement_summary_fields)):
			#print "\t\tfield "+str(f)+" : "+rearrangement_summary_fields[f]
			pass
		#print "\tValues:"
		#print vdjr_vals_list
		for v in range(len(vdjr_vals_list)):
			#print "\t\tval : "+str(v)+" : "+vdjr_vals_list[v]
			pass
		kvmap=makeMap(rearrangement_summary_fields,vdjr_vals_list[0])
		for k in kvmap:
			#print "\t"+k+"\t->\t"+kvmap[k]
			pass
		top_segs=["V","D","J"]
	if(len(junction_vals_list)>0):
		#print "\n\nJUNCTION summary:"
		#print "\tJunction Fields list : ",junction_fields
		#print "\tJunction vals list",junction_vals_list
		jmap=makeMap(junction_fields,junction_vals_list[0])
		#print "\tJunction map : "
		#printMap(jmap)
	return rb1
		


def getTopLineGivenSegment(line_data,segment):
	for l in range(len(line_data)):
		#print "Looking at ",line_data[l]
		pieces=line_data[l].split("\t")
		line_segment=pieces[0]
		#print "The line segment is ",line_segment
		if(line_segment==segment):
			#print "To return ",line_data[l]
			return line_data[l]
	return None



def writeModedHistogramFile(IncrementMapWrapper_count_map_dict,outfile):
	max_list=list()
	for d in IncrementMapWrapper_count_map_dict:
		#print "have to write for ",d
		max_d=IncrementMapWrapper_count_map_dict[d].getMaxKeyWithIntOrdering()
		#print "for ",d," the max is ",max_d
		max_list.append(max_d)
	max_overall=max(max_list)
	hist_writer=open(outfile,'w')
	hist_writer.write("CDR3_LENGTH")
	keys=IncrementMapWrapper_count_map_dict.keys()
	keys.sort()
	for d in keys:
		hist_writer.write("\t"+d)
	keys.sort()
	hist_writer.write("\n")
	cdr3_len=(-1)
	while(cdr3_len<=max_overall):
		if(cdr3_len!=0):
			hist_writer.write(str(cdr3_len))
			for d in keys:
				hist_writer.write("\t"+str(IncrementMapWrapper_count_map_dict[d].get_val(cdr3_len)))
			hist_writer.write("\n")
		cdr3_len+=1
	hist_writer.close()


def trimList(l):
	for idx, val in enumerate(l):
		#print idx, val
		l[idx]=l[idx].strip()
	return l





def scanOutputToVDJML(input_file,fact):
	INPUT=open(input_file,'r')
	line_num=0
	current_query=None
	current_query_id=(-1)
	vdjr_vals_list=list()
	junction_vals_list=list()
	summary_vals_list=list()
	hit_vals_list=list()
	#ddl=getListOfDescsFromDBList(db_fasta_list)
	#vlist=ddl[0]
	#dlist=ddl[1]
	#jlist=ddl[2]
	meta=None
	firstV=None
	firstJ=None
	igblast_rec_id=0
	for igblast_line in INPUT:
		line_num+=1
		igblast_line=igblast_line.strip()
		#print igblast_line
		if(re.compile('^#').match(igblast_line)):
			rem=re.search('^#\s+IGBLASTN\s([^\s]+)\s*$',igblast_line)
			if rem:
				igblast_version=rem.group(1)
				mode="IGB_VERSION"
			rem=re.search('^#\s+Query:\s(.*)',igblast_line)
			if(rem):# or ref):
				#reset values for new IGBLAST set of results
				current_query=rem.group(1)
				vdjr_vals_list=list()
				junction_vals_list=list()
				summary_vals_list=list()
				hit_vals_list=list()
				mode="query_retrieve"
				current_query_id+=1
				igblast_rec_id+=1
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
					headers=re.sub(r'\s*,\s*',',',headers)
					headers=re.sub(' ','_',headers,0)
					headers=re.sub('\-','_',headers,0)
					rearrangment_summary_fields=headers.split(',')
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
			rem=re.search('^#\s(\d+)\shits\sfound\s*',igblast_line)
			if rem:
				mode="hit_table"
				num_hits_found=int(rem.group(1))
				current_num_hits=num_hits_found
				if(num_hits_found==0):
					#since there are zero hits, do a serialization now.
					#there's no hit data to pickup
					#CALL SERIALIZER as the hits for this have been picked up
					getMapToo=True
					serialized=vdjml_read_serialize(
						current_query,					#current query name
						fact, 						#handle to VDJML factory
						hit_vals_list,					#hit values (list of lists)
						hit_fields,					#list of hit fields
						summary_fields,					#summary fields
						summary_vals_list,				#summary values list (list of lists)
						rearrangement_summary_fields,vdjr_vals_list,	#rearrangment values (list and list-of-lists)
						junction_fields,junction_vals_list,		#junction fields/values (list and list of lists)
						)
					rb1=serialized
					yield rb1.release()
		elif (not (re.compile('^\s*$').match(igblast_line))):
			if(mode=="vdj_rearrangement"):
				vdjr_vals_list.append(igblast_line)
			elif(mode=="vdj_junction"):
				junction_vals_list.append(igblast_line)
			elif(mode=="vdj_alignment_summary"):
				summary_vals_list.append(igblast_line)
			elif(mode=="hit_table"):
				hit_vals_list.append(igblast_line)
				if(len(hit_vals_list)==current_num_hits):
					#CALL SERIALIZER as the hits for this read/query been picked up
					getMapToo=True
					serialized=vdjml_read_serialize(
						current_query,					#current query name
						fact, 						#handle to VDJML factory
						hit_vals_list,					#hit values (list of lists)
						hit_fields,					#list of hit fields
						summary_fields,					#summary fields
						summary_vals_list,				#summary values list (list of lists)
						rearrangement_summary_fields,vdjr_vals_list,	#rearrangment values (list and list-of-lists)
						junction_fields,junction_vals_list,		#junction fields/values (list and list of lists)
						)
					rb1=serialized
					yield rb1.release()
			else:
				print "unhandled mode=",mode



#given the arguments, make a VDJML
def makeVDJMLFileFromArgs(args):
	#build a factory for initialization
	meta = vdjml.Results_meta()
	fact=vdjml.Result_factory(meta)
	fact.set_default_aligner(
			"IgBLAST",			#aligner name
			args.igblast_version,		#aligner ver
			args.igblast_params,		#aligner params
			args.igblast_runid)		#run id
	fact.set_default_gl_database(
		args.db_name,		#db name
		args.db_ver,		#db ver
		args.db_species,	#db species
		args.db_uri		#db uri
		)
	dom_class=args.igblast_dc
	if(dom_class=="imgt"):
		#print "setting imgt for factory...."
		fact.set_default_num_system(vdjml.Num_system.imgt)
	elif(dom_class=="kabat"):
		fact.set_default_num_system(vdjml.Num_system.kabat)
	#create a result store and begin to write output
	#result_store = vdjml.Result_store(meta)
	print "vdjml is ",str(args.vdjml[0])
	rrw = vdjml.Result_writer(str(args.vdjml[0]), meta)
	for read_result in scanOutputToVDJML(args.igblast_in[0],fact):
		#result_store.insert(read_result)
		rrw(read_result)
		


#return an argument parser
def makeParserArgs():
	parser = argparse.ArgumentParser(description='Parse IGBLAST output, write a VDJML output file')
	parser.add_argument('igblast_in',type=str,nargs=1,help="path to igblast analysis file of results, hits, segments, etc. *** NOTE *** SHOULD BE RUN WITH OUTFMT AS SEEN HERE -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'")
	parser.add_argument('vdjml',type=str,nargs=1,help="the path to the output vdjml")
	parser.add_argument('-igblast_version',type=str,nargs=1,default="2.2.27+",help="the version of igblast used")
	parser.add_argument('-igblast_params',type=str,nargs=1,default="IG_BLAST_PARAMETERS",help="the parameters passed to igblast")
	parser.add_argument('-igblast_runid',type=int,nargs=1,default=0,help="the ID assigned to the run")
	parser.add_argument('-db_name',type=str,nargs=1,default="DB_NAME",help="the name of the IgBLAST database")
	parser.add_argument('-db_ver',type=str,nargs=1,default="DB_VER",help="a version string associated with the database")
	parser.add_argument('-db_species',type=str,nargs=1,default="human",help="species of the db")
	parser.add_argument('-db_uri',type=str,nargs=1,default="http://www.vdjserver.org",help="a URI associated with the database")
	parser.add_argument('-igblast_dc',type=str,nargs=1,default="imgt",help="the domain classification system used",choices=["imgt","kabat"])	
	return parser





if (__name__=="__main__"):
	parser=makeParserArgs()
	args = parser.parse_args()
	if(args): 
		makeVDJMLFileFromArgs(args)
	else:
		#print "error in args!"
		parser.print_help()

