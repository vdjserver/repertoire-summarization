#!/usr/bin/env python

from bs4 import BeautifulSoup
from collections import defaultdict
import pickle
import urllib2
import re
import pprint
from imgt_utils import get_loci_list,imgt_db
import glob
#from utils import biopythonTranslate
from utils import *
from alignment import alignment
import os





#pretty print a tree
def prettyPrintTree(t):
	pp = pprint.PrettyPrinter()
	pp.pprint(dicts(t))
	

#from a URL of an IMGT gene table generate a dict() of dict() of dict()
#that indicates the hierarchy
#have a list of ALLOWABLE names that (if not None) acts as a 'white' list so that items must be in it to be returned in the hierarchy
#have an existing hierarchy that is added to (if not None)
#return an array of TWO elements [tree,cloneNamesMap]
#the tree is the resulting hierarchy
#the clone names maps allele names to clones (may have information in MOUSE cases)
def hierarchyTreeFromGenetableURL(url,locus,listOfAllowableNames=None,existingHierarchy=None):
	f = urllib2.urlopen(url)
	html=f.read()
	soup=BeautifulSoup(html)
	clone_names_map=dict()
	table_div=soup.find(id="genetable")
	if table_div is None:
		if(existingHierarchy==None):
			return [tree(),clone_names_map]
		return [existingHierarchy,clone_names_map]
	gene_table=table_div.find('table')
	if gene_table is None:
		if(existingHierarchy==None):
			return [tree(),clone_names_map]
		return [existingHierarchy,clone_names_map]
	rows = gene_table.findAll('tr')
	if rows is None:
		if(existingHierarchy==None):
			return [tree(),clone_names_map]
		return [existingHierarchy,clone_names_map]
	subgroup=""
	gene_name=""
	allele_name=""
	group=locus

	if existingHierarchy is None:
		my_hier=tree()
	else:
		my_hier=existingHierarchy
	for tr in rows:
		#print "\tEXPLORING A TR...."
		cols = tr.findAll('td')
		col_num=0
		col_note_num=0
		if not (cols is None):
			for td in cols:
				#print "\t\tEXPLORING A TD (col_num=",col_num,")... type=",type(td)
				if('class' in td.attrs):
					#print "class in td.attrs"
					class_val=td.attrs['class']
					#print "The class value is ",class_val
					class_val_first=str(class_val[0])
					#print "class_val_first is ",class_val_first
					class_val_first_str=str(class_val_first)
					if(class_val_first_str=="subgroup_middle_note" or class_val_first_str=="subgroup_note"):
						text = td.find(text=True)# + ';'
						if text is None:
							text="UNDEFINED_SUBGROUP;"
						else:
							text+=";"
						#print "GOT A SUBGROUP=",text
						subgroup=str(text)
						subgroup=removeTerminatingSemicolonIfItExists(subgroup)
						subgroup=subgroup.strip()
					elif(class_val_first_str=="gene_note"):
						text = td.find(text=True)# + ';'
						if text is None:
							text="UNDEFINED_GENE;"
						else:
							text+=";"
						#print "GOT A GENE=",text
						gene_name=str(text)
						gene_name=removeTerminatingSemicolonIfItExists(gene_name)
						gene_name=gene_name.strip()
					elif(class_val_first_str=="col_note"):
						if(col_note_num==12):
							text=td.find(text=True)
							#print "col_note_num =",col_note_num,"and text=",text
							if(text is not None):
								clone_names_map[allele_name]=re.sub(r'\s+','',text)
								#print "clone_names_map[",allele_name,"] =",clone_names_map[allele_name]
						col_note_num+=1
					elif(class_val_first_str=="allele_note"):
						text = td.find(text=True)# + ';'
						if not (text is None):
							#print "GOT AN ALLELE=",text
							allele_name=str(text)
							alleleEndRegex=re.compile(r'\*\d+$')
							if(not(alleleEndRegex.search(allele_name))):
								#if there is no allele, make it allele *01
								allele_name+="*01"
							allele_name=allele_name.strip()
							#SOMETIMES SUBGROUP ISN'T A COLUMN IN THE TABLE SO INHERIT IT FROM THE GENE NAME!
							#if(subgroup==""):
							#	print "Setting a blank subgroup to equal gene_name=",gene_name
							#	subgroup=gene_name
							#elif(not(gene_name.startswith(subgroup))):
							#	print "Setting an old subgroup to equal gene_name=",gene_name
							#	subgroup=gene_name
							if(listOfAllowableNames==None or ((listOfAllowableNames!=None) and (allele_name in listOfAllowableNames))):
								#print "proceeding with setting : gene_name="+gene_name+", allele_name="+allele_name+", and subgroup="+subgroup
								if(subgroup=="" or (subgroup==gene_name and not(subgroup==""))):
									#if there is no subgroup or 
									#its the same as the gene (and non-empty), then store this way
									my_hier[gene_name][allele_name]
								else:
									#otherwise, store this way, using the complete path/resolution
									my_hier[subgroup][gene_name][allele_name]
							else:
								print "NOTE : skipping the addition of ",allele_name," into the hierarchy  (url=",url," and locus=",locus,")  cause it's not in the list of alleles from the Reference Fasta!"
	return [my_hier,clone_names_map]
	









#these two defs (tree and dicts) from https://gist.github.com/hrldcpr/2012250
#"One-line Tree in Python"
def tree(): return defaultdict(tree)
def dicts(t): return {k: dicts(t[k]) for k in t}




#given a base dir for a database
#and any specified organism and/or locus
#scan the reference directory set FASTAs
#and see what entries of them are in the gene tables (or not)
#present a "side-by-side" diff showing any differences thus 
#allowing user to best know how to "patch" the gene 
#tables so that all fasta entries have a place!
#Also, obtain clone name data
#Finally, return all the hierarchy and clone data!
def analyze_download_dir_forVDJserver(base_dir,countsMap=None,specifiedOrganism=None,specifiedLocus=None):
	myDB=imgt_db(base_dir)
	organisms=myDB.getOrganismList()
	org_hierarchy=tree()
	clone_names_by_org=dict()
	if(not(os.path.exists(base_dir) or not(os.path.isdir(base_dir)))):
		#if it's not a path and not a directory, then skip!
		print "WARNING, DIRECTORY TO ANALYZE,",base_dir,"DOES NOT EXIST!, SKIPPING ANALYSIS!"
		return org_hierarchy
	for organism in organisms:
		clone_names_by_org[organism]=dict()
		loci=get_loci_list()
		org_hierarchy[organism]=tree()
		for locus in loci:
			if((specifiedOrganism is not None) and (specifiedLocus is not None)):
				if((not(locus==specifiedLocus)) or (not(organism==specifiedOrganism))):
					continue
			print "Analyzing for o=",organism," and l=",locus
			###########################################
			###########################################
			# HERE , FIND THE CORRESPONDING HTML FILES FOR THE GENE TABLES
			htmlGeneTablesGlob=base_dir+"/"+organism+"/GeneTables/*"+locus+"*.html"
			geneTableHTMLFiles=glob.glob(htmlGeneTablesGlob)
			geneTableHTMLFiles.sort() #sort here so that the HTML file is first and the orphon file is second
			###########################################
			###########################################
			# HERE , FIND THE CORRESPONDING FNA FILES in the Reference Directory Set Data
			fastaFilesGlob=base_dir+"/"+organism+"/ReferenceDirectorySet/*"+locus+"*.fna"
			fastaFiles=glob.glob(fastaFilesGlob)
			##################################################
			##################################################
			#HERE extract the IMGT-style names (e.g. IGHV4-1*02)				
			fastaMainMap=dict()
			for fastaFile in fastaFiles:
				fastaString=readFileIntoString(fastaFile)
				fastaList=read_fasta_string(fastaString)
				fastaMap=read_fasta_into_map(fastaList)	
				for fastaKey in fastaMap:
					fastaMainMap[fastaKey]=fastaMap[fastaKey]
			#print "done loading into main map..."
			#printMap(fastaMainMap)
			fastaListOfNames=getIMGTNameListFromFastaMap(fastaMainMap)
			fastaListOfNames.sort()
			##################################################
			##################################################
			#HERE "allelify" the names (add "*01") if no allele information is on the name
			#this is because the desired graphics/hierarchy requires an allele at the leaves of the hierarchy
			#and because I've seen a few cases of names without alleles!!!! 
     			fastaAlleleList=allelifyList(fastaListOfNames)
			#################################################
			#################################################
			#HERE, actually load into a "tree" dict the hierachy structures
			#also obtain the clone names from the HTML
			#This is done for regular genes and for orphons
			html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[0],locus,fastaAlleleList)
			locusHierarchyData=html_data[0]
			clone_names_plain=html_data[1]
			html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[1],locus,fastaAlleleList,locusHierarchyData)
			locusHierarchyData=html_data[0]
			#print "GOT HIERARCHY FROM HD0 ORPH:"
			#prettyPrintTree(locusHierarchyData)
			#print "DONE SHOWING PRETTY PRINT GOT HIERARCHY FROM HD0 ORPH"
			clone_names_orph=html_data[1]	
			clone_names=merge_maps(clone_names_plain,clone_names_orph)
			write_map_to_file(clone_names,geneTableHTMLFiles[0]+".clone_names.map")
			#print "GOT FASTA MAP KEY ALLELES :"
			#printList(fastaAlleleList)
			#print "GOT HIERARCHY:"
			#prettyPrintTree(locusHierarchyData)
			#print "DONE SHOWING PRETTY PRINT GOT HIERARCHY"
			#print "SHOWING CLONE NAME MAP:"
			#printMap(clone_names)
			print "\n\n\n\n"
			#from
			#http://stackoverflow.com/questions/82831/how-do-i-check-if-a-file-exists-using-python
			patchPath=geneTableHTMLFiles[0]+".patch"
			foundPatchFile=False
			if(os.path.isfile(patchPath)):
				foundPatchFile=True
				print "Found a patch file (",patchPath,") found, so now patching is being performed....."
				print "\n\n############################################"
				print "BEFORE PATCHING, THE HIERARCHY IS  :"
				prettyPrintTree(locusHierarchyData)
				locusHierarchyData=patchlocusHierarchyData(locusHierarchyData,patchPath)
				print "\n\n############################################"
				print "AFTER PATCHING, THE HIERARCHY IS  :"
				prettyPrintTree(locusHierarchyData)
				#print "THIS IS THE PRETTY PRINT FOR PATCHED ORG HIERARCHY DATA"
				org_hierarchy[organism][locus]=locusHierarchyData
			else:
				print "No patch file (",patchPath,") found, so no patching being performed....."
			#################################################
			#################################################
			#extract the leaves/alleles from the tree/hierarchy structure and perform a "side-by-side" diff			
			treeAlleles=get_list_of_alleles_appearing_in_tree(locusHierarchyData)
			setSameStats=(set(fastaAlleleList) == set(treeAlleles))
			print "For organism=",organism,"locus=",locus
			#################################################
			#################################################
			#based on the side-by-side diff result and based
			#on the patching result, 
			if setSameStats:
				print "The fasta/RefDirNames ARE the same as the hierarchy allele names!!! :)"
				if(foundPatchFile):
					#no need to re-print the hierarchy cause it was printed before and after patching
					pass
				else:
					#since no patch file was found, go ahead and print the hierarchy
					print "\n\n############################################"
					print "The hierarchy is : "
					prettyPrintTree(locusHierarchyData)						
			else:
				print "The fasta/RefDirNames ARE different from the hierarchy allele names!!! :("
				if(not(foundPatchFile)):
					#since no patch file was found, AND the stats are bad, print the hierarchy to 
					#help the user understand how to patch					
					print "\n\n############################################"
					print "The hierarchy is : "
					prettyPrintTree(locusHierarchyData)				
			briefSetDiff(fastaAlleleList,treeAlleles,"fasta alleles "+organism+"_"+locus,"tree alleles "+organism+"_"+locus)			
			extendedSortedSetDiff(fastaAlleleList,treeAlleles,"fasta alleles "+organism+"_"+locus,"tree alleles "+organism+"_"+locus)
			#print "THIS IS THE PRETTY PRINT FOR LOCUS HIERARCHY DATA, o=",organism,"l=",locus
			#prettyPrintTree(locusHierarchyData)
			#print "THIS IS THE PRETTY PRINT FOR ORG HIERARCHY DATA"
			org_hierarchy[organism][locus]=locusHierarchyData
			#prettyPrintTree(org_hierarchy)
			clone_names_by_org[organism]=merge_maps(clone_names_by_org[organism],clone_names)
			print "\n\n\n"
	#prettyPrintTree(org_hierarchy)
	pickleFilePath=base_dir+"/hierarchy_data.pkl"
	tbr=[org_hierarchy,clone_names_by_org]
	pickleWrite(pickleFilePath,tbr)
	return tbr







#given a gene table hierarchy and a patch file path
#patch the hierarchy
#format is add H1 H2,....HN-1,HN
#format is del H1 H2,....HN-1,HN
#to add and remove entries
#NOTE that the hierarchy and path must be at the same level (organism->locus)
def patchlocusHierarchyData(locusHierarchyData,patchfilepath):
	patch_file=open(patchfilepath,'r')
	joinStr="']['"
	for line in patch_file:
		line=line.strip()
		pieces=line.split('\t')
		if(re.search("'",line)):
			hasApos=True
		else:
			hasApos=False
		if(len(pieces)>=2 and not(hasApos)):
			afterCmd=pieces[1:]
			joinedStr=joinStr.join(afterCmd)
			stringToExec="locusHierarchyData['"
			if(pieces[0]=="add"):
				stringToExec=stringToExec+joinedStr+"']"
			elif(pieces[0]=="del"):
				stringToExec="del "+stringToExec+joinedStr+"']"
			if(not(stringToExec=="locusHierarchyData['")):
				#print "GOT A COMMAND TO EXECUTE : ",stringToExec
				exec(stringToExec)
	return locusHierarchyData







#write an object to a pickle file
def pickleWrite(pPath,o):
	output = open(pPath, 'wb')
	# Pickle dictionary using protocol 0.
	pickle.dump(o, output)
	output.close()







#read an object from a pickle file
def pickleRead(pPath):
	pkl_file = open(pPath, 'rb')
	data1 = pickle.load(pkl_file)
	pkl_file.close()
	return data1






#from a tree get the total number of counts
#using the counts map	
def get_total_tree(tree,this_name,counts_map):
	subtree_sum=0
	if this_name in counts_map:
		subtree_sum+=counts_map[this_name]
	for child in tree:
		subtree_sum+=get_total_tree(tree[child],child,counts_map)
	#print "For subtree rooted at :",this_name,"returning the sum",subtree_sum
	return subtree_sum






#zero pad a string e
#up to length tot
def zeropad(e,tot):
	str_e=str(e)
	if(len(str_e)>=tot):
		return str_e
	else:
		numZ=tot-len(str_e)
		zeroes=""
		for n in range(numZ):
			zeroes+="0"
		return zeroes+str_e


#zero pad a string (of digits)
#up to a length
def zeroPadDigitsTo(s,num=8):
	if(re.search(r'\d',s)):
		#proceed here
		reg=re.compile(r'\d+')
		digits_list=reg.findall(s)
		nondigits_list=reg.split(s)
		#print "got digits list : ",digits_list
		#print "got nondigits list : ",nondigits_list
		ret=""
		for i in range(len(nondigits_list)-1):
			ret+=nondigits_list[i]
			ret+=zeropad(digits_list[i],num)
		if(len(nondigits_list)>len(digits_list)):
			ret+=nondigits_list[len(nondigits_list)-1]
		#uniq_map=dict()
		#for digits in digits_list:
		#	uniq_map[digits]=digits
		#digits_list=uniq_map.keys()
		#for digits in digits_list:
		#	existing=str(digits)
		#	padded=zeropad(existing,num)
		#	s=re.sub(digits,padded,s)
		return ret
	else:	
		return s




#JSONIFY from a hierarchy with a counts map
#the "count_string" is used to make the label for the count in the JSON
#the "labelString" is used to make the label for the names
def jsonify_hierarchy(hier_map,name,counts_map,count_string,labelString="label"):
	#print "\n\n\n\n\n"
	hier_map_keys=hier_map.keys()
	original_to_padded=dict()
	padded_to_original=dict()
	padded_key_list=list()
	for hier_map_key in hier_map_keys:
		padded_key=zeroPadDigitsTo(hier_map_key)
		original_to_padded[hier_map_key]=padded_key
		padded_to_original[padded_key]=hier_map_key
		padded_key_list.append(padded_key)
	#print "A key list is :",hier_map_keys
	#print "A padded list is :",padded_key_list
	padded_key_list.sort()
	#print "A sorted padded list is :",padded_key_list
	#print "\n\n\n"
	JSON=""
	JSON+="{\n"
	JSON+="\""+labelString+"\":\""+name+"\",\n"
	actual_value=0
	actual_value=get_total_tree(hier_map,name,counts_map)
	#JSON+="\"value\":"+str(actual_value)+",\n"
	#JSON+="\"size\":"+str(actual_value)+"\n"
	JSON+="\""+count_string+"\":"+str(actual_value)+"\n"
	num_kids=len(hier_map)
	if(num_kids>=1):
		JSON+=",\n"
		JSON+="\"children\": [\n"
		kid_num=0
		#for child in hier_map:
		for c in range(len(hier_map_keys)):
			#child=hier_map_keys[c]
			child=padded_key_list[c]
			original_child=padded_to_original[child]
			JSON+=jsonify_hierarchy(hier_map[original_child],original_child,counts_map,count_string)
			if(kid_num<(num_kids-1)):
				JSON+=" , \n"
			kid_num+=1
		JSON+=" ] \n"
	else:
		JSON+=""
		#num_kids<1
	JSON+=" }\n"
	return JSON


#from a BLAST hit string, attempt to parse out 
#the name the allele and use that to return
#the descriptor from the fastalistofdescs
def getSegmentName(blast_hit,fastaListOfDescs):
	gnl_re=r'^(gnl\|)?BL_ORD_ID[\|:](\d+)$'
	if(blast_hit in fastaListOfDescs):
		return blast_hit.strip()
	elif(re.match(gnl_re,blast_hit)):
		reo=re.match(gnl_re,blast_hit)
		digits=reo.group(2)
		index=int(digits)
		return fastaListOfDescs[index].strip()
	else:
		raise Exception("Error, unhandled blast segment rearrangment descriptor '"+blast_hit+"' !") 




#return a list of descriptors from a file
def getFastaListOfDescs(fasta_file):
	fr=open(fasta_file,'r')
	descs=list()
	for line in fr:
		if(line.startswith(">")):
			to_append=line[1:]
			to_append=to_append.strip()
			descs.append(to_append)
	return descs





#given an IGBLAST output file, scan it and extract the name
#from the list of fasta files
#and increment counts based on the extractions
def get_segment_count_map_from_blast_output(blast_out,fasta_file_list):
	vlist=getFastaListOfDescs(fasta_file_list[0])
	dlist=getFastaListOfDescs(fasta_file_list[1])
	jlist=getFastaListOfDescs(fasta_file_list[2])
	blast_data=open(blast_out,'r')
	line_num=0
	counts_map=dict()
	capture_flag=False
	is_heavy=True
	for line in blast_data:
		if(capture_flag):
			pieces=line.split('\t')
			if(is_heavy):
				max_range=3
			else:
				max_range=2
			for i in range(max_range):
				pieces[i]=pieces[i].strip()
				list_to_use=vlist
				if(i==1 and is_heavy):
					list_to_use=dlist
				elif(i==1 and not is_heavy):
					list_to_use=jlist
				elif(i==2):
					list_to_use=jlist
				if(not(pieces[i]=="N/A")):
					segment=getSegmentName(pieces[i],list_to_use)
					#print "GOT A SEGMENT",segment
					if(segment in counts_map):
						counts_map[segment]+=1
					else:
						counts_map[segment]=1
			capture_flag=False
		if(line.startswith("# V(D)J rearrangement")):
			if(line.find('Top D')!=(-1)):
				is_heavy=True
			else:
				is_heavy=False
			capture_flag=True
		line_num+=1
		#print "at line="+str(line_num)
	return counts_map



#given the JALLELE name and IMGTDB object (and organism and mode)
#find the position (1-based) of the CDR3 end in the indicated sequece
#use the LOOKUP tables in the database 
#IMGT data use the imgt.dat .  Once data are looked up, store then in a
#dict structure for subsequent faster lookups
cdr3_end_lookup=dict()
def getADJCDR3EndFromJAllele(jallele,imgtdb_obj,org="human",mode="imgt"):
	#print "in getADJCDR3EndFromJAllele. jallele=",jallele," org=",org,"mode=",mode
	global cdr3_end_lookup
	if(org in cdr3_end_lookup):
		if(jallele in cdr3_end_lookup[org]):
			if(mode in cdr3_end_lookup[org][jallele]):
				#the lookup is there!  use it!
				return cdr3_end_lookup[org][jallele][mode]
			else:
				#the lookup ain't there, look it up , then cache it!
				pass
		else:
			#add the allele if it ain't there
			cdr3_end_lookup[org][jallele]=dict()		
	else:
		#add the organism if it ain't there
		#add the jallele cause it won't be there either
		cdr3_end_lookup[org]=dict()
		cdr3_end_lookup[org][jallele]=dict()
	val_to_return=(-1)
	if(mode=="imgt"):
		#print "mode=imgt"
		jdescriptor=imgtdb_obj.extractDescriptorLine(jallele,org)
		cdr3_end_raw=getCDR3EndFromJData(jallele,imgtdb_obj,org)
		if(cdr3_end_raw==None):
			cdr3_end_lookup[org][jallele][mode]=(-1)
			return cdr3_end_lookup[org][jallele][mode]
		cdr3_end_raw=int(cdr3_end_raw)
		#print "Got a descriptor ",jdescriptor
		#print "got a raw cdr3 end=",cdr3_end_raw
		desc_pieces=jdescriptor.split("|")
		if(desc_pieces[14]=="rev-compl"):
			desc_pieces[5]=swapIMGTDescInterval(desc_pieces[5])
			sep="|"
			jdescriptor=sep.join(desc_pieces)
		#print "got a corrected jdescriptor ",jdescriptor
		interval_re=re.compile(r'(\d+)\.+(\d+)')
		desc_range=desc_pieces[5]
		mr=re.search(interval_re,desc_range)
		if(mr):
			start=int(mr.group(1))
			stop=int(mr.group(2))
			#print "got start=",start
			#print "got stop=",stop
			if(min(start,stop)<=cdr3_end_raw and cdr3_end_raw<=max(start,stop)):
				adjusted=cdr3_end_raw-min(start,stop)+1
				#print "RETURNING VALID CDR3 END = "+str(adjusted)
				cdr3_end_lookup[org][jallele][mode]=adjusted
			else:
				#print "OUT OF RANGE????"
				cdr3_end_lookup[org][jallele][mode]=(-1)
		else:
			#print "FAILED TO MATCH ON INTERVAL REGEX!!!!\n"
			cdr3_end_lookup[org][jallele][mode]=(-1)
		return cdr3_end_lookup[org][jallele][mode]
	elif(mode=="kabat"):
		kabat_file=imgtdb_obj.getBaseDir()+"/"+org+"/ReferenceDirectorySet/KABAT/Jlookup.tsv"
		reader=open(kabat_file,'r')
		cdr3_end_adj_kabat=(-1)
		for line in reader:
			line=line.strip()
			if(line.startswith(jallele)):
				pieces=line.split("\t")
				cdr3_end_adj_kabat=int(pieces[1])
		reader.close()
		val_to_return=cdr3_end_adj_kabat
	cdr3_end_lookup[org][jallele][mode]=val_to_return
	return 	cdr3_end_lookup[org][jallele][mode]
		



#use the IMGT dat FILE to find the CDR3 end for a JALLELE
def getCDR3EndFromJData(allele,imgtdb_obj,org="human"):
	#print "*********************************\nEntering getCDR3EndFromJData\n*******************************"
	biopyrec=imgtdb_obj.getIMGTDatGivenAllele(allele,True,org)
	records=[biopyrec]
	reg_start=None
	reg_end=None
	extracted_descriptor=imgtdb_obj.extractDescriptorLine(allele,org)
	#print "The extracted descriptor from inputs ",allele," and ",org," is (extracted) : "+str(extracted_descriptor)
	extracted_descriptor_pieces=imgtdb_obj.extractIMGTDescriptorPieces(extracted_descriptor)
	extracted_descriptor_interval=imgtdb_obj.getStartStopFromIMGTDescPieces(extracted_descriptor_pieces)
	reg_start=extracted_descriptor_interval[0]
	reg_end=extracted_descriptor_interval[1]

	if(not(reg_start==None)):
		for record in records:
			feature_list=record.features
			for feature in feature_list:
				#print "got a feature : ",feature
				ftype=feature.type
				#print "the type is ",ftype	
				qualifiers=feature.qualifiers
				#print "qualifiers : ",qualifiers
				if(ftype=="J-TRP" or ftype=="J-PHE"):
					c_start=int(re.sub(r'[^0-9]','',str(feature.location.start)))
					c_end=int(re.sub(r'[^0-9]','',str(feature.location.end)))
					#print "found a jtrp"
					if(reg_start<=c_end and c_end<=reg_end):
						#this is it!
						#subtract 3 becuase we want to ignore the TRP or PHE residue and look at the residue immediately preceding
						return c_end-3
	else:
		print "failed to get a start!"
		pass
	#os.remove(tmp_file_path)
	print "getCDR3EndFromJData returning 'None' for allele="+str(allele)+" and org="+str(org)
	return None


#swap the start/stop in an interval descriptor 123..456 -> 456..123
def swapIMGTDescInterval(i):
	interval_re=re.compile(r'(\d+)(\.+)(\d+)')
	mr=re.search(interval_re,i)	
	if(mr):
		p1=mr.group(1)
		p2=mr.group(3)
		dots=mr.group(2)
		return p2+dots+p1
	else:
		raise Exception("Error, bad interval in IMGT descriptor "+i)






cdr3_adj_map=dict()
cdr3_adj_map["human"]=dict()
cdr3_adj_map["Mus_musculus"]=dict()
cdr3_adj_map["human"]["imgt"]=dict()
cdr3_adj_map["human"]["kabat"]=dict()
cdr3_adj_map["Mus_musculus"]["kabat"]=dict()
cdr3_adj_map["Mus_musculus"]["imgt"]=dict()
#return the CDR3 start using a cache mechanism
#LOOKUPs are used
#IMGT-HIGHVQUEST for IMGT
#KABAT TSV values from IGBLAST for KABAT
def getAdjustedCDR3StartFromRefDirSetAllele(allele,imgtdb_obj,organism="human",mode="imgt"):
	#use this map as a caching mechanism!
	global cdr3_adj_map
	if(mode=="kabat" and alleleIsTR(allele)):
		#no kabat for CDR3 for TR!
		return (-1)
	if(organism in cdr3_adj_map):
		if(allele in cdr3_adj_map[organism]["imgt"]):
			return cdr3_adj_map[organism]["imgt"][allele]
	else:
		cdr3_adj_map[organism]=dict()
	cdr3_int=getVRegionStartAndStopGivenRefData(allele,organism,imgtdb_obj,"CDR3",mode)
	cdr3_adj_map[organism]["imgt"][allele]=cdr3_int[0]
	return cdr3_int[0]






#given a V segment alignment extract from it the sub-portion for a given
#region.  Return a 4-length array with subject and query alignment data  [region_alignment,frame_mask,r_q_start,r_q_end] (region_alignment has query first, then sub)
#return None if alignment too short or doesn't cover region
def getRegionAlignmentFromLargerVAlignment(sub_info_map,org,mode,region_name,imgtdb_obj,wholeOnly=False):
	#print "THE SIM is ",sub_info_map
	if(wholeOnly==True):
		raise Exception("Error, wholeOnly not yet implemented!!!!!!!!!")
	#print "\n\n\n\nusing region=",region_name
	valid_regions=getVRegionsList()
	if(not(region_name in valid_regions)):
		print region_name," is an invalid region!!!!"
		return None
	sub_name=sub_info_map['subject ids']
	ref_region_interval=getVRegionStartAndStopGivenRefData(sub_name,org,imgtdb_obj,region_name,mode)
	ref_region_transcript_start=getVRegionStartAndStopGivenRefData(sub_name,org,imgtdb_obj,region_name,"imgt")[0]
	#print "For reference=",sub_name," for org=",org," region=",region_name," got (mode=",mode,")region=",ref_region_interval
	if(ref_region_interval[0]==(-1) or ref_region_interval[1]==(-1)):
		#print "Can't get region info, start or stop of ref region is neg 1...."
		return None
	else:
		#create an alignment from the whole V
		aln_obj=alignment(
			sub_info_map['query seq'],
			sub_info_map['subject seq'],
			sub_info_map['q. start'],
			sub_info_map['q. end'],
			sub_info_map['s. start'],
			sub_info_map['s. end']
			)
		#extract a sub-alignment
		region_aln=aln_obj.getSubAlnInc(ref_region_interval[0],ref_region_interval[1],"subject")
		#aquire frame information
		region_frame_start=(ref_region_interval[0]-ref_region_transcript_start)%3
		#use/set the frame information in the alignment object
		region_aln.setSFM(region_frame_start)
		return region_aln











	
#get an empty map as a placeholder
def getEmptyRegCharMap():
	char_map=dict()
	char_map['insertions']=0
	char_map['deletions']=0
	char_map['base_sub']=0
	char_map['synonymous_bsb']=0
	char_map['nonsynonymous_bsb']=0
	char_map['mutations']=0
	char_map['pct_id']=0
	char_map['bsb_freq']=0
	char_map['ns_rto']=0
	char_map['stp_cdn']=False
	return char_map


#can be 'None', True, or False
#uses frame information from begining of J alignment, end of V alignment to determine
#if the arrangment is "productive".  It's productive if the first J nuc
#is in consistent frame with the last V nuc
#If no CDR3 is found then NONE is returned
#otherwise, get the frames where V/J end/begin aligning (using the CDR3 start and ends)
#and verify that the frames are consistent (ie if no BP in between the frames are 0->1, if 1, then 0->2, if 2, then 0->0, if 3 0->1, and so on
def getNAProductiveRearrangmentFlagFromVJHitLineData(firstVMap,firstJMap,organism,imgtdb_obj):
	###########################################
	#some code to look at the 'productive rearrangement
	productive_flag=None
	if(firstVMap!=None and firstJMap!=None):
		v_s_fr3_imgt_start=getVRegionStartAndStopGivenRefData(firstVMap['subject ids'],organism,imgtdb_obj,"FR3","imgt")[0]
		if(v_s_fr3_imgt_start!=(-1)):
			v_s_end=int(firstVMap['s. end'])
			v_s_start=int(firstVMap['s. start'])
			v_s_cdr3Start=getAdjustedCDR3StartFromRefDirSetAllele(firstVMap['subject ids'],imgtdb_obj,organism,"imgt")
			v_s_end_frame=getTheFrameForThisReferenceAtThisPosition(firstVMap['subject ids'],organism,imgtdb_obj,v_s_end)
			j_s_start=int(firstJMap['s. start'])
			j_s_imgt_cdr3_end=int(getADJCDR3EndFromJAllele(firstJMap['subject ids'],imgtdb_obj,organism,"imgt"))
			if(j_s_imgt_cdr3_end==(-1)):
				#can't get CDR3 end
				return productive_flag
			if(j_s_start>j_s_imgt_cdr3_end):
				#algignment to J starts after CDR3 end...
				return productive_flag
			j_s_start_frame=(2-(j_s_start-j_s_imgt_cdr3_end))%3  
			#the 2 (at the beginning) is here because the CDR3 end of J is in frame2 cause its the 3rd NA in a codon
			num_bp_between_V_and_J=int(firstJMap['q. start'])-int(firstVMap['q. end'])-1
			if(num_bp_between_V_and_J<0):
				#unhandled case where V/J overlap!!!
				v_s_end_frame=ofRightmostVNucleotideNotAligningToJGetItsFrameAssumingJunctionSegmentOverlap(vMap,jMap,v_s_end_frame)
				j_s_start_frame=ofLeftmostJNucleotideNotAligningToJGetItsFrameAssumingJunctionSegmentOverlap(vMap,jMap,v_s_end_frame)
				num_bp_between_V_and_J=num_bp_between_V_and_J*(-1)
			#print "num_bp_between_V_and_J",num_bp_between_V_and_J
			#print "v_s_end_frame=",v_s_end_frame
			#print "j_s_start_frame",j_s_start_frame
			if((v_s_end_frame+num_bp_between_V_and_J+j_s_start_frame)%3==1):
				#print "TESTED PRODUCTIVE BY ALIGNMENT"
				productive_flag=True
			else:
				productive_flag=False
		else:
			#couldn't get fr3, so can't get CDR3 either....
			return productive_flag
	return productive_flag


#of the rightmost V nucleotide in the query not also aligned to J get its frame!
def ofRightmostVNucleotideNotAligningToJGetItsFrameAssumingJunctionSegmentOverlap(vMap,jMap,frameAtVSEnd):
	#assume NO indels in a junction overlap
	q_v_end=int(vMap['q. end'])
	q_j_bgn=int(jMap['q. start'])
	overlap_size=q_v_end-q_j_bgn+1
	desired_frame=(frameAtVSEnd-overlap_size)%3
	return desired_frame


#use the frame at the end of V to get the frame when J alignment begins
def ofLeftmostJNucleotideNotAligningToJGetItsFrameAssumingJunctionSegmentOverlap(vMap,jMap,frameAtVSEnd):
	#assume no indels in a junction overlap
	desired_frame=(frameAtVSEnd+1)%3
	return desired_frame



#at the indicated position find the frame by using the
#imgt delineation/region starts as a "base" or "frame of reference"
def getTheFrameForThisReferenceAtThisPosition(refName,organism,imgtdb_obj,refPos):
	regions=getVRegionsList(True)
	for r in range(len(regions)-1):
		#print "refName=",refName
		#print "refOrg=",organism
		#print "r=",r
		#print "region=",regions[r]
		interval_r=getVRegionStartAndStopGivenRefData(refName,organism,imgtdb_obj,regions[r],"imgt")
		interval_r_next=getVRegionStartAndStopGivenRefData(refName,organism,imgtdb_obj,regions[r+1],"imgt")
		#print "interval_r=",interval_r
		#print "interval_r_next=",interval_r_next
		if(interval_r!=None and interval_r_next!=None):
			interval_r_start=int(interval_r[0])
			interval_r_stop=int(interval_r[1])
			interval_r_next_start=int(interval_r_next[0])
			interval_r_next_stop=int(interval_r_next[1])
			if(
				#interval_r_start!=(-1) and 
				interval_r_stop!=(-1) and 
				interval_r_next_start!=(-1)#  and
				#interval_r_next_stop!=(-1)
				):
				return (refPos-interval_r_next_start)%3
			else:
				pass
	eMsg="ERROR, COULD NOT FIND ANY FRAME POSITION AT ALL FOR "+refName+" with organism="+organism+"!"
	eMsg+="\nIs frame data avaialble in the database?!?!? Consider Re-run but skipping annotation ('-skip_char')"
	print eMsg
	raise Exception(eMsg)




def getVRegionsList(includeCDR3=False):
	regions=["FR1","CDR1","FR2","CDR2","FR3"]
	if(includeCDR3):
		regions.append("CDR3")
	return regions



def alleleIsTR(an):
	if(an.startswith("TR")):
		return True
	else: 
		return False



#have a cache map for region information for reference data
reg_adj_map=dict()
reg_adj_map["human"]=dict()
reg_adj_map["Mus_musculus"]=dict()
reg_adj_map["human"]["imgt"]=dict()
reg_adj_map["human"]["kabat"]=dict()
reg_adj_map["Mus_musculus"]["kabat"]=dict()
reg_adj_map["Mus_musculus"]["imgt"]=dict()
def getVRegionStartAndStopGivenRefData(refName,refOrg,imgtdb_obj,region,mode):
	global reg_adj_map
	regions=getVRegionsList(True)
	if(refOrg in reg_adj_map):
		if(mode in reg_adj_map[refOrg]):
			if(refName in reg_adj_map[refOrg][mode]):
				if(region in reg_adj_map[refOrg][mode][refName]):
					#print "using cache for ",refName," o=",refOrg," region=",region," and mode=",mode
					#print "to return", reg_adj_map[refOrg][mode][refName][region]
					#sys.exit(0)
					return reg_adj_map[refOrg][mode][refName][region]
				else:
					pass
					#region start/stop not in reg_adj_map[refOrg][mode][refName]
			else:
				#refname not in org->mode map, so add it
				#print "adding refname key=",refName
				reg_adj_map[refOrg][mode][refName]=dict()
		else:
			#mode not there! a bad mode!
			raise Exception("unknown mode ",mode," (expected kabat or imgt)!")
			#sys.exit(0)
	else:
		#organism not in there!
		raise Exception("Unknown organism "+refOrg)
	#first, form a path to a lookup.
	#this depends on organism and mode
	seq_type=refName[0:2]
	lookupPath=imgtdb_obj.getBaseDir()+"/"+refOrg+"/ReferenceDirectorySet/"+mode.upper()+"/Vlookup."+seq_type+".tsv"
	if(not(os.path.exists(lookupPath)) or not(os.path.isfile(lookupPath))):
		raise Exception("Error, failed to find lookup file "+lookupPath+" for "+refName+" for organism = "+refOrg)
	idx_num=None
	for i in range(len(regions)):
		if(regions[i]==region):
			idx_num=i
	if((idx_num is None) or not(region in regions)):
		print "ERROR, UNKNOWN REGION : ",region
		sys.exit(0)
	region_reader=open(lookupPath,'r')
	for line in region_reader:
		line=line.strip()
		line_pieces=line.split('\t')
		line_segment=line_pieces[0]
		#print "looking at ",line
		if(line_segment==refName):
			col_num=1+idx_num*2
			region_interval=[int(line_pieces[col_num]),int(line_pieces[col_num+1])]
			reg_adj_map[refOrg][mode][refName][region]=region_interval
			region_reader.close()
			return reg_adj_map[refOrg][mode][refName][region]
	region_reader.close()
	print "ERROR, FAILED TO FIND "+mode.upper()+" REGION FOR REFERENCE NAMED "+refName+" in "+lookupFile
	#return [(-1),(-1)]
	sys.exit(0)	



#class for count map
class IncrementMapWrapper():

	#have a count map
	count_map=None

	#have a default value
	default_value=0

	#constructor
	def __init__(self):
		self.count_map=dict()
		self.default_value=0

	#treating the keys as ints, find the max key
	#return -1 if no keys
	def getMaxKeyWithIntOrdering(self):
		keys=self.count_map.keys()
		int_keys=list()
		for k in keys:
			int_keys.append(int(k))
		int_keys.sort()
		if(len(int_keys)>0):
			return max(int_keys)
		else:
			return (-1)

	#be able to change the default value
	def setDefault(self,d):
		self.default_value=d

	#have a subroutine to return the map
	def get_map(self):
		return  self.count_map


	#have a getter method for the default value
	def get_val(self,key):
		if(key in self.count_map):
			return self.count_map[key]
		else:
			return self.default_value
	
	#increment the given value
	def increment(self,val):
		if(val in self.count_map):
			self.count_map[val]+=1
		else:
			self.count_map[val]=1

	#print it
	def printMap(self):
		print "I am an increment map."
		for k in self.count_map:
			print str(k)+"\t"+str(self.count_map[k])

	#zero out at a key
	def zeroOut(self,key):
		self.count_map[key]=0

	#given another increment map wrapper, merge its counts into this/self 
	def mergeInto(self,other):
		other_map=other.get_map()
		for k in other_map:
			k_count=other_map[k]
			#print "GOT k=",k," and k_count=",k_count," from other map"
			for i in range(k_count):
				self.increment(k)

	#write JSON to file
	def JSONIFYToFile(self,db_base,organism,filePath,filterbyFastaAlleles=False):
		JSON=self.JSONIFYIntoHierarchy(db_base,organism,filterbyFastaAlleles)
		writer=open(filePath,'w')
		writer.write(JSON)
		writer.close()


	#JSONIFY into hierarchy
	def JSONIFYIntoHierarchy(self,db_base,organism,filterbyFastaAlleles=False):
		hierarchy=getHierarchyByOrganism(db_base+"/"+organism+"/GeneTables/",organism,filterbyFastaAlleles)
		JSON=jsonify_hierarchy(hierarchy,organism,self.count_map,"value")
		return JSON









#the lean right means that if the subject pos aligns to a gap, then choose the bp position to the right
#otherwise if the subject aligns to a gap at the position, lean "left" and choose the previous bp positoin
#RIGHT is the default!
def getQueryIndexGivenSubjectIndexAndAlignment(query_aln,subject_aln,q_start,q_stop,s_start,s_stop,subject_pos,lean="right"):
	#print "in getQueryIndexGivenSubjectIndexAndAlignment"
	#import hashlib
	#print "query_aln=",query_aln
	#print "query_aln_digest=",hashlib.md5(query_aln).hexdigest()
	#print "subct_aln=",subject_aln
	#print "subct_aln_digest=",hashlib.md5(subject_aln).hexdigest()	
	#print "q_start=",q_start
	#print "q_stop=",q_stop
	#print "s_start=",s_start
	#print "s_stop=",s_stop
	#print "desired=",subject_pos
	firstCond=(s_start<=subject_pos)
	secondCond=(subject_pos<=s_stop)
	#print "firstCond (s_start<=subject_pos) is ",firstCond
	#print "secondCond (subject_pos<=s_stup)   is ",secondCond
	if(not(firstCond) or not(secondCond)):
		#invalid range!
		#print "INVALID RANGE!"
		#print "s_start=",s_start,"dbus=",subject_pos,"s_stop=",s_stop
		#print "getQueryIndexGivenSubjectIndexAndAlignment to return (neg) (-1)"
		return (-1)
	if(firstCond and secondCond):
		#if within the boundaries!
		s_counter=s_start
		q_counter=q_start
		delta=0
		q_val=query_aln[delta]
		s_val=subject_aln[delta]
		#print "s counter=",s_counter
		#print "s_val=",s_val
		#print "q counter=",q_counter
		#print "q_val=",q_val
		#print "alignment s, then q:"
		#print subject_aln
		#print query_aln
		while(delta<len(query_aln)):
			q_val=query_aln[delta]
			s_val=subject_aln[delta]
			#print "in while"
			#print "delta=",delta
			#print "len(query_aln[0])=",len(query_aln)
			#print "s_counter=",s_counter
			#print "q_counter=",q_counter
			#print "q_val=",q_val
			#print "s_val=",s_val
			#print "delta+s_counter="+str(delta+s_counter)
			#print "desired=",subject_pos
			if(subject_pos<=s_counter):#: and s_val!="-"):
				#print "got it"
				#return 1-based index
				#return q_counter
				if(s_val!="-"):
					if(lean=="right"):
						return q_counter
					elif(lean!="right" and q_val=="-"):
						return q_counter-1
					else:
						return q_counter
					#print "getQueryIndexGivenSubjectIndexAndAlignment to return (qcounter) ",q_counter
					#if(q_val!="-"):
					#	return q_counter
					#else:
					#	return q_counter+1
			if(not(q_val=="-")):
				q_counter+=1
			if(not(s_val=="-")):
				s_counter+=1
			delta+=1
			#print "\n\n"



#rooted at an organism get the hierarchy from the gene tables
def getHierarchyByOrganism(geneTablesDirectoryOfHTMLFiles,org_name,filterbyFastaAlleles=False):
	loci=get_loci_list()
	hierarchy=tree()
	for locus in loci:
		htmlGeneTablesGlob=geneTablesDirectoryOfHTMLFiles+"/*"+locus+"*.html"
		geneTableHTMLFiles=glob.glob(htmlGeneTablesGlob)
		fastaAlleleList=None
		if(filterbyFastaAlleles):
			fastaFilesGlob=geneTablesDirectoryOfHTMLFiles+"/../ReferenceDirectorySet/*"+locus+"*.fna"
			print "the fasta glob is ",fastaFilesGlob
			fastaFiles=glob.glob(fastaFilesGlob)
			fastaMainMap=dict()
			for fastaFile in fastaFiles:
				fastaString=readFileIntoString(fastaFile)
				fastaList=read_fasta_string(fastaString)
				fastaMap=read_fasta_into_map(fastaList)	
				for fastaKey in fastaMap:
					fastaMainMap[fastaKey]=fastaMap[fastaKey]
			print "done loading into main map..."
			fastaListOfNames=getIMGTNameListFromFastaMap(fastaMainMap)
			print "FASTA LIST OF NAMES:"
			fastaListOfNames.sort()
			printList(fastaListOfNames)
			print "got name list.... its length is",len(fastaListOfNames)
			fastaAlleleList=allelifyList(fastaListOfNames)
			print "The Fasta Allele list is ",fastaAlleleList
		html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[0],locus,fastaAlleleList)
		locusHierarchyData=html_data[0]
		html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[1],locus,fastaAlleleList,locusHierarchyData)
		locusHierarchyData=html_data[0]
		hierarchy[locus]=locusHierarchyData
	return hierarchy



if (__name__=="__main__"):
	#getQueryIndexGivenSubjectIndexAndAlignment
	#getQueryIndexGivenSubjectIndexAndAlignment(query_aln,subject_aln,q_start,q_stop,s_start,s_stop,subject_pos):
	query_aln="AC-XN-GT"
	sbjct_aln="-GAT-ACA"
	query_from=2
	query_to=7
	sbjct_f=3
	sbjct_t=8
	for i in range(3,8+1):
		print "\n\n\n"
		print "Now analyzing (s at top, q at bottom):"
		print sbjct_f,sbjct_aln,sbjct_t
		print query_from,query_aln,query_to
		print "q_from=",query_from," query_to=",query_to,",s_from=",sbjct_f,", s_to=",sbjct_t," At subject position=",i
		#print "The corresponding query position is "+str(getQueryIndexGivenSubjectIndexAndAlignment(query_aln,sbjct_aln,query_from,query_to,sbjct_f,sbjct_t,i))
		#print "The corresponding (lean_left) query position is "+str(getQueryIndexGivenSubjectIndexAndAlignment(query_aln,sbjct_aln,query_from,query_to,sbjct_f,sbjct_t,i,"left"))







