#!/usr/bin/env python

from bs4 import BeautifulSoup
from collections import defaultdict
import pickle
import urllib2
import re
import pprint
from imgt_utils import get_loci_list,imgt_db
import glob


def removeTerminatingSemicolonIfItExists(s):
	s=re.sub(r';+$',"",s)
	return s
	

def prettyPrintTree(t):
	pp = pprint.PrettyPrinter()
	pp.pprint(dicts(t))
	






def hierarchyTreeFromGenetableURL(url,locus,listOfAllowableNames=None,existingHierarchy=None):
	f = urllib2.urlopen(url)
	html=f.read()
	soup=BeautifulSoup(html)
	clone_names_map=dict()
	table_div=soup.find(id="genetable")
	if table_div is None:
		if(existingHierarchy==None):
			return tree()
		return [existingHierarchy,clone_names_map]
	gene_table=table_div.find('table')
	if gene_table is None:
		if(existingHierarchy==None):
			return tree()
		return [existingHierarchy,clone_names_map]
	rows = gene_table.findAll('tr')
	if rows is None:
		if(existingHierarchy==None):
			return tree()
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
					if(class_val_first_str=="subgroup_middle_note"):
						text = td.find(text=True)# + ';'
						if text is None:
							text="UNDEFINED_SUBGROUP;"
						else:
							text+=";"
						#print "GOT A SUBGROUP=",text
						subgroup=str(text)
						subgroup=removeTerminatingSemicolonIfItExists(subgroup)
					elif(class_val_first_str=="gene_note"):
						text = td.find(text=True)# + ';'
						if text is None:
							text="UNDEFINED_GENE;"
						else:
							text+=";"
						#print "GOT A GENE=",text
						gene_name=str(text)
						gene_name=removeTerminatingSemicolonIfItExists(gene_name)
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
							#SOMETIMES SUBGROUP ISN'T A COLUMN IN THE TABLE SO INHERIT IT FROM THE GENE NAME!
							#if(subgroup==""):
							#	print "Setting a blank subgroup to equal gene_name=",gene_name
							#	subgroup=gene_name
							#elif(not(gene_name.startswith(subgroup))):
							#	print "Setting an old subgroup to equal gene_name=",gene_name
							#	subgroup=gene_name
							if(listOfAllowableNames==None or ((listOfAllowableNames!=None) and (allele_name in listOfAllowableNames))):
								#print "proceeding with setting : gene_name="+gene_name+", allele_name="+allele_name+", and subgroup="+subgroup
								if(subgroup==""):
									my_hier[gene_name][allele_name]
								else:
									my_hier[subgroup][gene_name][allele_name]
							else:
								print "NOTE : skipping the addition of ",allele_name," into the hierarchy cause it's not in the list of allowable alleles!"
	return [my_hier,clone_names_map]
	









#these two defs (tree and dicts) from https://gist.github.com/hrldcpr/2012250
#"One-line Tree in Python"
def tree(): return defaultdict(tree)
def dicts(t): return {k: dicts(t[k]) for k in t}





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
			htmlGeneTablesGlob=base_dir+"/"+organism+"/GeneTables/*"+locus+"*.html"
			geneTableHTMLFiles=glob.glob(htmlGeneTablesGlob)
			print "Got html files"
			printList(geneTableHTMLFiles)
			fastaFilesGlob=base_dir+"/"+organism+"/ReferenceDirectorySet/*"+locus+"*.fna"
			fastaFiles=glob.glob(fastaFilesGlob)
			print "got fasta files :"
			printList(fastaFiles)
			fastaMainMap=dict()
			for fastaFile in fastaFiles:
				fastaString=readFileIntoString(fastaFile)
				fastaList=read_fasta_string(fastaString)
				fastaMap=read_fasta_into_map(fastaList)	
				for fastaKey in fastaMap:
					fastaMainMap[fastaKey]=fastaMap[fastaKey]
			#print "done loading into main map..."
			fastaListOfNames=getIMGTNameListFromFastaMap(fastaMainMap)
			#print "FASTA LIST OF NAMES:"
			fastaListOfNames.sort()
			#printList(fastaListOfNames)
			#print "got name list.... its length is",len(fastaListOfNames)
     			fastaAlleleList=allelifyList(fastaListOfNames)
			#print "got allele list....it is"
			#print "ALLELIFIED LIST:"
			#printList(fastaAlleleList)
			html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[0],locus,fastaAlleleList)
			locusHierarchyData=html_data[0]
			#print "GOT HIERARCHY FROM HD0:"
			#prettyPrintTree(locusHierarchyData)
			#print "DONE SHOWING PRETTY PRINT GOT HIERARCHY FROM HD0"
			clone_names_plain=html_data[1]
			html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[1],locus,fastaAlleleList,locusHierarchyData)
			locusHierarchyData=html_data[0]
			print "GOT HIERARCHY FROM HD0 ORPH:"
			prettyPrintTree(locusHierarchyData)
			print "DONE SHOWING PRETTY PRINT GOT HIERARCHY FROM HD0 ORPH"
			clone_names_orph=html_data[1]	
			clone_names=merge_maps(clone_names_plain,clone_names_orph)
			write_map_to_file(clone_names,geneTableHTMLFiles[0]+".clone_names.map")
			print "GOT FASTA MAP KEY ALLELES :"
			printList(fastaAlleleList)
			print "GOT HIERARCHY:"
			prettyPrintTree(locusHierarchyData)
			print "DONE SHOWING PRETTY PRINT GOT HIERARCHY"
			#print "SHOWING CLONE NAME MAP:"
			#printMap(clone_names)
			print "\n\n\n\n"
			treeAlleles=get_list_of_alleles_appearing_in_tree(locusHierarchyData)
			setSameStats=(set(fastaAlleleList) == set(treeAlleles))
			if setSameStats:
				print "The fasta/RefDirNames ARE the same as the hierarchy allele names!!! :)"
			else:
				print "SAD the fasta/RefDirNames ARE different from the hierarchy allele names!!! :("
			briefSetDiff(fastaAlleleList,treeAlleles,"fasta alleles "+organism+"_"+locus,"tree alleles "+organism+"_"+locus)			
			extendedSortedSetDiff(fastaAlleleList,treeAlleles,"fasta alleles "+organism+"_"+locus,"tree alleles "+organism+"_"+locus)
			print "THIS IS THE PRETTY PRINT FOR LOCUS HIERARCHY DATA, o=",organism,"l=",locus
			prettyPrintTree(locusHierarchyData)
			print "THIS IS THE PRETTY PRINT FOR ORG HIERARCHY DATA"
			org_hierarchy[organism][locus]=locusHierarchyData
			prettyPrintTree(org_hierarchy)
			#from
			#http://stackoverflow.com/questions/82831/how-do-i-check-if-a-file-exists-using-python
			patchPath=geneTableHTMLFiles[0]+".patch"
			if(os.path.isfile(patchPath)):
				print "Found a patch file (",patchPath,") found, so now patching is being performed....."
				locusHierarchyData=patchlocusHierarchyData(locusHierarchyData,patchPath)
				print "AFTER PATCHING, THE HIERARCHY IS NOW :"
				prettyPrintTree(locusHierarchyData)
				print "THIS IS THE PRETTY PRINT FOR PATCHED ORG HIERARCHY DATA"
				org_hierarchy[organism][locus]=locusHierarchyData
				prettyPrintTree(org_hierarchy)
			else:
				print "No patch file (",patchPath,") found, so no patching being performed....."
			clone_names_by_org[organism]=merge_maps(clone_names_by_org[organism],clone_names)
			print "\n\n\n"
	#prettyPrintTree(org_hierarchy)
	pickleFilePath=base_dir+"/hierarchy_data.pkl"
	tbr=[org_hierarchy,clone_names_by_org]
	pickleWrite(pickleFilePath,tbr)
	return tbr



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




def pickleWrite(pPath,o):
	output = open(pPath, 'wb')
	# Pickle dictionary using protocol 0.
	pickle.dump(o, output)
	output.close()




def pickleRead(pPath):
	pkl_file = open(pPath, 'rb')
	data1 = pickle.load(pkl_file)
	pkl_file.close()
	return data1


	
def get_total_tree(tree,this_name,counts_map):
	subtree_sum=0
	if this_name in counts_map:
		subtree_sum+=counts_map[this_name]
	for child in tree:
		subtree_sum+=get_total_tree(tree[child],child,counts_map)
	#print "For subtree rooted at :",this_name,"returning the sum",subtree_sum
	return subtree_sum


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





def jsonify_hierarchy(hier_map,name,counts_map,count_string):
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
	JSON+="\"label\":\""+name+"\",\n"
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

def getFastaListOfDescs(fasta_file):
	fr=open(fasta_file,'r')
	descs=list()
	for line in fr:
		if(line.startswith(">")):
			to_append=line[1:]
			to_append=to_append.strip()
			descs.append(to_append)
	return descs






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


def getADJCDR3EndFromJAllele(jallele,imgtdb_obj,org="human",mode="imgt"):
	if(mode=="imgt"):
		jdescriptor=imgtdb_obj.extractDescriptorLine(jallele,org)
		cdr3_end_raw=getCDR3EndFromJData(jallele,imgtdb_obj,org)
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
			if(start<=cdr3_end_raw and cdr3_end_raw<=stop):
				return cdr3_end_raw-start
			else:
				#out of range
				return (-1)
		else:
			print "FAILED TO MATCH ON INTERVAL REGEX!!!!\n"
			return None
	elif(mode=="kabat"):
		kabat_file="/home/data/DATABASE/01_22_2014/"+org+"/ReferenceDirectorySet/KABAT/Jlookup.tsv"
		reader=open(kabat_file,'r')
		cdr3_end_adj_kabat=(-1)
		for line in reader:
			line=line.strip()
			if(line.startswith(jallele)):
				pieces=line.split("\t")
				cdr3_end_adj_kabat=int(pieces[1])
		reader.close()
		return cdr3_end_adj_kabat
		
		




def getCDR3EndFromJData(allele,imgtdb_obj,org="human"):
	#tmp_file_path="/dev/shm/"+re.sub(r'[^0-9\.a-zA-Z]','',allele)
	#print "tpath is ",tmp_file_path
	#writer=open(tmp_file_path,'w')
	#writer.write(jdata)
	#writer.close()
	#from Bio import SeqIO
	#records=SeqIO.parse(tmp_file_path,"imgt")
	biopyrec=imgtdb_obj.getIMGTDatGivenAllele(allele,True,org)
	records=[biopyrec]
	reg_start=None
	reg_end=None
	for record in records:
		feature_list=record.features
		for feature in feature_list:
			#print "got a feature : ",feature
			ftype=feature.type
			#print "the type is ",ftype	
			qualifiers=feature.qualifiers
			#print "qualifiers : ",qualifiers
			if(ftype=="J-REGION"):
				if("IMGT_allele" in qualifiers):
					#print "found allele in qualifiers!"
					#print "the value of the allele is ",qualifiers["IMGT_allele"]
					allele_qualifier_list=qualifiers["IMGT_allele"]
					actual_value=allele_qualifier_list[0]
					#print "the actual is ",actual_value
					if(actual_value==allele):
						#print "THIS IS THE RIIIIIIIIIIIIIIIIIIIIIGHT ONE!"
						#print "the start is ",feature.location.start
						reg_start=int(re.sub(r'[^0-9]','',str(feature.location.start)))
						#print "fetched is ",reg_start
						reg_end=int(re.sub(r'[^0-9]','',str(feature.location.end)))
	#records=SeqIO.parse(tmp_file_path,"imgt")
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
						#os.remove(tmp_file_path)
						return c_end
	else:
		print "failed to get a start!"
	#os.remove(tmp_file_path)
	return None



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



def isValBetweenTwoNumsInc(n1,n2,val):
	#print "in testing func"
	if(int(n1)<=int(n2)):
		#print "in first isv"
		if(n1<=val and val<=n2):
			return True
	if(int(n2)<=int(n1)):
		#print "in second isv"
		if(n2<=val and val<=n1):
			return True
	#print "n1=",n1,"n2=",n2,"val=",val
	return False


def adjustCDR3StartToSubseq(cdr3_start,descriptor,imgtdb_obj,organism="human"):
	print "in adjustCDR3StartToSubseq"
	print "raw cdr3 start is ",cdr3_start
	desc_pieces=descriptor.split("|")
	print "before rev:"
	print desc_pieces
	reversedFlag=False
	if(desc_pieces[14]=="rev-compl"):
		#desc_pieces[5]=swapIMGTDescInterval(desc_pieces[5])
		#sep="|"
		#descriptor=sep.join(desc_pieces)
		reversedFlag=True
	desc_range=desc_pieces[5]
	interval_re=re.compile(r'(\d+)\.+(\d+)')
	mr=re.search(interval_re,desc_range)
	print "after rev"
	print desc_pieces
	if(mr):
		start=int(mr.group(1))
		stop=int(mr.group(2))
		print "start and stop are ",start," and ",stop
		if(isValBetweenTwoNumsInc(start,stop,cdr3_start)):
			if(reversedFlag):			
				cdr3_adj=cdr3_start-min(start,stop)
				cdr3_adj=max(start,stop)-cdr3_adj
				cdr3_adj-=start
			else:
				cdr3_adj=start-cdr3_start
				cdr3_adj*=(-1)
				cdr3_adj+=1
			print "rev status : ",str(reversedFlag)
			print "Adjusted : ",cdr3_adj
			#sys.exit(0)
			return cdr3_adj
		else:
			#this ref dir seq sequence doesn't cover the CDR3 start! :(
			return (-1)
	else:
		#ref dir set lacks the XX..YY interval markers
		#print "raw cdr3=",cdr3_start
		#print "descriptor=",descriptor
		allele=extractIMGTNameFromKey(descriptor)
		#print "allele\n"
		vdata=imgtdb_obj.getIMGTDatGivenAllele(allele,organism)
		#print "vdata ",vdata
		biopyrec=imgtdb_obj.extractIMGTDatRecordUsingRefDirSetDescriptor(descriptor,True)
		mybiopyseq=str(biopyrec.seq)
		mybiopyseq=mybiopyseq.upper()
		#print "the biopy seq :",mybiopyseq
		refdirsetseq=imgtdb_obj.getRefDirSetFNAGivenCompleteDescriptor(descriptor,organism)
		refdirsetseq=refdirsetseq.upper()
		#print "The ref dir set seq is ",refdirsetseq
		refdirsetseq_noperiods=re.sub(r'[^ACTG]','',refdirsetseq)
		#print "The ref dir set seq (no periods)  is ",refdirsetseq_noperiods
		found_res=mybiopyseq.find(refdirsetseq_noperiods)
		if(found_res==(-1)):
			raise Exception("Error, failed to get interval (start..stop) from "+descriptor," despite using biopython and subseq!")
		else:
			desc_start=found_res
			if(not(desc_start<=cdr3_start and cdr3_start<=desc_end)):
				#not in range!
				return (-1)
			else:
				#print "cdr3_reg=",cdr3_start
				cdr3_adj=desc_start-cdr3_start
				cdr3_adj*=(-1)
				#print "cdr3_adj=",cdr3_adj
				return cdr3_adj	



cdr3_adj_map=dict()
cdr3_adj_map["human"]=dict()
cdr3_adj_map["Mus_musculus"]=dict()
cdr3_adj_map["human"]["imgt"]=dict()
cdr3_adj_map["human"]["kabat"]=dict()
cdr3_adj_map["Mus_musculus"]["kabat"]=dict()
cdr3_adj_map["Mus_musculus"]["imgt"]=dict()

def getAdjustedCDR3StartFromRefDirSetAllele(allele,imgtdb_obj,organism="human",mode="imgt"):
	#use this map as a caching mechanism!
	global cdr3_adj_map
	#if(organism=="blah"):
	#	#print "getAdjustedCDR3StartFromRefDirSetAllele called with a=",allele," and org=",organism
	#	vdata=imgtdb_obj.getIMGTDatGivenAllele(allele,False,organism)
	#	#print "getAdjustedCDR3StartFromRefDirSetAllele vdata is \n",vdata
	#	cdr3_start=getCDR3StartFromVData(vdata,allele,imgtdb_obj,organism)
	#	print "getAdjustedCDR3StartFromRefDirSetAllele raw cdr3_start is ",cdr3_start
	#	if(cdr3_start==(-1)):
	#		#print "getAdjustedCDR3StartFromRefDirSetAllele returning -1 because it is!\n"
	#		return (-1)
	#	else:
	#	
	#		descriptor=imgtdb_obj.extractDescriptorLine(allele,organism)
	#		#print "The extracted descriptor is ",descriptor
	#		cdr3_adjusted=adjustCDR3StartToSubseq(cdr3_start,descriptor,imgtdb_obj,organism)
	#		return cdr3_adjusted
	#elif(organism=="Mus_musculus" ):
	if(1==1):
		if(mode=="imgt"):
			#read the HIGH-VQUEST ANNOTATION
			#print "need to handle IMGT mode for get CDR3 start"
			if(organism in cdr3_adj_map):
				if(allele in cdr3_adj_map[organism]["imgt"]):
					return cdr3_adj_map[organism]["imgt"][allele]
			else:
				cdr3_adj_map[organism]=dict()
			if(organism=="Mus_musculus"):
				baseDir=imgtdb_obj.getBaseDir()+"/Mus_musculus/ReferenceDirectorySet/MOUSE/IMGT_HighV-QUEST_individual_files_folder"
			elif(organism=="human"):
				baseDir=imgtdb_obj.getBaseDir()+"/human/ReferenceDirectorySet/HUMAN_REF/IMGT_HighV-QUEST_individual_files_folder"
			fglobstr=baseDir+"/*"
			imgt_files=glob.glob(fglobstr)
			fileToUse=None
			annotationStartLine=None
			for imgt_file in imgt_files:
				if(fileToUse==None):
					reader=open(imgt_file,'r')
					lineNum=0
					for line in reader:
						sline=line.strip()
						if(sline==">"+allele.strip()):
							fileToUse=imgt_file
						if(sline.startswith("13. Annotation by IMGT/Automat") and not(fileToUse==None)):
							annotationStartLine=lineNum
						lineNum+=1
			if(fileToUse==None):
				#print "Found no file to use"
				pass
			else:
				#print "Found file ",fileToUse,"to use to start at line=",annotationStartLine
				pass
			cdr3_reader=open(fileToUse,'r')
			lineNum=0
			cdr3_re=re.compile(r'^CDR3(\-IMGT)?\s+(\d+)[^0-9]+(\d+)',re.IGNORECASE)
			cdr3_start=None
			for line in cdr3_reader:
				if(lineNum>annotationStartLine):
					sline=line.strip()
					search_result=re.search(cdr3_re,sline)
					if(cdr3_start!=None):
						#print "Looking at ",sline
						pass
					if(search_result):
						#print "FOUND A CDR3 line: ",
						leftInt=int(search_result.group(2))
						rightInt=int(search_result.group(3))
						cdr3_start=leftInt
				lineNum+=1
			if(cdr3_start==None):
				cdr3_start=(-1)
			cdr3_adj_map[organism]["imgt"][allele]=cdr3_start
			return cdr3_adj_map[organism]["imgt"][allele]
		elif(mode=="kabat"):
			#print "need to handle KABAT mode for CDR3 start!"
			if(organism in cdr3_adj_map):
				if(allele in cdr3_adj_map[organism]["kabat"]):
					return cdr3_adj_map[organism]["kabat"][allele]			
			kabat_file=imgtdb_obj.getBaseDir()+"/"+organism+"/ReferenceDirectorySet/KABAT/Vlookup.tsv"
			reader=open(kabat_file,'r')
			for line in reader:
				line=line.strip()
				if(line.startswith(allele)):
					pieces=line.split("\t")
					cdr3_start=pieces[11]
					cdr3_adj_map[organism]["kabat"][allele]=int(cdr3_start)
			reader.close()
			toReturn=None
			if(allele in cdr3_adj_map[organism]["kabat"]):
				toReturn=cdr3_adj_map[organism]["kabat"][allele]
			else:
				cdr3_adj_map[organism]["kabat"][allele]=(-1)
				toReturn=cdr3_adj_map[organism]["kabat"][allele]				
			#print "for kabat (",allele,"), to return ",toReturn
			return toReturn
		else:
			print "NON EXISTENT MODE : ",mode
			sys.exit(0)

		
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
	regions=["FWR1","CDR1","FWR2","CDR2","FWR3","CDR3"]
	if(refOrg in reg_adj_map):
		if(mode in reg_adj_map[refOrg]):
			if(refName in reg_adj_map[refOrg][mode]):
				if(region in reg_adj_map[refOrg][mode][refName]):
					print "using cache for ",refName," o=",refOrg," region=",region," and mode=",mode
					return reg_adj_map[refOrg][mode][refName][region]
			else:
				#refname not in it, so add it
				reg_adj_map[refOrg][mode][refName]=dict()
		else:
			#mode not there! a bad mode!
			print "unknown mode ",mode,"!"
			sys.exit(0)
	else:
		#organism not in there!
		print "unknown organism :",refOrf,"!"
		sys.exit(0)
	#first, form a path to a lookup.
	#this depends on organis and mode
	#for KABAT its a file
	#for IMGT it starts as a directory but later becomes a file!
	lookupFile=None
	if(mode=="kabat"):
		if(refOrg=="human"):
			lookupFile=imgtdb_obj.getBaseDir()+"/human/ReferenceDirectorySet/KABAT/Vlookup.tsv"
		elif(refOrg=="Mus_musculus"):
			lookupFile=imgtdb_obj.getBaseDir()+"/Mus_musculus/ReferenceDirectorySet/KABAT/Vlookup.tsv"
		else:
			print "ERROR, UNKNOWN ORGANISM ",refOrg
			sys.exit(0)
	elif(mode=="imgt" or mode=="IMGT"):
		if(refOrg=="human"):
			lookupFile=imgtdb_objgetBaseDir()+"/human/ReferenceDirectorySet/HUMAN_REF/IMGT_HighV-QUEST_individual_files_folder"
		elif(refOrg=="Mus_musculus"):
			lookupFile=imgtdb_objgetBaseDir()+"/Mus_musculus/ReferenceDirectorySet/MOUSE/IMGT_HighV-QUEST_individual_files_folder"
		else:
			print "ERROR, UNKNOWN ORGANISM ",refOrg
			sys.exit(0)
	else:
		print "ERROR, undefined mode : ",mode
	########################################
	#now that a lookup is selected do an actual lookup!
	print "TO USE LOOKUP : ",lookupFile
	if(mode=="kabat"):
		idx_num=None
		for i in range(len(regions)):
			if(regions[i]==region):
				idx_num=i
		if((idx_num is None) or not(region in regions)):
			print "ERROR, UNKNOWN REGION : ",region
			sys.exit(0)
		kabat_reader=open(lookupFile,'r')
		for line in kabat_reader:
			line=line.strip()
			line_pieces=line.split('\t')
			line_segment=line_pieces[0]
			print "looking at ",line
			if(line_segment==refName):
				col_num=1+idx_num*2
				region_interval=[int(line_pieces[col_num]),int(line_pieces[col_num+1])]
				reg_adj_map[refOrg]["kabat"][refName][region]=region_interval
				kabat_reader.close()
				return reg_adj_map[refOrg]["kabat"][refName][region]
				#return reg_adj_map[refOrg]["kabat"][refName]
		kabat_reader.close()
		print "ERROR, FAILED TO FIND KABAT REGION FOR REFERENCE NAMED "+refName+" in "+lookupFile
		sys.exit(0)
	print "got to bad end...."
			
			
		


def getCDR3StartFromVData(vdata,allele,imgtdb_obj,organism):
	pyobj=imgtdb_obj.getIMGTDatGivenAllele(allele,True,organism)
	descriptor=imgtdb_obj.extractDescriptorLine(allele,organism)
	desc_pieces=descriptor.split("|")
	revCompFlag=False
	if(desc_pieces[14]=="rev-compl"):
		revCompFlag=True
	#find the v region with the allele name
	#and get the range in the v region.
	#find the CDR3 inside the v region.  Return it
	record=pyobj
	vregion_start=None
	vregion_stop=None
	#first find the region start, stop with the matching allele
	for feature_list in record.features:
		#print "got feature list : ",feature_list
		ftype=feature_list.type
		qualifiers=feature_list.qualifiers
		if(ftype=="V-REGION"):
			if("IMGT_allele" in qualifiers):
				allele_val=str(qualifiers["IMGT_allele"][0])
				#print "comparison allele : ",allele
				#print "allele_val is ",allele_val
				if(allele_val==allele):
					#print "\n\n*****MATCH!*****\n\n"
					vregion_start=int(re.sub(r'[^0-9]','',str(feature_list.location.start)))
					#print "fetched is ",vregion_start
					vregion_stop=int(re.sub(r'[^0-9]','',str(feature_list.location.end)))
					#sys.exit(0)
				else:
					pass
					#print "NO MATCH!"
					#sys.exit(0)
	
	if(vregion_start==None or vregion_stop==None):
		#print "found an invalid region....!"
		return (-1)
	else:
		#print "found a region start and stop, so look for cdr3 in it...."
		for feature_list in record.features:
			#print "got a feature list"
			#print feature_list
			ftype=feature_list.type
			if(re.search(r'CDR3',ftype,re.IGNORECASE)):
				print "found a CDR3!"
				cdr3_start=int(re.sub(r'[^0-9]','',str(feature_list.location.start)))
				cdr3_stop=int(re.sub(r'[^0-9]','',str(feature_list.location.end)))
				qualifiers=feature_list.qualifiers
				#print feature_list
				#print qualifiers
				if("complement" in qualifiers):
					revCompFlag=True
				else:
					revCompFlag=False
				if(vregion_start<=cdr3_start and cdr3_start<=vregion_stop):
					#print "returning ",cdr3_start
					#return cdr3_start
					#sys.exit(0)
					if(not(revCompFlag)):
						return cdr3_start
					else:
						return cdr3_stop
	return (-1)

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
		JSON=jsonify_hierarchy(hierarchy,"human",self.count_map,"value")
		return JSON



def getQueryIndexGivenSubjectIndexAndAlignment(query_aln,subject_aln,q_start,q_stop,s_start,s_stop,subject_pos):
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
					#print "getQueryIndexGivenSubjectIndexAndAlignment to return (qcounter) ",q_counter
					return q_counter
			if(not(q_val=="-")):
				q_counter+=1
			if(not(s_val=="-")):
				s_counter+=1
			delta+=1
			#print "\n\n"



def looksLikeAlleleStr(a):
	are=re.compile(r'\*\d+$')
	res=re.search(are,a)
	if(res):
		return True
	else:
		return False



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
	import sys
	refOrg="human"
	imgtdb_obj=imgt_db("/home/data/DATABASE/01_22_2014/")
	print "MY BASE : ",imgtdb_obj.getBaseDir()
	mode="kabat"
	refName="IGHV1-69*13"
	for i in range(3):
		regions=["FWR1","CDR1","FWR2","CDR2","FWR3","CDR3"]
		for region in regions:
			result=getVRegionStartAndStopGivenRefData(refName,refOrg,imgtdb_obj,region,mode)
			print "THE ACTUAL RETURN IS ",result
	#print "***JSON",JSON	
	#s_from=10
	#s_to=20
	#q_from=30
	#q_to=43
	#s_pos=15
	#       0 2  45 67 89"
	#S_ALN="-GATT-AC-AA-TT-"
	#Q_ALN="TGA-TTACACCGTTC"
	#      0   3 5    0 
	# 10 31
	# 11 32
	# 12 33
	# 13 33
	# 14 35
	# 15 36
	# 16 38
	# 17 39
	# 18 41
	# 19 42

	#for sub_desired in range(20-10):
	#	sub_desired+=s_from
	#	print "DESIRED BEFORE PASS IN : ",sub_desired
	#	res=getQueryIndexGivenSubjectIndexAndAlignment(Q_ALN,S_ALN,q_from,q_to,s_from,s_to,sub_desired)
	#	print "***************result (for ",sub_desired,") : ",res
	#	print "\n\n"
	sys.exit(0)
		



	#allele="IGHV4-4*01"
	#allele="IGHJ5*02"
	imgtdb_obj=imgt_db("/home/data/DATABASE/01_22_2014/")
	#print "indexing...\n"
	#imgtdb_obj.indexIMGTDatFile()
	#print "done indexing...."
	#sys.exit(0)
	#data_rec=imgtdb_obj.getIMGTDatGivenAllele(allele)
	#print "the rec is ",data_rec
	#s_start=getCDR3StartFromVDataa(data_rec)
	#print "the s_start is ",s_start
	#jdata=data_rec
	#cdr3_end=getCDR3EndFromJData(jdata,allele)
	#print "got end=",cdr3_end
	#from Bio import SeqIO
	#records=SeqIO.parse("/dev/shm/X86355.embl","imgt")
	#for record in records:
	#	print record
	#	feature_list=record.features
	#	for feature in feature_list:
	#		print "got a feature : ",feature
	#		ftype=feature.type
	#		print "the type is ",ftype
	#dce=getCDR3EndFromJData(data_rec,allele,1,2)
	#print "desired c_end=",dce
	#print "ANALYSIS RESULTS : \n\n\n"
	#analyze_download_dir_forVDJserver("/home/data/DATABASE/01_22_2014/")
	#heavy
	#acc_test="X92266"
	#all_test="IGHV4-55*05"
	#my_res=vamd=imgtdb_obj.getIMGTDatGivenAllele(all_test,"human")
	#print "\n\n"
	#print "the res : "
	#print my_res
	#sys.exit(0)
	vlist="IGHV","IGKV","IGLV"
	jlist="IGHJ","IGKJ","IGLJ"
	organisms=["human","Mus_musculus"]
	#organisms=["human"]
	#organisms=["Mus_musculus"]
	flag=False
	accession_noCDR3=dict()
	allele_miss=dict()
		
	for organism in organisms:
		accession_noCDR3[organism]=dict()
		allele_miss[organism]=0
		if(flag):
			break;
		for p in range(len(vlist)):
			if(flag):
				break;
			vlocus=vlist[p]
			jlocus=jlist[p]
			print "working on ",vlocus," for organism ",organism
			fna_v_glob="/home/data/DATABASE/01_22_2014/"+organism+"/ReferenceDirectorySet/"+vlocus+"*html.fna"
			fna_j_glob="/home/data/DATABASE/01_22_2014/"+organism+"/ReferenceDirectorySet/"+jlocus+"*html.fna"
			fna_v_files=glob.glob(fna_v_glob)
			fna_j_files=glob.glob(fna_j_glob)
			if(len(fna_v_files)!=1 or len(fna_j_files)!=1):
				print "BAD"
				flag=True
				break;
			vmap=read_fasta_file_into_map(fna_v_files[0])
			jmap=read_fasta_file_into_map(fna_j_files[0])
			for vdesc in vmap:
				print "Now have ",vdesc
				v_allele=extractIMGTNameFromKey(vdesc)
				print "got allele ",v_allele
				print "working with ",v_allele," for ",organism
				cdr3_start_adj=getAdjustedCDR3StartFromRefDirSetAllele(v_allele,imgtdb_obj,organism)
				if(cdr3_start_adj!=(-1)):
					print "For rds allele=",v_allele," with vdesc=",vdesc," got adjusted cdr3start=",cdr3_start_adj
					coding_seq=re.sub(r'\.','',vmap[vdesc])[cdr3_start_adj-1:]
					#print "Coding seq (preX) =",coding_seq
					#coding_seq=add_necessary_trailing_N_for_mul3(coding_seq)
					#coding_seq=re.sub(r'N','X',coding_seq)
					print "Coding seq (post) =",coding_seq
					translation=biopythonTranslate(coding_seq)
					print "FOUND Translation=",translation
				else:
					print "Adjusted CDR3 start is (-1)!"
					print "FOUND Translation=NONE"
				print "\n\n\n\n\n\n\n\n\n\n\n"

	#print "***********************************\nEXAMINED ACCESSIONS\n*********************************"
	#for org in organisms:
	#	for va in accession_noCDR3[org]:
	#		print va,"HUM missing CDR3 start !"
	#		allele_to_use=accession_noCDR3[org][va]
	#		print "Using allele=",allele_to_use
	#		print "Using organism=",org
	#		vamd=imgtdb_obj.getIMGTDatGivenAllele(allele_to_use,org)
	#		print "\n",vamd
	#		print "\n\n\n\n"
	#for va in accession_noCDR3['human']:
	#	print va,"HUM missing CDR3 start !"
	#	allele_to_use=accession_noCDR3['human'][va]
	#	print "Using allele=",allele_to_use
	#	vamd=imgtdb_obj.getIMGTDatGivenAllele(allele_to_use,"human")
	#	print "\n\n\n\n\n"
	#for va in accession_noCDR3['Mus_musculus']:
	#	print va,"MUS missing CDR3 start !"
	#	allele_to_use=accession_noCDR3['Mus_musculus'][va]
	#	print "Using allele=",allele_to_use
	#	vamd=imgtdb_obj.getIMGTDatGivenAllele(allele_to_use,"Mus_musculus")
	#	print "\n\n\n\n\n"

























