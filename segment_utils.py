#!/usr/bin/env python

from bs4 import BeautifulSoup
from collections import defaultdict
import pickle
import urllib2
import re
import pprint
from imgt_utils import get_loci_list,imgt_db
import glob
from utils import biopythonTranslate
from igblast_utils import printNiceAlignment
from utils import printMap


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
				return cdr3_end_raw-start+1
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





#given a V segment alignment extract from it the sub-portion for a given
#region.  Return a 2-length array with subject and query alignment data
#return None if alignment too short or doesn't cover region
def getRegionAlignmentFromLargerVAlignment(sub_info_map,org,mode,region_name,imgtdb_obj,wholeOnly=False):
	#print "\n\n\n\nusing region=",region_name
	valid_regions=["FR1","CDR1","FR2","CDR2","FR3"]
	if(not(region_name in valid_regions)):
		print region_name," is an invalid region!!!!"
		return None
	sub_name=sub_info_map['subject ids']
	ref_region_interval=getVRegionStartAndStopGivenRefData(sub_name,org,imgtdb_obj,region_name,mode)
	ref_region_transcript_start=getVRegionStartAndStopGivenRefData(sub_name,org,imgtdb_obj,region_name,"imgt")[0]
	#print "For reference=",sub_name," for org=",org," region=",region_name," got (mode=",mode,")region=",ref_region_interval
	if(ref_region_interval[0]==(-1) and ref_region_interval[1]==(-1)):
		print "Can't get region info, start and stop of ref region is neg 1...."
		return None
	else:
		reg_start=int(ref_region_interval[0])
		if(reg_start==(-1)):
			reg_start=int(sub_info_map['s. start'])
		reg_end=int(ref_region_interval[1])
		if(reg_end==(-1)):
			reg_end=int(sub_info_map['s. end'])
		s_start=int(sub_info_map['s. start'])
		s_end=int(sub_info_map['s. end'])
		q_start=int(sub_info_map['q. start'])
		q_end=int(sub_info_map['q. end'])
		#print "s_start=",s_start," and s_end=",s_end
		s_aln=sub_info_map['subject seq']
		q_aln=sub_info_map['query seq']
		region_alignment=["",""]
		temp_index=0
		temp_index_sbjct=s_start
		temp_index_qury=q_start
		frame_mask=list()
		r_q_start=(-1)
		r_q_end=(-1)
		refName=sub_info_map['subject ids']
		while(temp_index<len(s_aln)):
			if(reg_start<=temp_index_sbjct and temp_index_sbjct<=reg_end):
				#in region
				#subject at 0, query at 1
				if(r_q_start==(-1) and q_aln[temp_index]!="-"):
					r_q_start=temp_index_qury
				if(q_aln[temp_index]!=(-1)):
					r_q_end=temp_index_qury
				region_alignment[0]+=s_aln[temp_index]
				region_alignment[1]+=q_aln[temp_index]
				if(s_aln[temp_index]!="-"):
					#if the frame is knowable, use it
					frame_mask.append(getTheFrameForThisReferenceAtThisPosition(refName,org,imgtdb_obj,temp_index_sbjct))
				else:
					#if the frame is not knowable, put 1 everywhere for the frame to trigger plain sub cts
					frame_mask.append(1)					
			else:
				#not in region
				pass
			if(s_aln[temp_index]!="-"):
				temp_index_sbjct+=1
			if(q_aln[temp_index]!="-"):
				temp_index_qury+=1
			temp_index+=1
		return [region_alignment,frame_mask,r_q_start,r_q_end]



#from a codon list, get an amino list
def getAminosFromCodonSpace(codon_space):
	amino_list=list()
	for codon in codon_space:
		amino_list.append(str(biopythonTranslate(codon)))
	return amino_list


#given a codon, if there are gaps in it
#return all possible codons that are the 
#same but have all 4 nucleotides at the gaps
def getCodonSpace(codon):
	if(codon.find("-")==(-1)):
		return [codon]
	codon_space=list()
	bs=['A','C','G','T']
	for first in bs:
		for second in bs:
			for third in bs:
				temp_codon=first+second+third
				if(isInCodonSpace(codon,temp_codon)):
					codon_space.append(temp_codon)
	return codon_space
	

#it's in the space if it only differs at gap positions
def isInCodonSpace(codon,proposal):
	non_gap_differ=0
	gap_differ=0
	for pos in range(len(codon)):
		if(codon[pos]=="-"):
			gap_differ+=1
		elif(codon[pos]!=proposal[pos]):
			non_gap_differ+=1
	if(non_gap_differ==0):
		return True
	else:
		return False



#change a list of codons to a set
def getCodonSpaceAsSet(codon_space):
	css=set()
	for codon in codon_space:
		css.add(str(codon))
	return css



#combine two maps (assuming same keysets!)
#by adding the values
#c[x]=a[x]+b[x]
def combineCharMaps(map_1,map_2):
	combo_map=dict()
	for k in map_1:
		combo_map[k]=map_1[k]+map_2[k]
	return combo_map



#characterize CDR3
def getCDR3RegionSpecificCharacterization(vData,dData,jData,organism,imgtdb_obj,dMode):
	print "###############################ENTER CDR3#################################"
	char_map=dict()
	num_ins=0
	num_del=0
	num_syn=0
	num_nsy=0
	num_bsb=0
	VrefName=vData['subject ids']
	############################################
	#get CDR3 from V characterization
	v_ref_cdr3_start=getAdjustedCDR3StartFromRefDirSetAllele(VrefName,imgtdb_obj,organism,dMode)
	cdr3_v_char_map=getEmptyRegCharMap()
	qry_cdr3_start=(-1)
	qry_cdr3_start_last_frame=(-1)
	lastVPos=(-1)
	v_tot=0
	v_same=0
	cdr3_q_start_for_return=(-1)
	if(vData!=None and v_ref_cdr3_start!=(-1)):
		v_s_aln=vData['subject seq']
		v_q_aln=vData['query seq']
		vRefTo=int(vData['s. end'])
		vRefFrom=int(vData['s. start'])
		print "V call is ",VrefName
		qry_cdr3_start=getQueryIndexGivenSubjectIndexAndAlignment(v_q_aln,v_s_aln,int(vData['q. start']),int(vData['q. end']),vRefFrom,vRefTo,v_ref_cdr3_start)
		print "The CDR3 start (mode=",dMode,") is ",v_ref_cdr3_start
		temp_v=0
		temp_v_pos=int(vData['s. start'])
		temp_v_q_pos=int(vData['q. start'])
		print "sub start and end are ",temp_v_pos," and ",vRefTo
		cdr3_s_aln=""
		cdr3_q_aln=""
		frame_mask=list()
		while(temp_v<len(v_s_aln)):
			if(temp_v_pos>=v_ref_cdr3_start):
				
				cdr3_s_aln+=v_s_aln[temp_v]
				cdr3_q_aln+=v_q_aln[temp_v]
				qry_cdr3_start_last_frame=getTheFrameForThisReferenceAtThisPosition(VrefName,organism,imgtdb_obj,temp_v_pos)
				if(v_s_aln[temp_v]!="-" and v_q_aln[temp_v]!="-"):
					v_tot+=1
					if(v_s_aln[temp_v]==v_q_aln[temp_v]):
						v_same+=1
				frame_mask.append(qry_cdr3_start_last_frame)
			if(v_q_aln[temp_v]!="-" and cdr3_q_start_for_return!=(-1)):
				cdr3_q_start_for_return=temp_v_q_pos
			if(v_q_aln[temp_v]!="-"):
				temp_v_q_pos+=1
			if(v_s_aln[temp_v]!="-"):
				lastVPos=temp_v_pos
				temp_v_pos+=1
			temp_v+=1
		cdr3_v_char_map=getRegionSpecifcCharacterization(cdr3_s_aln,cdr3_q_aln,"CDR3",frame_mask,dMode)
	if(cdr3_q_start_for_return==(-1)):
		cdr3_q_start_for_return=int(vData['q. end'])+1
	print "THE V CDR3 CHAR MAP is "
	printMap(cdr3_v_char_map)
	##########################################	
	#D CDR3
	cdr3_d_char_map=getEmptyRegCharMap()
	d_tot=0
	d_same=0
	if(dData!=None):
		d_s_aln=dData['subject seq']
		print "D CALL is ",dData['subject ids']
		d_q_aln=dData['query seq']
		#non_frame=[1]*len(d_s_aln)
		d_frame=[]
		q_d_start=int(dData['q. start'])
		temp_q_d_pos=q_d_start
		temp_d=0
		while(temp_d<len(d_s_aln)):
			if(d_s_aln[temp_d]!="-" and d_q_aln[temp_d]!="-"):
				d_tot+=1
				if(d_s_aln[temp_d]==d_q_aln[temp_d]):
					d_same+=1
			if lastVPos!=(-1):
				d_frame.append((temp_q_d_pos-qry_cdr3_start_last_frame)%3)
			else:
				d_frame.append(1)
			if(d_q_aln[temp_d]!="-"):
				temp_q_d_pos+=1
			temp_d+=1
		print "Some d info:"
		cdr3_d_char_map=getRegionSpecifcCharacterization(d_s_aln,d_q_aln,"CDR3",d_frame,dMode)
		print "the d alignment"
		print printNiceAlignment(d_q_aln,d_s_aln)
		print d_frame
		print "d map "
		printMap(cdr3_d_char_map)
	###############################################
	#J CDR3
	cdr3_j_char_map=getEmptyRegCharMap()
	firstJFrame=(-1)
	j_tot=0
	j_same=0
	j_cdr3_reg_end=(-1)
	if(jData!=None):
		jDataRefName=jData['subject ids']
		print "J CALL is ",jDataRefName
		ref_cdr3_end=getADJCDR3EndFromJAllele(jDataRefName,imgtdb_obj,organism,dMode)
		if(ref_cdr3_end!=(-1)):
			j_s_aln=jData['subject seq']
			j_q_aln=jData['query seq']
			j_from=int(jData['s. start'])
			j_to=int(jData['s. end'])
			i_pos=0
			sub_pos=j_from
			qry_pos=int(jData['q. start'])
			cdr3_s_aln=""
			cdr3_q_aln=""
			j_frame=[]
			while(i_pos<len(j_s_aln)):
				if(j_q_aln[i_pos]!="-" and firstJFrame==(-1)):
					firstJFrame=(ref_cdr3_end+1-sub_pos)%3
				if(sub_pos<=ref_cdr3_end):
					j_frame.append((ref_cdr3_end+1-sub_pos)%3)
					cdr3_s_aln+=j_s_aln[i_pos]
					cdr3_q_aln+=j_q_aln[i_pos]
					if(j_q_aln[i_pos]!="-" and j_s_aln[i_pos]!="-"):
						j_tot+=1
						if(j_q_aln[i_pos]==j_s_aln[i_pos]):
							j_same+=1
					if(j_q_aln[i_pos]!="-"):
						j_cdr3_reg_end=qry_pos
				if(j_s_aln[i_pos]!="-"):
					sub_pos+=1
				if(j_q_aln[i_pos]!="-"):
					qry_pos+=1
				i_pos+=1
			print "jcdr3 end is ",ref_cdr3_end," but s_from and to are ",j_from," and ",j_to
			cdr3_j_char_map=getRegionSpecifcCharacterization(cdr3_s_aln,cdr3_q_aln,"CDR3",j_frame,dMode)
			print "THE J CDR3 CHAR MAP and frame info:"
			printMap(cdr3_j_char_map)
			print j_frame
		if(j_cdr3_reg_end==(-1)):
			j_cdr3_reg_end=int(jData['q. start'])-1
	if(j_cdr3_reg_end==(-1)):
		j_cdr3_reg_end=cdr3_q_start_for_return
	############################################
	#cumulation
	combo_map=combineCharMaps(cdr3_v_char_map,cdr3_d_char_map)
	combo_map=combineCharMaps(combo_map,cdr3_j_char_map)
	cum_tot=v_tot+d_tot+j_tot
	print "v,d,j tots",v_tot,d_tot,j_tot
	print "v,d,j smes",v_same,d_same,j_same
	cum_sme=v_same+d_same+j_same
	print "cum same and tot : ",cum_sme," ",cum_tot
	cum_pct_id=float(cum_sme)/float(cum_tot)
	combo_map['pct_id']=cum_pct_id*100.0
	print "THE COMBO MAP is "
	printMap(combo_map)
	#temp_j=0
	#temp_j_pos=int(jData['s. start'])
	#j_s_aln=vData['subject seq']
	#j_q_aln=vData['query seq']
	#cdr3_s_aln=""
	#cdr3_q_aln=""	
	#while(temp_j<len(j_s_aln)):
	#	if(temp_j_pos<j_ref_cdr3_end):
	#		cdr3_s_aln+=j_s_aln[temp_j]
	#		cdr3_q_aln+=j_q_aln[temp_j]
	#	if(j_s_aln[temp_j]!="-"):
	#		temp_j_pos+=1
	#Jrefname=jData['subject ids']
	#jRefTo=int(jData['s. end'])
	#jRefFrom=int(jData['s. start'])
	#j_ref_cdr3_end=getADJCDR3EndFromJAllele(Jrefname,imgtdb_obj,organism,dMode)
	return [combo_map,cdr3_q_start_for_return,j_cdr3_reg_end]
	
	



#given a sub alignment with a region, compute :
#A) synonymous mutations, B) non-synonymous mutations
#C) insertions, D) deletions, E) number stop codons
#F) mutation(sum A-D)
#NOTE "A" and "B" are 'base substitutions'
def getRegionSpecifcCharacterization(s_aln,q_aln,reg_name,frame_mask,mode):
	char_map=dict()
	num_ins=0
	num_del=0
	num_syn=0
	num_nsy=0
	num_bsb=0
	print "\n\nNice alignment (region=",reg_name,") query top, subject bottom : \n",printNiceAlignment(q_aln,s_aln),"\n"
	if(len(s_aln)!=len(q_aln)):
		raise Exception("ERROR, Q_ALN LENGTH NOT EQUAL TO S_ALN LENGTH!?!?!")
	#do counts independent of codons/translations
	for i in range(len(s_aln)):
		if(s_aln[i]=="-"):
			num_ins+=1
		elif(q_aln[i]=="-"):
			num_del+=1
	#do counts with codon information
	temp_index=0
	#print "Now doing codon-based counting..."
	while(temp_index<len(s_aln)):
		#print "temp_index=",temp_index," and frame mask is ",frame_mask[temp_index]
		if(
			temp_index<len(s_aln)-2 and 
			frame_mask[temp_index]==0 and 
			frame_mask[temp_index+1]==1 and 
			frame_mask[temp_index+2]==2
			):
			#print "encountered a codon...."
			s_codon=s_aln[temp_index:(temp_index+3)]
			s_codon=re.sub(r'\-','N',s_codon)
			#print "s codon is ",s_codon
			q_codon=q_aln[temp_index:(temp_index+3)]
			q_codon=re.sub(r'\-','N',q_codon)
			#print "q codon is ",q_codon
			if(s_codon.find("-")==(-1) and q_codon.find("-")==(-1)):
				#no gaps, perform analysis
				s_amino=str(biopythonTranslate(s_codon))
				#print "subject amino :",s_amino
				q_amino=str(biopythonTranslate(q_codon))
				#print "query amino ",q_amino
				for cp in range(3):
					if(s_codon[cp]!=q_codon[cp] and s_amino==q_amino):
						num_syn+=1
						num_bsb+=1
					elif(s_codon[cp]!=q_codon[cp] and s_amino!=q_amino):
						num_nsy+=1
						num_bsb+=1
			else:
				#for codons with gaps
				for cp in range(3):
					if(s_codon[cp]!="-" and q_codon[cp]!="-" and s_codon[cp]!=q_codon[cp]):
						num_bsb+=1
			#q_codon_space=getCodonSpace(q_codon)
			#s_codon_space=getCodonSpace(s_codon)
			#q_codon_set=getCodonSpaceAsSet(q_codon_space)
			#s_codon_set=getCodonSpaceAsSet(s_codon_space)
			#print "The codon space from query codon ",q_codon," is ",q_codon_space," and the set is ",q_codon_set
			#print "The subject cd space from subject ",s_codon," is ",s_codon_space," and the set is ",s_codon_set
			# "The query amino space is ",getAminosFromCodonSpace(s_codon_space)
			temp_index+=3
		else:
			#print "encountered a single...."
			if(s_aln[temp_index]!=q_aln[temp_index] and s_aln[temp_index]!="-" and q_aln[temp_index]!="-"):
				num_bsb+=1
			temp_index+=1
	#have this loop below for pct id
	num_same=0
	num_tot=0
	temp_index=0
	while(temp_index<len(s_aln)):
		num_tot+=0
		if(s_aln[temp_index]==q_aln[temp_index] and s_aln[temp_index]!="-" and q_aln[temp_index]!="-"):
			num_same+=1
		if(s_aln[temp_index]!="-" and q_aln[temp_index]!="-"):
			num_tot+=1
		temp_index+=1
	pct_id=0
	if(num_tot==0):
		pct_id=float(0)
	else:
		pct_id=float(num_same)/float(num_tot)
	#print "done with codon-based counting...."
	char_map['insertions']=num_ins
	char_map['deletions']=num_del
	char_map['base_sub']=num_bsb
	char_map['synonymous_bsb']=num_syn
	char_map['nonsynonymous_bsb']=num_nsy
	char_map['mutations']=num_nsy+num_syn+num_bsb+num_del+num_ins
	char_map['pct_id']=pct_id*100.0
	return char_map



	
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
	return char_map





#at the indicated position find the frame by using the first imgt delineation start as a "base" or "frame of reference"
def getTheFrameForThisReferenceAtThisPosition(refName,organism,imgtdb_obj,refPos):
	regions=["FR1","CDR1","FR2","CDR2","FR3"]
	for region in regions:
		interval=getVRegionStartAndStopGivenRefData(refName,organism,imgtdb_obj,region,"imgt")
		if(not(interval is None)):
			start=int(interval[0])
			if(start!=(-1)):
				#assume start is in frame 0
				return (refPos-start)%3
	raise Exception("ERROR, COULD NOT FIND ANY FRAME POSITION AT ALL FOR ",refName," with organism=",organism)




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
	regions=["FR1","CDR1","FR2","CDR2","FR3","CDR3"]
	if(refOrg in reg_adj_map):
		if(mode in reg_adj_map[refOrg]):
			if(refName in reg_adj_map[refOrg][mode]):
				if(region in reg_adj_map[refOrg][mode][refName]):
					#print "using cache for ",refName," o=",refOrg," region=",region," and mode=",mode
					return reg_adj_map[refOrg][mode][refName][region]
			else:
				#refname not in it, so add it
				#print "adding refname key=",refName
				reg_adj_map[refOrg][mode][refName]=dict()
		else:
			#mode not there! a bad mode!
			print "unknown mode ",mode," (expect kabat or imgt)!"
			sys.exit(0)
	else:
		#organism not in there!
		#print "unknown organism :",refOrf,"!"
		sys.exit(0)
	#first, form a path to a lookup.
	#this depends on organis and mode
	#for KABAT its a file
	#for IMGT it st(<?arts as a directory but later becomes a file!
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
			lookupFile=imgtdb_obj.getBaseDir()+"/human/ReferenceDirectorySet/HUMAN_REF/IMGT_HighV-QUEST_individual_files_folder"
		elif(refOrg=="Mus_musculus"):
			lookupFile=imgtdb_obj.getBaseDir()+"/Mus_musculus/ReferenceDirectorySet/MOUSE/IMGT_HighV-QUEST_individual_files_folder"
		else:
			print "ERROR, UNKNOWN ORGANISM ",refOrg
			sys.exit(0)
	else:
		print "ERROR, undefined mode : ",mode
	########################################
	#now that a lookup is selected do an actual lookup!
	#print "TO USE LOOKUP : ",lookupFile
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
			#print "looking at ",line
			if(line_segment==refName):
				col_num=1+idx_num*2
				region_interval=[int(line_pieces[col_num]),int(line_pieces[col_num+1])]
				reg_adj_map[refOrg]["kabat"][refName][region]=region_interval
				kabat_reader.close()
				return reg_adj_map[refOrg]["kabat"][refName][region]
		kabat_reader.close()
		print "ERROR, FAILED TO FIND KABAT REGION FOR REFERENCE NAMED "+refName+" in "+lookupFile
		#return [(-1),(-1)]
		sys.exit(0)
	elif(mode=="IMGT" or mode=="imgt"):
		lookupBase=lookupFile
		filesToScan=lookupBase+"/*"
		filesToScan=glob.glob(filesToScan)
		for current_file in filesToScan:
			#print "Now to scan in "+current_file
			current_file_refname=None
			target_file_flag=False
			passed_annotation_flag=False
			imgt_reader=open(current_file,'r')
			for line in imgt_reader:
				line=line.strip()
				if(line.startswith(">"+refName.strip())):
					target_file_flag=True
				if(line.startswith("13. Annotation")):
					passed_annotation_flag=True
				if(passed_annotation_flag and target_file_flag and (line.startswith(region+"-IMGT"))):
					#print "Found target line "+line+" in file "+current_file
					regRE=re.compile(r'(\-?IMGT)?\s+(<?)(\d+)[^0-9]+(\d+)(>?)\s*',re.IGNORECASE)
					#note that this regex PROHIBITS the <,> signs
					matchRes=re.search(regRE,line)
					if(matchRes):
						lt=matchRes.group(2)
						gt=matchRes.group(5)
						reg_start=matchRes.group(3)
						if(lt=="<"):
							reg_start=(-1)
						reg_end=matchRes.group(4)
						if(gt==">"):
							reg_end=(-1)
						region_interval=[int(reg_start),int(reg_end)]
						reg_adj_map[refOrg]["imgt"][refName][region]=region_interval
						imgt_reader.close()
						return reg_adj_map[refOrg]["imgt"][refName][region]
			imgt_reader.close()
		reg_adj_map[refOrg]["imgt"][refName][region]=[(-1),(-1)]
		return reg_adj_map[refOrg]["imgt"][refName][region]


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
	mode="imgt"
	refName="IGHV1-69*13"
	for i in range(3):
		regions=["FR1","CDR1","FR2","CDR2","FR3","CDR3"]
		for region in regions:
			result=getVRegionStartAndStopGivenRefData(refName,refOrg,imgtdb_obj,region,mode)
			print "THE ACTUAL RETURN IS Allan Jones",result
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

























