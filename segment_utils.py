#!/usr/bin/python

from bs4 import BeautifulSoup
from collections import defaultdict
import pickle
import urllib2
import re
import pprint
from imgt_utils import *

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
							text="NONE;"
						else:
							text+=";"
						#print "GOT A SUBGROUP=",text
						subgroup=str(text)
						subgroup=removeTerminatingSemicolonIfItExists(subgroup)
					elif(class_val_first_str=="gene_note"):
						text = td.find(text=True)# + ';'
						if text is None:
							text="NONE;"
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



def getQueryIndexGivenSubjectIndexAndAlignment(query_aln,subject_aln,q_start,q_stop,s_start,s_stop,subject_pos):
	if(not(s_start<=subject_pos and subject_pos<=s_stop)):
		#invalid range!
		return None
	else:
		s_counter=s_start
		q_counter=q_start
		delta=0
		q_val=query_aln[delta]
		s_val=subject_aln[delta]
		#print query_aln+"\n"+subject_aln+"\n"
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
			if(((delta)+s_counter)>subject_pos):
				#print "got it"
				return delta+q_counter-1
			if(not(q_val=="-")):
				q_counter+=1
			if(not(s_val=="-")):
				s_counter+=1
			delta+=1
			#print "\n\n"




def analyze_download_dir_forVDJserver(base_dir,countsMap=None,specifiedOrganism=None,specifiedLocus=None):
	organisms=getOrganismList()
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



def getCDR3EndFromJData(jdata,allele):#,s_start,s_stop):
	#get line numbers where /IMGT_allel="allele" is found
	#get line numbers where J-TRP or J-PHEFT   J-PHE 
	#FT                       J-TRP 
	lines=jdata.split("\n")
	allele_subseq="/IMGT_allele=\""+allele+"\""
	trypRE=re.compile(r'^FT\s+J\-TRP\s+[<>]?(\d+)[^\d]')
	phenRE=re.compile(r'^FT\s+J\-PHE\s+[<>]?(\d+)[^\d]')
	allele_lines=list()
	regEnd_lines=list()
	for i in range(len(lines)):
		if(a_subseq_of_b(allele_subseq,lines[i])):
			allele_lines.append(i)
		elif(trypRE.match(lines[i])):
			regEnd_lines.append(i)
		elif(phenRE.match(lines[i])):
			regEnd_lines.append(i)
	#return the lowest JTRP or JPHE line number with line number greater than or equal to the allele_numbers 
	data=None
	for j in range(len(regEnd_lines)):
		geqFlag=False
		for a in range(len(allele_lines)):
			if(regEnd_lines[j]>=allele_lines[a]):
				geqFlag=True
		if(geqFlag):
			data=lines[regEnd_lines[j]]
			print "data is ",data
			res=trypRE.search(data)
			if(res):
				return int(res.group(1))
			res=phenRE.search(data)
			if(res):
				return int(res.group(1))
	if(data==None):
		return (-1)

			

def getCDR3StartFromVData(vdata):
	pieces=vdata.split("\n")
	cdr3re=re.compile("^FT\s+CDR3\-IMGT\s+[<>]?(\d+)[^0-9]")
	for i in range(len(pieces)):
		#print "got line #",i," : ",pieces[i]
		cdr3reMatchRes=cdr3re.match(pieces[i])
		if(cdr3reMatchRes):
			#print "MATCH!"
			start=cdr3reMatchRes.group(1)
			return start
		else:
			pass
	return (-1)


if (__name__=="__main__"):
	import sys
	q_aln="-ATTACA-ATTACA"
	s_aln="GATTACAGATTACA"
	s_f=10
	s_t=22
	q_f=30
	q_t=42
	d_sub=13
	res=getQueryIndexGivenSubjectIndexAndAlignment(q_aln,s_aln,q_f,q_t,s_f,s_t,d_sub)
	print "first res is ",res
	#allele="IGHV4-4*01"
	#allele="IGHJ5*02"
	#imgtdb_obj=imgt_db("/home/data/DATABASE/01_22_2014/")
	#data_rec=imgtdb_obj.getIMGTDatGivenAllele(allele)
	#print "the rec is ",data_rec
	#s_start=getCDR3StartFromVDataa(data_rec)
	#print "the s_start is ",s_start
	#jdata=data_rec
	#cdr3_end=getCDR3EndFromJData(jdata,allele)
	#print "got end=",cdr3_end


























