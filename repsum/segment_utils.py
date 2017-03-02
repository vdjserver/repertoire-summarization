#!/usr/bin/env python

from bs4 import BeautifulSoup
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
							if(text is not None):
								clone_names_map[allele_name]=re.sub(r'\s+','',text)
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
							if(listOfAllowableNames==None or ((listOfAllowableNames!=None) and (allele_name in listOfAllowableNames))):
								#print "proceeding with setting : gene_name="+gene_name+", allele_name="+allele_name+", and subgroup="+subgroup
								if(subgroup==""):
									#if there is no subgroup or 
									#store this way
									my_hier [gene_name][allele_name]
								else:
									#otherwise, store this way, using the complete path/resolution
									my_hier[subgroup][gene_name][allele_name]
							
							else:
								print "NOTE : skipping the addition of ",allele_name," into the hierarchy  (url=",url," and locus=",locus,")  cause it's not in the list of alleles from the Reference Fasta!"
	return [my_hier,clone_names_map]
	











#given a base dir for a database
#and any specified organism and/or locus
#scan the reference directory set FASTAs
#and see what entries of them are in the gene tables (or not)
#present a "side-by-side" diff showing any differences thus 
#allowing user to best know how to "patch" the gene 
#tables so that all fasta entries have a place!
#Also, obtain clone name data
#Finally, return all the hierarchy and clone data!
def analyze_download_dir_forVDJserver(base_dir,pickle_file_full_path=None,countsMap=None,specifiedOrganism=None,specifiedLocus=None):
	myDB=imgt_db(base_dir)
	#analyze_IGBlastLookupsVSIMGTDatLookups(myDB)
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
				print "PATCH FILE CONTENTS : "
				print readFileIntoString(patchPath)
				print "\n"
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



def twoListsMatchByLessEucDistThan(l1,l2,m):
	d1=(l1[0]-l2[0])**2
	d2=(l1[1]-l2[1])**2
	dist=(d1+d2)**(0.5)
	if(dist<=m):
		return True
	else:
		return False



def analyze_IGBlastLookupsVSIMGTDatLookups(imgtdb_obj):
	print "Now comparing IgBLAST region lookup information with imgt.dat lookup information...."
	vq_ig=0
	vq_ig_mm=0
	vq_dd=0
	vq_dd_mm=0
	ig_dd=0
	ig_dd_mm=0
	sdm=0
	sdm_mm=0
	organisms=organisms=imgtdb_obj.getOrganismList()
	lom=dict()
	for organism in organisms:
		print "Now anlyzing for organism ",organism
		seq_types=types=["IG","TR"]
		for seq_type in seq_types:
			print "\tNow analyzing for seq type ",seq_type
			v_gene_alleles=imgtdb_obj.getAlleles("V",organism,seq_type)
			for allele in v_gene_alleles:
				print "\t\tNow analyzing for allele = ",allele
				region_list=["FR1","CDR1","FR2","CDR2","FR3"]
				for region in region_list:
					print "\t\t\tNow analyzing for region ",region," for allele=",allele," for organism =",organism
					imgt_data_ss=imgtdb_obj.getRegionStartStopFromIMGTDat(allele,organism,region)
					print "\t\t\t\timgt ss is ",imgt_data_ss
					igblast_data_ss=getVRegionStartAndStopGivenRefData(allele,organism,imgtdb_obj,region,"imgt")
					print "\t\t\t\tigb ss is ",igblast_data_ss
					v_quest_ss=imgtdb_obj.getVQuestRegionInformation(organism,allele,region)
					print "\t\t\t\tVQ is ",v_quest_ss
					if(twoListsMatch(imgt_data_ss,igblast_data_ss)):
						ig_dd+=1
					else:
						ig_dd_mm+=1
					if(twoListsMatch(v_quest_ss,imgt_data_ss)):
						vq_dd+=1
					else:
						vq_dd_mm+=1
					if(twoListsMatch(v_quest_ss,igblast_data_ss)):
						vq_ig+=1
					else:
						if(twoListsMatchByLessEucDistThan(v_quest_ss,igblast_data_ss,2**(0.5))):
							sdm+=1
						else:
							sdm_mm+=1
						vq_ig_mm+=1
						if(allele+"_"+organism in lom):
							lom[allele+"_"+organism].append(region)
						else:
							lom[allele+"_"+organism]=list()
							lom[allele+"_"+organism].append(region)
	print "NUMBER IGB IMGT.DAT MATCH : ",ig_dd," MISMATCH : ",ig_dd_mm
	print "NUMBER VQ IMGT.DAT MATCH : ",vq_dd," MISMATCH : ",vq_dd_mm
	print "NUMBER VQ IGB MATCH : ",vq_ig," MISMATCH : ",vq_ig_mm
	print "NUMBER VQ IGB MISMATCH WITH EUC. DIST<=SQRT(2) ", sdm, ", NOT : ",sdm_mm
	print "\n\n\n"
	print lom
	for allele_o in lom:
		pieces=allele_o.split("_")
		allele=pieces[0]
		organism=pieces[1]
		if(len(pieces)>2):
			for i in range(2,len(pieces)):
				organism=organism+"_"+pieces[i]		
		print "\n\n\n"
		print "Examining allele : ",allele," and organism ",organism
		regions_to_examine=lom[allele_o]
		print "Now looking at regions ",regions_to_examine," for it\n"
		for region in regions_to_examine:
			print "Looking at region :",region
			imgt_data_ss=imgtdb_obj.getRegionStartStopFromIMGTDat(allele,organism,region)
			print "IMGT DAT reports : ",imgt_data_ss
			igblast_data_ss=getVRegionStartAndStopGivenRefData(allele,organism,imgtdb_obj,region,"imgt")
			print "IGB reports : ",igblast_data_ss
			v_quest_ss=imgtdb_obj.getVQuestRegionInformation(organism,allele,region)
			print "V QUEST reports : ",v_quest_ss
	sys.exit(0)
			




def twoListsMatch(ia1,ia2):
	if(type(ia1)!=type(ia2)):
		#can't be matching if the types aren't matching
		return False
	if((ia1 is None) and (ia2 is None)):
		#they're not lists, but they're both None?
		return True
	if(type(ia1)!=list or type(ia2)!=list):
		#they're not lists!
		return False
	if(len(ia1)!=len(ia2)):
		#lengths not equal
		return False
	for n in range(len(ia1)):
		if(ia1[n]!=ia2[n]):
			return False
	#never failed anywhere!!!!
	return True
	



			


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













def determineAllelicOrdering(a1,a2):
	if(a1==a2):
		return 0
	else:
		aa=[a1,a2]
		o_aa=orderAlleicArrayWithZeroPadding(aa)
		if(a1==o_aa[0]):
			return (-1)
		else:
			return 1	



def orderAlleicArrayWithZeroPadding(a_arr):
	#the reason for all this "padded" stuff is so that digits are zero-padded
	#so that ordering (based on ASCII) will put IGHV4 before IGHV40 (and not IGHV40 before IGHV4)
	original_to_padded=dict()
	padded_to_original=dict()	
	for a in a_arr:
		padded_key=zeroPadDigitsTo(hier_map_key)
		original_to_padded[a]=padded_key
		padded_to_original[padded_key]=hier_map_key
		padded_key_list.append(padded_key)
	padded_key_list.sort()
	ready_list=list()
	for p in range(len(padded_key_list)):
		ready_list.append(padded_to_original[padded_key_list[p]])
	return ready_list





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
#the lookup represents the position in the reference data IgBLASTED
#against that is just before the J-TRP or J-PHE (one bp which is in the last codon position of the immediately preceding amino residue)
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
			ellipsis_interval=desc_pieces[5]
			ellipsis_interval=ellipsis_interval.strip()
			if(ellipsis_interval==""):
				accession_num=desc_pieces[0]
				ss_arr=imgtdb_obj.getSegmentRegionStartStopFromIMGTDatGivenAlleleAndAccession(jallele,accession_num,"J")
				desc_pieces[5]=str(ss_arr[0])+".."+str(ss_arr[1])
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
#return the index of the bp immediately before the first bp of the codon of JTRP or JPHE
def getCDR3EndFromJData(allele,imgtdb_obj,org="human"):
	#print "*********************************\nEntering getCDR3EndFromJData\n*******************************"
	biopyrec=imgtdb_obj.getIMGTDatGivenAllele(allele,True,org)
	records=[biopyrec]
	reg_start=None
	reg_end=None
	extracted_descriptor=imgtdb_obj.extractDescriptorLine(allele,org)
	#print "The extracted descriptor from inputs ",allele," and ",org," is (extracted) : "+str(extracted_descriptor)
	extracted_descriptor_pieces=imgtdb_obj.extractIMGTDescriptorPieces(extracted_descriptor)
	#print "The pieces is ",extracted_descriptor_pieces
	extracted_descriptor_interval=imgtdb_obj.getStartStopFromIMGTDescPieces(extracted_descriptor_pieces)
	#print "extracted interval : ",extracted_descriptor_interval
	if(extracted_descriptor_interval==None):
		#print "COULDN'T GET INTERVAL FOR ALLELE=",allele," organism=",org
		#sys.exit(0)
		#print "getCDR3EndFromJData returning 'None' for allele="+str(allele)+" and org="+str(org)
		return None		
	if(extracted_descriptor_interval!=None and type(extracted_descriptor_interval)==list):
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
						#print "found a jtrp with c_start and c_end being ",c_start," and ",c_end
						#print "NOTE c_start already subtrats one...compared to imgt.dat....why???"
						if(reg_start<=c_end and c_end<=reg_end):
							#this is it! it's within the range
							return c_start
		else:
			#print "failed to get a start!"
			pass
	#os.remove(tmp_file_path)
	#print "getCDR3EndFromJData returning 'None' for allele="+str(allele)+" and org="+str(org)
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
		if(mode in cdr3_adj_map[organism]):
			#mode in the dict
			if(allele in cdr3_adj_map[organism][mode]):
				return cdr3_adj_map[organism][mode][allele]
			else:
				#do the actual looking up below!
				pass
		else:
			#mode not in the dict
			cdr3_ajd_map[organism][mode]=dict()
	else:
		cdr3_adj_map[organism]=dict()
		cdr3_adj_map[organism][mode]=dict()
	cdr3_int=getVRegionStartAndStopGivenRefData(allele,organism,imgtdb_obj,"CDR3",mode)
	cdr3_adj_map[organism][mode][allele]=cdr3_int[0]
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
	elif(
		(ref_region_interval[0]!=(-1) and  ref_region_interval[0]>int(sub_info_map['s. end'])) or 
		(ref_region_interval[1]!=(-1) and ref_region_interval[1]<int(sub_info_map['s. start']))
		):
		#this is for the case where the region isn't covered at all!
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



def getTheFrameForThisJReferenceAtThisPosition(refJName,organism,imgtdb_obj,refPos):
	ref_cdr3_end=getADJCDR3EndFromJAllele(refJName,imgtdb_obj,organism,"imgt")
	if(ref_cdr3_end==(-1)):
		return None
	else:
		#ref_cdr3_end in in FRAME #2 because it represents the last BP of the codon immediately before W or F of CDR3 end.
		if(refPos==ref_cdr3_end):
			return 2
		else:
			if(refPos<ref_cdr3_end):
				diff=ref_cdr3_end-refPos
				neg_diff=(-1)*(diff)
				frame=((neg_diff)-1)%3
				return frame
			else:
				diff=refPos-ref_cdr3_end
				frame=(diff-1)%3
				return frame
				
		




#at the indicated position find the frame by using the
#imgt delineation/region starts as a "base" or "frame of reference"
#The region used is the region whose start/stop are *entirely* within the bounds of
#the subject!  Otherwise, an uncovered region is used or a region start "off the read"
#is used which may give incorrect frame
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



class recombFreqManager():

	#have a frequency map
	freq_map=None

	def itemToStr(self,item):
		if(item is None):
			return "None"
		else:
			return str(item)

	#constructor
	def __init__(self):
		self.freq_map=dict()
	
	#increment the count for a particular combination	
	def addVDJRecombination(self,vName,dName,jName,maskAllele=True):
		if(maskAllele):
			#print "v=",vName,"d=",dName,"j=",jName
			vName=deAllelifyName(self.itemToStr(vName))
			dName=deAllelifyName(self.itemToStr(dName))
			jName=deAllelifyName(self.itemToStr(jName))
		comb=[self.itemToStr(vName),self.itemToStr(dName),self.itemToStr(jName)]
		comb_str=",".join(comb)
		if(comb_str in self.freq_map):
			self.freq_map[comb_str]=self.freq_map[comb_str]+1
		else:
			self.freq_map[comb_str]=1

	#write JSON to a file
	def writeJSONToFile(self,outPath):
		writer=open(outPath,'w')
		JSON=self.makeJSON()
		writer.write(JSON)
		writer.close()


	#write data to CSV
	def writeDataToCSV(self,maskAllele=False):
		csv_dict=dict()
		if(maskAllele):
			#merge over alleles into genes
			for triple in self.freq_map:
				vdj=triple.split(',')
				v=deAllelifyName(vdj[0])
				d=deAllelifyName(vdj[1])
				j=deAllelifyName(vdj[2])
				gene_comb=str(v)+","+str(d)+","+str(j)
				if(not(gene_comb in csv_dict)):
					csv_dict[gene_comb]=1
				else:
					csv_dict[gene_comb]+=1
		else:
			csv_dict=self.freq_map
		ret_val="V gene,J gene,D gene\n"
		for triple in sorted(csv_dict):
			ret_val+=triple+","+str(csv_dict[triple])+"\n"
		return ret_val


	#create JSON from combos
	def makeJSON(self):
		json="{\ncombinations : [\n"
		combo_array=list()
		for combination in self.freq_map:
			combo_json="{"
			combo_pieces=combination.split(",")
			combo_json+=" v : '"+combo_pieces[0]+"',"
			combo_json+=" d : '"+combo_pieces[1]+"',"
			combo_json+=" j : '"+combo_pieces[2]+"',"
			combo_json+=" count : "+str(self.freq_map[combination])
			combo_json+="}"
			combo_array.append(combo_json)
		all_combos=",\n".join(combo_array)
		json+=all_combos+"\n]}"
		return json



	
		







#
# This function translates a position from germline coordinates into query coordinates
# given an alignment between the two. It takes into account indels.
#
# It may be the case that the position is beyond the alignment, in which case it
# assumes a straight translation and shifts by remaining bp.
#
#the lean right means that if the subject pos aligns to a gap, then choose the bp position to the right
#otherwise if the subject aligns to a gap at the position, lean "left" and choose the previous bp positoin
#RIGHT is the default!
def getQueryIndexGivenSubjectIndexAndAlignment(query_aln,subject_aln,q_start,q_stop,s_start,s_stop,subject_pos,lean="right"):
	#print "in getQueryIndexGivenSubjectIndexAndAlignment"
	#print "query_aln=",query_aln
	#print "subct_aln=",subject_aln
	#print "q_start=",q_start
	#print "q_stop=",q_stop
	#print "s_start=",s_start
	#print "s_stop=",s_stop
	#print "desired=",subject_pos
        s_counter=s_start
        q_counter=q_start
        delta=0
        q_val=query_aln[delta]
        s_val=subject_aln[delta]

        # if the position is before the alignment even starts
        # then best we can do is shift by the bp diff
        if (subject_pos < s_start):
                return q_counter - (s_start - subject_pos)

        # walk through the alignment
        while(delta<len(query_aln)):
                q_val=query_aln[delta]
                s_val=subject_aln[delta]
                # if we reach the desired position
                if(subject_pos<=s_counter):
                        if(s_val!="-"):
                                if(lean=="right"):
                                        return q_counter
                                elif(lean!="right" and q_val=="-"):
                                        return q_counter-1
                                else:
                                        return q_counter
                # - indicates an indel
                if(not(q_val=="-")):
                        q_counter+=1
                if(not(s_val=="-")):
                        s_counter+=1
                delta+=1
		#print "q_counter=",q_counter,"s_counter=",s_counter,"delta=",delta

        # if we drop out of loop, then likely the desired position is beyond the end
        # of the alignment, so shift by the remaining bp
        return (subject_pos - s_counter) + q_counter

#rooted at an organism get the hierarchy from the gene tables
def getHierarchyBy(geneTablesDirectoryOfHTMLFiles,org_name,filterbyFastaAlleles=False,fullPklPath=None):
	print "Trying to use " ,fullPklPath
	if(not(fullPklPath==None)):
		#print "not none"
		if(os.path.exists(fullPklPath)):
			#print "it exists"
			unpickled_data=pickleRead(fullPklPath)
			#print unpickled_data
			#print "len of data is ",len(unpickled_data)
			subset_for_extraction=unpickled_data[0]
			for item in subset_for_extraction:
				if(item==org_name):
					return subset_for_extraction[item]
				#print item
			#return unpickled_data[org_name]
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



def testAcquireRegionInfo():
	imgtdb_obj=imgt_db("/home/data/DATABASE/10_29_2014/")
	orgs=imgtdb_obj.getOrganismList()
	for organism in orgs:
		#vsegs
		types=["IG","TR"]
		for seq_type in types:
			v_segs=getAlleles("V",organism,seq_type)
			regions="FR1","FR2","FR3","CDR1","CDR2"
			




if (__name__=="__main__"):
	pass
	#getQueryIndexGivenSubjectIndexAndAlignment
	#getQueryIndexGivenSubjectIndexAndAlignment(query_aln,subject_aln,q_start,q_stop,s_start,s_stop,subject_pos):
	#query_aln="AC-XN-GT"
	#sbjct_aln="-GAT-ACA"
	#query_from=2
	#query_to=7
	#sbjct_f=3
	#sbjct_t=8
	#for i in range(3,8+1):
	#	print "\n\n\n"
	#	print "Now analyzing (s at top, q at bottom):"
	#	print sbjct_f,sbjct_aln,sbjct_t
	#	print query_from,query_aln,query_to
	#	print "q_from=",query_from," query_to=",query_to,",s_from=",sbjct_f,", s_to=",sbjct_t," At subject position=",i
	#	#print "The corresponding query position is "+str(getQueryIndexGivenSubjectIndexAndAlignment(query_aln,sbjct_aln,query_from,query_to,sbjct_f,sbjct_t,i))
	#	#print "The corresponding (lean_left) query position is "+str(getQueryIndexGivenSubjectIndexAndAlignment(query_aln,sbjct_aln,query_from,query_to,sbjct_f,sbjct_t,i,"left"))






