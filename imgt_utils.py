#!/usr/bin/env python

from sets import Set
import subprocess
import urllib2
from bs4 import BeautifulSoup
from collections import defaultdict
import pprint
import re
from utils import *
import glob
import json
from subprocess import call
from Bio.Blast import NCBIXML
import sys
import traceback
import pickle
from Bio import SeqIO
import glob


#find the number of leaves
def countNumTerminalEntriesInHierarchy(h):
	num_kids=0
	for k in h:
		if(len(k.keys)==0):
			num_kids+=1
		else:
			num_kids+=countNumTerminalEntriesInHierarchy(h[k])
	return num_kids












#from an IMGT descriptor, extract the IMGT name (eg IGHV4-1*01)
def extractIMGTNameFromKey(k):
	pieces=k.split("|")
	#print "the pieces are :",pieces
	#sys.exit(0)
	return pieces[1]




#get the IMGT URL base
#used in the downloading process
def getIMGTURLBase():
	return "http://www.imgt.org/"





#get the list of LOCI
def get_loci_list():
	loci=["IGHV","IGHD","IGHJ","IGKV","IGKJ","IGLV","IGLJ","TRAV","TRAJ","TRBV","TRBD","TRBJ","TRDV","TRDD","TRDJ","TRGV","TRGJ"]
	return loci





#of the 17 loci, return the ones defined as "heavy"
def get_heavy_loci():
	defined_as_heavy=["IGHD","IGHJ","IGHV","TRBD","TRBJ","TRBV","TRDD","TRDJ","TRDV"]
	return defined_as_heavy





#get loci that are light (not heavy)
def get_light_loci():
	all_loci=get_loci_list()
	heavy_loci=get_heavy_loci()
	light_loci=list()
	for locus in all_loci:
		if(locus in heavy_loci):
			pass
		else:
			light_loci.append(locus)
	return light_loci




#is the given loci one of the 17???
def isLegitimateLoci(locus):
	loci_list=get_loci_list()
	locus=locus.upper()
	return (locus in loci_list)

#given a species and locus
#form a URL to download the FNA data
#that is the REFERENCE DIRECTORY SET!
#note only "allowed" species and loci are used
#in the formation of URLs
#unallowed triggers an exception
def formRefDirURL(species,locus):
	base=getIMGTURLBase()
	allowed=dict()
	allowed['Homo+sapiens']=1
	allowed['Mus_musculus']=1
	allowed['Mus_spretus']=1
	allowed['Rattus+norvegicus']=1
	allowed['Oryctolagus+cuniculus']=1
	allowed['Oncorhynchus+mykiss']=1
	allowed['Macaca+mulatta']=1
	allowed['Danio+rerio']=1
	allowed['Sus+scrofa']=1
	if((species in allowed) and (isLegitimateLoci(locus))):
		#THESE URLS ARE FROM http://www.imgt.org/genedb/html/directlinks.html
		#base+="IMGT_GENE-DB/GENElect?query=7.2+"+locus+"&species="+species
		#http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.2+TRDV&species=Homo+sapiens   
		#IMGT/GENE-DB reference sequences in FASTA format:
		#Nucleotide sequences for F+ORF+all P alleles
		#for F+ORF+in-frame P alleles, including orphons

		#http://www.imgt.org/IMGT_GENE-DB/GENElect?query=7.14+TRDV&species=Homo+sapiens 
		#IMGT/V-QUEST reference sequences in FASTA format:  
		#Nucleotide sequences with gaps according to the IMGT unique numbering 
		#for F+ORF+in-frame P alleles, including orphons
		#see also URLS/links here : http://www.imgt.org/vquest/refseqh.html
		#IMGT/V-QUEST reference directory sets
		#    The IMGT/V-QUEST reference directory sets are constituted by sets of sequences which contain the V-REGION, D-REGION and J-REGION alleles, isolated from the Functional (F), ORF and in-frame pseudogene (P) allele IMGT reference sequences. By definition, these sets contain one sequence for each allele. Allele names of these sequences are shown in red in Alignments of alleles.
		#    The IMGT/V-QUEST reference directory sets also include orphons.
		#    The human and the mouse IG and TR sets are exhaustive and correspond to the IMGT/GENE-DB content in terms of F+ORF+in-frame P genes and alleles.
		#    IMGT/V-QUEST reference directory sets comprise IG and TR sets (for one or several species) which are used in IMGT/V-QUEST for the identification and alignment of the input sequence. 
		base+="IMGT_GENE-DB/GENElect?query=7.14+"+locus+"&species="+species
		return base

		

	else:
		exceptionStatus="Invalid locus or invalid species for reference directory set URL formation!"
		exceptionStatus+="\nLocus ("+locus+") legitimacy : "+str(isLegitimateLoci(locus))
		exceptionStatus+="\nSpecies ("+species+") legitimacy : "+str((species in allowed))
		raise Exception(exceptionStatus)




#given HTML data, count the number
#of appearances of "&gt;", which,
#when decoded, is ">" which is the fasta descriptor signal
def countNumApparentFastaRecsInStr(s):
	if(s is None):
		return (-1)
	lines=s.split("\n")
	n=0
	for l in lines:
		if(l.startswith("&gt;")):
			n+=1
		elif(l.startswith(">")):
			n+=1
		else:
			#print "Couldn't find a fasta in ",l[0:10],"...."
			pass
	return n




#from <pre>...</pre> HTML tag
#extract the fasta data within it
#and return it as a string
def filterOutFastaStringFromFastaPre(p):
	dnaRE=re.compile(r'[actg\.]')
	lines=p.split("\n")
	good_data=list()
	for line in lines:
		if(line.startswith("&gt;") or dnaRE.match(line)):
			good_data.append(line.replace("&gt;",">"))
	nl="\n"
	result=nl.join(good_data)
	#print "Returning the result:",result
	return result





#given a locus and species download the corresponding
#reference directory set data from IMGT and return
#it as a string
#use BEAUTIFUL SOUP for parsing
def downloadRefDirFasta(locus,species,URLOverRide=None):
	if(URLOverRide is None):
		url=formRefDirURL(species,locus)
	else:
		url=URLOverRide
	#url="file:///home/data/vdj_server/igblast_routines/GENElect?query=7.2+IGHD&species=Homo+sapiens"
	#Number of results=162
	print "THE URL IS ",url
	#f = urllib2.urlopen(url)
	#html=f.read()
	#print html
	html= readAURL(url)
	#print "THE HTML READ IS ",html
	matchObj = re.search( r'Number\s+of\s+results\s*=\s*(\d+)', html, re.M|re.I)
	#matchObj = re.match( r'=(\d+)', html, re.M|re.I)
	if(matchObj):
		num_expected=int(matchObj.group(1))
		print "THE NUM EXPECTED IS ",num_expected
		soup=BeautifulSoup(html)
		y=soup.find_all('pre') #returns data between <pre> tags
		pre_num=0
		for a in y:
			print "GOT A PRE ",pre_num
			#print a[0:10]
			pre_num+=1
		print "Total number pre found : ",pre_num
		pre_num=0			
		for a in y:
			print "GOT A PRE ",pre_num
			pre_num+=1
			#print a
			print "the type : ",type(a)
			if a is not None:
				z_num=int(countNumApparentFastaRecsInStr(str(a)))
				print "An apparent is ",z_num
				#print "The source apparent is ",a
				if(z_num==num_expected):
					#print "got expected fasta: ",a
					fastaString=a
					#print fastaString
					fastaString=filterOutFastaStringFromFastaPre(str(fastaString))
					return str(fastaString)
				else:
					#print "Error, expected ",num_expected," FASTA records, but found ",z_num," records instead!"
					pass
			else:
				print "SOUP failed to find a <pre> tag from which FASTA can be extracted!"
		
	else:
		print "\nFAILURE TO MATCH! from URL=",url





#given a species and a locus
#form a URL for subsequent download 
#from IMGT
def formGeneTableURLs(species,locus):
	#http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=human&group=IGHV
	allowed=dict()
	allowed['human']=1
	allowed['Mus_musculus']=1
	if( (species in allowed)  and (isLegitimateLoci(locus))):
		base=getIMGTURLBase()
		base+="IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species="+species+"&group="+locus
		withOrph=base+"&orphon"
		stuff=list()
		stuff.append(base)
		stuff.append(withOrph)
		return stuff	
	else:
		exceptionStatus="Invalid locus or invalid species for gene table set URL formation!"
		exceptionStatus+="\nLocus ("+locus+") legitimacy : "+isLegitimateLoci(locus)
		exceptionStatus+="\nSpecies ("+species+") legitimacy : "+(species in allowed)
		raise Exception(exceptionStatus)




#given a list, return it back
#but with each element having been
#applied with "strip"
def trimList(l):
	for idx, val in enumerate(l):
		#print idx, val
		l[idx]=l[idx].strip()
	return l







#Given a map (parent/child) and a name
#get the lineage
def hier_look(pmap,name,max_iter):
	#here pmap is a map keys are children, values are parents
	#name is a child whose lineage to root is desired
	#max_iter is a maximum number of iterations to use in the recursive lookup
	if(max_iter<=0):
		return "ERR_STACKOVERFLOW_"+name
	if(name in pmap):
		parent=pmap[name]
		if(parent==name):
			return (name)
		else:
			return name+"->"+hier_look(pmap,parent,max_iter-1)

	else:
		print "ERROR, name=",name," not a key!"
		return "ERR"















#from a tree object, return a map of parent/child relationships
#keys are kids
#values are their parents
def getPMapFromTree(t,emap,currentParent):
	print "getPMapFromTree called with currentParent="+currentParent
	for k in t:
		print "setting k/v pair with k=",str(k),"and v=",str(currentParent)
		emap[str(k)]=str(currentParent)
		emap=getPMapFromTree(t[k],emap,k)
	return emap	











#the IMGT class
#object container for data and methods
#for database FNA and GENETABLE files
class imgt_db:
	####################
	#data members
	org_allele_name_desc_map=None
	org_allele_desc_gapless_map=None
	db_base=None
	db_idx_extension=".acc_idx"
	accession_start_stop_map=None
	accession_dat_file_map=None
	imgt_dat_rel_path="www.imgt.org/download/LIGM-DB/imgt.dat"
	imgt_genedb_rel_path="www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP"
	imgt_dat_path=None
	indexPath=None
	ref_dir_set_desc_seqs_map=None
	ol=["human","Mus_musculus"]
	pickle_file_name="hierarchy_data.pkl"

	####################
	#constructor(s)
	def __init__(self,init_db_base):
		self.db_base=init_db_base
		self.imgt_dat_path=self.db_base+"/"+self.imgt_dat_rel_path

	####################
	#function members

	def getPickleFullPath(self):
		return self.getBaseDir()+"/"+self.pickle_file_name



	#for each locus, for each organism, dump it into a FASTA to be blast formatted
	def prepareFASTAForBLASTFormatting(self):
		organisms=self.getOrganismList()
		for organism in organisms:
			loci=get_loci_list()
			for locus in loci:
				rds_base=self.db_base+"/"+organism+"/ReferenceDirectorySet/"
				source_html_glob_str=rds_base+locus+"*fna"
				glob_res=glob.glob(source_html_glob_str)
				if(len(glob_res)!=1):
					no_file_msg="ERROR, failed to find files ",source_html_glob_str," for BLAST formatting preparation!"
					raise Exception(no_file_msg)
				source_html_fna=glob_res[0]
				target_fna=rds_base+"/"+organism+"_"+locus[0:2]+"_"+locus[3]+".fna"
				fna_map=read_fasta_file_into_map(source_html_fna)
				blast_fna_map=dict()
				for desc in fna_map:
					blast_desc=getIMGTNameFromRefDirSetDescriptor(desc)
					blast_fna_map[blast_desc]=fna_map[desc]
				blast_writer=open(target_fna,'a')
				num_descs=len(fna_map)
				rec_num=0
				for blast_desc in blast_fna_map:
					desc_str=">"+blast_desc
					dna_str=blast_fna_map[blast_desc]
					dna_str=re.sub(r'\.','',dna_str)
					blast_writer.write(desc_str+"\n"+dna_str)
					if(rec_num==num_descs-1):
						pass
					else:
						blast_writer.write("\n")
				blast_writer.close()
		#downloading of some databases have shown that some sequences can be duplicated
		#which leads to BLAST formatting errors....so the code below reads a fasta into a
		#map and writes it back....this way effectively removing duplicates
		#esalina2@eddiecomp:/home/data/DATABASE/04_22_2014/Mus_musculus/ReferenceDirectorySet$ grep -P 'TRAV15\-1/DV6\-1\*02' *.fna
		#Mus_musculus_TR_V.fna:>TRAV15-1/DV6-1*02
		#Mus_musculus_TR_V.fna:>TRAV15-1/DV6-1*02
		#TRAV.html.fna:>X63939|TRAV15-1/DV6-1*02|Mus musculus|(F)|V-REGION|16..303|292 nt|1| | | | |292+42=334| | |
		#TRDV.html.fna:>X63939|TRAV15-1/DV6-1*02|Mus musculus|(F)|V-REGION|16..303|292 nt|1| | | | |292+42=334| | |
		for organism in organisms:
			loci=get_loci_list()
			for locus in loci:
				segments=['V','D','J']
				for segment in segments:
					seq_types=["IG","TR"]
					for seq_type in seq_types:
						fasta=self.db_base+"/"+organism+"/ReferenceDirectorySet/"+organism+"_"+seq_type+"_"+segment+".fna"
						fasta_map=read_fasta_file_into_map(fasta)
						fasta_writer=open(fasta,'w')
						for desc in fasta_map:
							fasta_writer.write(">"+desc+"\n"+fasta_map[desc]+"\n")
						fasta_writer.close()

					
				
	#format the BLAST database!
	def blastFormatFNAInRefDirSetDirs(self,makeblastdbbin):
		organisms=self.getOrganismList()
		for organism in organisms:
			rds_base=self.db_base+"/"+organism+"/ReferenceDirectorySet/"
			segments=['V','D','J']
			for segment in segments:
				seq_types=["IG","TR"]
				for seq_type in seq_types:
					fasta=self.db_base+"/"+organism+"/ReferenceDirectorySet/"+organism+"_"+seq_type+"_"+segment+".fna"
					blastFormatFLEX(fasta,makeblastdbbin,True)



	#return the organism list
	def getOrganismList(self,fromHardCode=True):
		if(fromHardCode):
			#get from a hard-coded list in this file
			return self.ol
		else:
			#get downloaded organism!
			thedir=self.db_base
			org_list=list()
			# http://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python
			total_list=[ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
			#print "total_list is ",total_list
			search_and_avoid=["\.py","down","\.pkl","\.sh","\.zip"]
			for tl in total_list:
				tls=str(tl.strip())
				matched_no_regex=True
				for s in search_and_avoid:
					if(re.search(s,tls,re.IGNORECASE)):
						matched_no_regex=False
				if(matched_no_regex):
					org_list.append(tl.strip())
			return org_list



	#return base directory
	def getBaseDir(self):
		return self.db_base

	#return base directory alias
	def getDirBase(self):
		return self.getBaseDir()




	#get the alleles in the database (by segment, organism, and seq_type)
	def getAlleles(self,segment="V",org="human",seq_type="IG"):
		fastaPath=self.getDirBase()+"/"+org+"/ReferenceDirectorySet/"+org+"_"+seq_type+"_"+segment+".fna"
		if(not(os.path.isfile(fastaPath))):
			raise Exception("Error, fasta file "+fastaPath+" does not exists!  Is an invalid organism, seq type or segment used???")
		fastaDict=read_fasta_file_into_map(fastaPath)
		fastaKeys=fastaDict.keys()
		fastaKeysList=list(fastaKeys)
		return fastaKeysList

	
	#get a list of genes, withouth the alleles designation
	def getGenes(self,alleles=None,segment="V",org="human",seq_type="IG"):
		if(alleles==None):
			alleles=self.getAlleles(segment,org,seq_type)
		gene_set=set()
		for allele in alleles:
			if(looksLikeAlleleStr(allele)):
				#good!
				gene_name=deAllelifyName(allele)
				if(not(gene_name in gene_set)):
					gene_set.add(gene_name)
			else:
				#huh?
				raise Exception("Error, the allele "+allele+" doesn't look allelic by name !")
		gene_list=list(gene_set)
		gene_list.sort()
		return gene_list

		


	#download GENE tables only
	def download_GeneTables(self,unconditionalForceReplace=False):
		print "in download_imgt_RefDirSeqs_AndGeneTables_HumanAndMouse"
		base=self.db_base
		organisms=self.getOrganismList(True)
		print "To download for ",organisms
		for organism in organisms:
			#do all organisms
			loci=get_loci_list()
			print "To download for ",loci
			for locus in loci:
				#do all 17 groups
				print "Downloading",locus,"for",organism,"at",formatNiceDateTimeStamp(),"..."
				#first download and save the gene tables
				geneTableOrgName=organism
				GeneTablesURLs=formGeneTableURLs(geneTableOrgName,locus)		
				#download each regular table and orphon table and write to file
				regularURL=GeneTablesURLs[0]
				orphonURL=GeneTablesURLs[1]
				geneTablesBase=base+"/"+organism+"/GeneTables"
				if(not(os.path.isdir(geneTablesBase))):
					os.makedirs(geneTablesBase)
				regularTablePath=geneTablesBase+"/"+locus+".html"
				if(not(os.path.exists(regularTablePath)) or unconditionalForceReplace==True):
					print "Gene table",locus,"for organism",organism,"not found or force replace set to true...so downloading it...."
					print "Downloading gene table",locus,"for organism",organism,"from URL=",regularURL," saving to",regularTablePath
					downloadURLToLocalFileAssumingDirectoryExists(regularURL,regularTablePath)
				else:
					print "Gene table",locus,"for organism",organism," found and force replace set to False...so not downloading it...."
				orphonTablePath=regularTablePath+".orphons.html"
				if(not(os.path.exists(orphonTablePath)) or unconditionalForceReplace==True):
					print "Orphon gene table",locus,"for organism",organism,"not found, so downloading it..."
					print "Downloading orphon gene table",locus,"for organism",organism,"from URL=",orphonURL,"and saving to",orphonTablePath
					downloadURLToLocalFileAssumingDirectoryExists(orphonURL,orphonTablePath)
				else:
					print "Orphon Gene table",orphonTablePath,"for organism",organism," found and force replace set to False...so not downloading it...."

		
	#download the GeneDB and LIGM-DB
	def download_IMGT(self,unconditionalForceReplace=False):
		self.buildAndExecuteWGETDownloadScript()



	#build the Ref. Dir Set files (locus-based files to IgBLAST against )
	#from the downloaded GENEDB files
	def buildRefDirSetsFromGENEDB(self):
		base_fna=self.getBaseDir()+"/"+self.imgt_genedb_rel_path
		if(not(os.path.isfile(base_fna)) or not (os.path.exists(base_fna))):
			gmm="ERROR, GENEDB file not found!  Has it been downloaded?!?!?"
			print gmm
			raise Exception(gmm)
		else:
			organisms=self.getOrganismList(True)
			org_imgt_regex_map=dict()
			org_imgt_regex_map['human']='^homo\s*sapien'
			org_imgt_regex_map['Mus_musculus']='^mus\smusculus'
			genedb_map=read_fasta_file_into_map(base_fna)
			#FOR EACH ORGANISM
			for organism in organisms:
				loci=get_loci_list()
				org_re=re.compile(org_imgt_regex_map[organism],re.IGNORECASE)
				#FOR EACH LOCUS
				for locus in loci:
					segments=get_segment_list()
					#FOR EACH SEGMENT
					#from the base_fna (GENEDB), gather all the key/value pairs for this organism and locus
					org_loci_map=dict()
					for descriptor in genedb_map:
						pieces=descriptor.split("|")
						allele_name=pieces[1]
						org_name=pieces[2]
						reg_name=pieces[4]
						#MATCH BY ORGANISM AND THE REG-NAME (WHICH IS BASED ON SEGMENT/REGION)
						if(allele_name.startswith(locus) and reg_name==locus[3]+"-REGION"    ):
							#okay it matches on the locus and region name, now let's check on the organism
							re_res=re.match(org_re, org_name)
							if(re_res):
								#matched on the organism, so add it to the map!
								org_loci_map[descriptor]=genedb_map[descriptor]
					#write the gathered descriptor/sequence pairs to FASTA files
					dest_dir=self.getBaseDir()+"/"+organism+"/ReferenceDirectorySet/"
					dest_file=dest_dir+"/"+locus+".fna"
					if(not(os.path.isdir(dest_dir))):
						 os.makedirs(dest_dir)
					writer=open(dest_file,'w')
					for descriptor in org_loci_map:
						writer.write(">"+descriptor.strip()+"\n"+org_loci_map[descriptor].strip()+"\n")
					writer.close()
					
			

	#download from IMGT
	def buildAndExecuteWGETDownloadScript(self):
		if(not(os.path.isdir(self.getBaseDir()))):
			print "ERROR, failed to find directory ",self.getBaseDir(),"!"
			print "Skipping download of imgt.org annotations/database files!"
			return
		refDBURL="http://www.imgt.org/download/"
		wgetCMD="cd "+self.db_base+"\n"
		wgetCMD+="/usr/bin/wget -r -np "+refDBURL+"GENE-DB/ "+refDBURL+"/LIGM-DB/\n"
		#we're not using VBASE any more!
		#wgetCMD+="/usr/bin/wget -m http://www.vbase2.org\n"
		#vbase_down="www.vbase2.org"
		#don't download these fetch fasta from the data record queries
		#down_human_vbase_cmd="/usr/bin/wget -O "+vbase_down+"/humanall.fasta http://www.vbase2.org/download/humanall.fasta"
		#down_mouse_vbase_cmd="/usr/bin/wget -O "+vbase_down+"/mouseall.fasta http://www.vbase2.org/download/mouseall.fasta"
		#wgetCMD+="if [ -d \""+vbase_down+"\" ] ; then "+down_human_vbase_cmd+" ; "+down_mouse_vbase_cmd+" ; else echo \"APPARENT ERROR IN DOWNLOADING VBASE!\" ; fi ;\n"
		uncomp_cmd="echo \"Now searching for .Z compressed files to uncompress ...\" ; for COMPRESSED in `find "+str(self.db_base)+"|grep -P '\.Z$'` ; do UNCOMPRESSED=`echo $COMPRESSED|sed -r \"s/\.Z//gi\"` ;    echo \"Found compressed file $COMPRESSED ... to uncompress it to $UNCOMPRESSED ...\" ; echo \"USING command uncompress -c $COMPRESSED > $UNCOMPRESSED\" ; uncompress -c $COMPRESSED > $UNCOMPRESSED ; done ;\n"
		wgetCMD+=uncomp_cmd
		wgetScriptPath=self.db_base+"/wgetscript.sh"
		wgetScriptLogBase=wgetScriptPath+".log"
		wgetScriptOutLog=wgetScriptLogBase+".out"
		wgetScriptErrLog=wgetScriptLogBase+".err"
		write_temp_bash_script(wgetCMD,wgetScriptPath)
		print "Proceeding to download "+refDBURL+"GENE-DB/ and "+refDBURL+"/LIGM-DB/ ...  (see logs at "+wgetScriptLogBase+".*)"
		execute_bash_script(wgetScriptPath,outPath=wgetScriptOutLog,errPath=wgetScriptErrLog)
	


	#given a full fasta descriptor, get the record from IMGT.dat
	#note that the descriptor should not have the ">" at the beginning
	def extractIMGTDatRecordUsingRefDirSetDescriptor(self,descriptor,biopythonRec=False):
		if(self.db_base==None):
			raise Exception("Error, db_base is not, must initialize first!")
		pieces=descriptor.split("|")
		#BN000872|IGHV5-9-1*02|Mus musculus_C57BL/6|F|V-REGION|2334745..2335041|294 nt|1| | | | |294+24=318| | |
		#print "got pieces ",pieces
		accession=pieces[0]
		#print "accession=",accession
		ss=self.getStartStopFromIndexGivenAccession(accession)
		#print "got start/stop : ",ss
		if(len(ss)==2):
			#regular accession
			start=ss[0]
			stop=ss[1]
			#return self.fetchBioPythonRecFromDat(start,stop,biopythonRec)
			return self.fetchRecFromDat(start,stop,biopythonRec)
		else:
			#irregular accession, use descriptor and index to find the correct accession
			raise Exception("Error :  irregular descriptor : "+descriptor+" failed to obain imgt.dat start/stop of accession  "+accession)





	#read index into dict/map
	def cacheIndex(self,indexPath=None,existingCache=None):
		if(not(self.accession_start_stop_map==None)):
			#it's already initialized
			return
		if(self.indexPath is None):
			self.indexPath=str(str(self.imgt_dat_path)+str(self.db_idx_extension))
		#print "Cacheing index for imgt.dat file ",self.imgt_dat_path
		if(not(os.path.exists(str(self.indexPath)))):
			#print "Creating the index",str(self.indexPath)," first because it doesn't exist..."
			self.indexIMGTDatFile(self.imgt_dat_path,self.indexPath)
		else:
			#print "The index file ",str(self.indexPath)," was found, so no need to re-create it...."
			pass
		if(not(existingCache==None)):
			if(self.accession_start_stop_map==None):
				self.accession_start_stop_map=dict()
			self.accession_start_stop_map=merge_maps(self.accession_start_stop_map,existingCache)
		else:
			self.accession_start_stop_map=dict()
		idxReader=open(self.indexPath,'r')
		for line in idxReader:
			line=line.strip()
			pieces=line.split("\t")
			accession=pieces[0]
			if(len(pieces)==3):
				ss=[pieces[1],pieces[2]]
			else:
				ss=[pieces[1]]
			self.accession_start_stop_map[accession]=ss
		idxReader.close()


	#get the IMGT record given an allele name
	def getIMGTDatGivenAllele(self,a,biopythonrec=False,org="human"):
		#print "passed : ",a," with org="+str(org)
		descriptor=self.extractDescriptorLine(a,org)
		#print "got descriptor = ",descriptor
		imgtDAT=self.extractIMGTDatRecordUsingRefDirSetDescriptor(descriptor,biopythonrec)
		return imgtDAT


	#given the allele name and organism, get the descriptor
	def getIMGTDescriptor(self,allele_name,org):
		locus=allele_name[0:4]
		fasta_path=self.getDirBase()+"/"+org+"/ReferenceDirectorySet/"+locus+".fna"
		fasta_dict=read_fasta_file_into_map(fasta_path)
		desired_descriptor=None
		for desc in fasta_dict:
			desc_pieces=desc.split("|")
			desc_allele=desc_pieces[1]
			if(desc_allele==allele_name):
				desired_descriptor=desc
		return desired_descriptor

	#get the IMGT region start and stop from IMGT given the allele_name, organism, and region
	def getRegionStartStopFromIMGTDat(self,allele_name,org,region):
		desired_descriptor=self.getIMGTDescriptor(allele_name,org)
		if(not(desired_descriptor is None)):
			imgt_dat_rec=self.getIMGTDatGivenAllele(allele_name,True,org)
			start_stop=self.getStartStopFromIMGTDesc(desired_descriptor)
			print "The desired descriptor is ",desired_descriptor
			reg_int=self.getRegionStartStopFromIMGTDatRecAndRange(imgt_dat_rec,True,start_stop[0],start_stop[1],region)
			return reg_int

	#given an interval range and an imgtdata record, fetch
	#the start/stop of the specified region
	def getRegionStartStopFromIMGTDatRecAndRange(self,imgt_dat,isBiopython,range_start,range_end,region_name):
		if(not(isBiopython)):
			lines=imgt_dat.split("\n")
			reg_regex=re.compile("^FT\s+"+region_name+"[^\s]*\s+<?(\d+)\.+(\d+)>?\s*$")
			for line in lines:
				search_res=re.search(reg_regex,line)					
				if(search_res):
					start=min(int(search_res.group(1)),int(search_res.group(2)))
					end=max(int(search_res.group(1)),int(search_res.group(2)))
					if(range_start<=start and end<=range_end):
						return [start,end]
		else:
			feature_list=imgt_dat.features
			for feature in feature_list:
				#print "got a feature : ",feature
				ftype=feature.type
				#print "the type is ",ftype	
				qualifiers=feature.qualifiers
				#print "qualifiers : ",qualifiers
				location=feature.location
				#print "location : ",location
				l_start=int(location.start)+1	#add 1 cause it's 0-based
				l_end=int(location.end)+1	#add 1 cause it's 0-based
				strand=location.strand
				#print "The strand is ",strand
				#print "The l_start and l_end are : ",l_start," and ",l_end
				r_start=min(l_start,l_end)
				r_end=max(l_start,l_end)
				r_len=abs(r_start-r_end)
				if(ftype.startswith(region_name)):
					if(range_start-1<=r_start and r_end<=range_end+1):
						ret_start=r_start-range_start+1
						ret_end=ret_start+r_len-1
						if(strand==1):
							return [ret_start,ret_end]
						else:
							print "NEED RETURN REVERSE!"
							range_len=abs(range_start-range_end)
							new_start=range_len-ret_end+2
							new_end=new_start+r_len-1
							return [new_start,new_end]
							#sys.exit(0)
					else:
						pass
						#print "Though it's the right region it falls out of range .  The range is (",range_start,",",range_end,")"
						#print "location : ",r_start,",",r_end
						#print "Start test : ",range_start,"<=",r_start, "???"
						#if(range_start<=r_start):
						#	print "start test success"
						#else:
						#	print "start test fail"
						#print "End test : ",r_end,"<=",range_end,"????"
						#if(r_end<=range_end):
						#	print "end test success"
						#else:
						#	print "end test fail"
						
				else:
					#print "type ",ftype," isn't the desired region...."
					pass
		return [(-1),(-1)]			



	#get start/stop from index given accession
	def getStartStopFromIndexGivenAccession(self,a):
		if(self.accession_start_stop_map==None):
			self.cacheIndex(self.imgt_dat_path+self.db_idx_extension)
		if(a in self.accession_start_stop_map):
			return self.accession_start_stop_map[a]
		else:
			pass 
		

	#given a complete descriptor and organism, fetch the corresponding reference directory set sequence
	def getRefDirSetFNAGivenCompleteDescriptor(self,descriptor,organism):
		#do a dict lookup to return the data if the routine has been used before
		#NOTE that using the descriptor as the dict key addresses the ORGANISM not being used cause the descriptor contains the organism!
		if(self.ref_dir_set_desc_seqs_map==None):
			self.ref_dir_set_desc_seqs_map=dict()
		if(descriptor in self.ref_dir_set_desc_seqs_map):
			return self.ref_dir_set_desc_seqs_map[descriptor]


		#if the routine hasn't been used before, load all the data!
		myloci=get_loci_list()
		for locus in myloci:
			html_fna_path=self.db_base+"/"+organism+"/ReferenceDirectorySet/"+locus+".fna"
			fasta_recs=read_fasta_file_into_map(html_fna_path)
			for fasta_desc in fasta_recs:
				self.ref_dir_set_desc_seqs_map[fasta_desc]=fasta_recs[fasta_desc].upper()
		if(descriptor in self.ref_dir_set_desc_seqs_map):
			return self.ref_dir_set_desc_seqs_map[descriptor]
		else:
			raise Exception("Error, descriptor "+descriptor+" points to a non-existent reference directory set sequence under "+organism+"???")

	#get the sequence optionally removing gaps
	def getRefDirSetFNASeqGivenOrgAndAllele(self,allele_name,organism,removeGaps=True):
		#get the descriptor
		subject_descriptor=self.extractDescriptorLine(allele_name,organism)
		#use the descriptor to get the sequence
		subject_sequence=self.getRefDirSetFNAGivenCompleteDescriptor(subject_descriptor,organism)
		if(removeGaps):
			if(self.org_allele_desc_gapless_map!=None):
				if(allele_name in self.org_allele_desc_gapless_map):
					if(organism in self.org_allele_desc_gapless_map[allele_name]):
						return self.org_allele_desc_gapless_map[allele_name][organism]
					else:
						#do the lookup below then add it!
						pass
				else:
					#allele_name not here yet....
					self.org_allele_desc_gapless_map[allele_name]=dict()
			else:
				self.org_allele_desc_gapless_map=dict()
				self.org_allele_desc_gapless_map[allele_name]=dict()
			repPattern=re.compile(r'[^A-Za-z]')
			subject_sequence=re.sub(repPattern,"",subject_sequence)
			self.org_allele_desc_gapless_map[allele_name][organism]=subject_sequence
		return subject_sequence



	#given an allele name and an organism string, extract the fasta descriptor with the specified allele name
	#use dictionary for cache purposes
	def extractDescriptorLine(self,allele_name,org="human"):
		#print "inside EDL allele_name=",allele_name," org=",org
		#sys.exit(0)
		#this little code here does a cache lookup
		if(not(self.org_allele_name_desc_map==None)):
			#print "using EDL dict to see if org=",org," is in it...."
			if(org in self.org_allele_name_desc_map):
				if(allele_name in self.org_allele_name_desc_map[org]):
					#print "Using cache lookup for EDL with org=",org,"and to return ",self.org_allele_name_desc_map[org][allele_name]
					#sys.exit(0)
					return self.org_allele_name_desc_map[org][allele_name]
				else:
					#print "INNermost EDL cache test fails"
					pass
			else:
				#print "org ain't in EDL map!"
				pass
		else:
			#print "initing EDL dict..."
			self.org_allele_name_desc_map=dict()
		#sys.exit(0)
		if(self.db_base==None):
			#self.db_base=self.db_base
			raise Exception("Error, db_base not set! Did you initialize???")
		org_dir=self.db_base+"/"+org
		to_be_returned=None
		if(os.path.isdir(org_dir)):
			fna_glob_str=org_dir+"/ReferenceDirectorySet/*.fna"
			fna_files=glob.glob(fna_glob_str)
			for fna_file in fna_files:
				fna_reader=open(fna_file,'r')
				for fna_line in fna_reader:
					if(fna_line.startswith(">")):
						descriptor=fna_line[1:]
						if(not(org in self.org_allele_name_desc_map)):
							self.org_allele_name_desc_map[org]=dict()
						pieces=descriptor.split("|")
						if(len(pieces)>1):
							descriptor_allele=pieces[1]
							if(descriptor_allele.strip()==allele_name.strip()):
								to_be_returned=descriptor.strip()
							self.org_allele_name_desc_map[org][descriptor_allele.strip()]=descriptor.strip()
			if(not(to_be_returned==None)):
				return to_be_returned
			else:
				raise Exception("Error, descriptor with allele name = '"+str(allele_name)+"' not found under "+str(self.db_base)+" for organism = "+str(org))
		else:
			raise Exception("Error, invalid organism="+str(org)+", its directory doesn't exist under"+str(self.db_base)+"!")



	#from an IMGT descriptor (">sfsdf|sdklfdjsf|dklsfjlsf|....") extract the pices
	def extractIMGTDescriptorPieces(self,dp):
		if(dp[0]==">"):
			new_desc=dp[1:]
		else:
			new_desc=dp
		pieces=new_desc.split("|")
		return pieces



	#from IMGT descriptor return an array of start/stop where start<=stop
	def getStartStopFromIMGTDesc(self,desc):
		pieces=desc.split("|")
		return self.getStartStopFromIMGTDescPieces(pieces)

	#from IMGT descriptor pieces return an array of start/stop where start<=stop
	def getStartStopFromIMGTDescPieces(self,desc_pieces):
		#>M13911|IGHV1-NL1*01|Homo sapiens|P|V-REGION|125..420|296 nt|1| | | | |296+24=320| |rev-compl|
		#>D87017|IGLJ5*02|Homo sapiens|ORF|J-REGION|11386..11423|38 nt|2| | | | |38+0=38| |rev-compl|
		interval_piece=desc_pieces[5]
		ss_re='(\d+)\.+(\d+)'
		interval_re=re.compile(ss_re)
		search_res=re.search(interval_re,interval_piece)
		if(search_res):
			first=int(search_res.group(1))
			second=int(search_res.group(2))
			start=min(first,second)
			end=max(first,second)
			t=[start,end]
			return t
		else:
			raise Exception("Error, from pieces "+str(desc_pieces)+" (piece='"+interval_piece+"' with regex="+ss_re+") unable to retrieve start/stop!")




	#fetch a record given position interval from the imgt.dat file
	def fetchRecFromDat(self,start,stop,biopython=False,idxpath=None):
		if(idxpath==None):
			idxpath=self.imgt_dat_path
		#print "i want to open ",idxpath
		if(not(stop>start)):
			return ""
		reader=open(idxpath,'r')
		reader.seek(int(start))
		if(biopython):
			records=SeqIO.parse(reader,"imgt")
			for record in records:
				my_rec=record
				reader.close()
				return my_rec
		data=reader.read(int(stop)-int(start))
		reader.close()
		return data

	#given a record from imgt.dat, find the first region interval
	def getRefVRegionInterval(self,data,region_name):
		lines=data.split("\t")
		for line in lines:
			#FT   CDR3-IMGT           371..412
			reg_regex="^FT\s+"+region_name+"[^\s]*\s+<?(\d+)\.+(\d+)>?\s*$"
			search_res=re.search(reg_regex,line)
			if(search_res):
				start=search_res.group(1)
				end=search_res.group(2)
				return [start,end]
		return None



	#index the imgt.dat file
	def indexIMGTDatFile(self,filepath=None,indexfile=None):
		if(filepath==None):
			filepath=self.db_base+"/"+self.imgt_dat_rel_path
		if(indexfile==None):
			indexfile=filepath+self.db_idx_extension
		print "Creating index file from ",filepath," ... writing index to ",indexfile
		reader=open(filepath,'r')
		acc_re=re.compile(r'^ID\s+([A-Z0-9]+)[^A-Z0-9]+')
		embl_tpa_re=re.compile(r'^DR\s+EMBL\-TPA;\s+([^\s]+)\.\s*$')
		current_embl_tpa=None
		current_accession=None
		rec_start=None
		rec_end=None
		index_file=open(indexfile,'w')
		rec_num=0
		flag=True
		while(flag):
			line=reader.readline()
			if(line):
				#print "Got ",line.strip()
				rs=re.search(acc_re,line)
				es=re.search(embl_tpa_re,line)
				if(rs):
					current_accession=rs.group(1)
					rec_start=reader.tell()-len(line)
					if(self.accession_dat_file_map==None):
						self.accession_dat_file_map=dict()
					self.accession_dat_file_map[current_accession]=filepath
				elif(es):
					current_embl_tpa=es.group(1)
				elif(line.startswith("//")):
					rec_end=reader.tell()-1
					index_file.write(current_accession+"\t"+str(rec_start)+"\t"+str(rec_end)+"\n")
					if(not(current_embl_tpa==None)):
						#index_file.write(current_embl_tpa+"\t"+str(rec_start)+"\t"+str(rec_end)+"\n")
						index_file.write(current_embl_tpa+"\t"+current_accession+"\n");
					current_embl_tpa=None
				#print "current accession = ",current_accession
			else:
				flag=False
		index_file.close()	
	



