#!/usr/bin/env python
from datetime import datetime
import time
import pprint, pickle
import os
import subprocess
from subprocess import Popen
#import urllib2
#from urlgrabber import urlopen
import re
from subprocess import call
#import nwalign as nw
import pickle
import glob
import sys


#used with argparse
def extractAsItemOrFirstFromList(t):
	if(type(t)==list):
		if(len(t)>0):
			return t[0]
	else:
		return t
	return None






#turn an alignment (in the form of an array)
#into a printabel nice string
def getNiceAlignment(aln):
	a=str("")
	for idx, val in enumerate(aln):
		a+=val+"\n"
	return a


#repeat a string (s) n times
def repeatString(s,n):
	if(n<=0):
		return ""
	else:
		r=str("")
		for i in range(n):
			r+=s
		return r
		




#build an alignment printable string
#from the query and subject, but include a MIDLINE
def buildAlignmentWholeSeqsDirect(q,s):
	aln=["","",""]
	aln[0]=str(q)
	aln[2]=str(s)
	for b in range(len(aln[0])):
		qbase=aln[0][b]
		sbase=aln[2][b]
		if(qbase=="-" or sbase=="-"):
			aln[1]+=" "
		elif(not(qbase==sbase)):
			aln[1]+="X"
		else:
			aln[1]+="|"
	aln[0]="QURY:"+aln[0]
	aln[1]="MDLN:"+aln[1]
	aln[2]="SBCT:"+aln[2]
	return aln		


#returns a string given query first, then subject
def printNiceAlignment(q,s):
	joiner="\n"
	return joiner.join(buildAlignmentWholeSeqsDirect(q,s))




#given a BTOP string build a dummy alignment, putting N
#where sequence cannot be known without source data
def buildDummyQAndSSeqsFromBTop(btop):
	#query, then subject
	q_s=["",""]
	if(len(btop)<=0):
		return q_s
	else:
		while(len(btop)>0):
			adm=re.search('^(\d+)$',btop)
			dm=re.search('^(\d+)[^0-9]',btop)
			firstTwoLetters=re.search('^([a-z\\-])([a-z\\-])',btop,re.IGNORECASE)
			if adm:
				#all digital
				q_s[0]+=repStr("N",int(adm.group(1)))
				q_s[1]+=repStr("N",int(adm.group(1)))
				return q_s
			elif(dm):
				#just starts with digits
				q_s[0]+=repStr("N",int(dm.group(1)))
				q_s[1]+=repStr("N",int(dm.group(1)))
				num_digits=len(dm.group(1))
				btop=btop[num_digits:]
			else:
				#2 characters!
				q_s[0]+=firstTwoLetters.group(1)
				q_s[1]+=firstTwoLetters.group(2)
				btop=btop[2:]
		return q_s
			
		




#returns array of query, then match/midline, then subject based on btop alignment specification
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
		#BTOP is 100% digital
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
		#BTOP only starts with digits...
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
		#BTOP starts with a dash and a letter, a letter and a dash, or two letters
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








#given an alignment, and a position of a sequence IN the alignment
#return the porition of the alignment whose bp are <= or >= that position in the alignment
def getAlnAtAndCond(q_aln,s_aln,q_from,q_to,s_from,s_to,a_pos,st="query",cond="geq"):
	if(not(st=="query")):
		st="subject"
	aln=["",""] #Q, then S
	if(not(cond=="leq")):
		cond="geq"
	s_pos=s_from
	q_pos=q_from
	temp=0
	while(temp<len(q_aln)):
		if(cond=="leq"):
			if(st=="query"):
				if(a_pos>=q_pos):
					aln[0]+=q_aln[temp]
					aln[1]+=s_aln[temp]
			else:
				if(a_pos>=s_pos):
					aln[0]+=q_aln[temp]
					aln[1]+=s_aln[temp]
		else:
			if(st=="query"):
				if(a_pos<=q_pos):
					aln[0]+=q_aln[temp]
					aln[1]+=s_aln[temp]
			else:
				if(a_pos<=s_pos):
					aln[0]+=q_aln[temp]
					aln[1]+=s_aln[temp]				
		if(s_aln[temp]!="-"):
			s_pos+=1
		if(q_aln[temp]!="-"):
			q_pos+=1
		temp+=1
	return aln #Q then S
			
	
#does a position in a sequence align to a gap???
#return true or false
def doesPosInSeqAlignToGap(q_aln,s_aln,q_from,q_to,s_from,s_to,a_pos,st="query"):
	if(st!="query"):
		st="subject"
	if(st=="subject" and not(s_from<=a_pos and a_pos<=s_to)):
		return None
	elif(st=="query" and not(q_from<=a_pos and a_pos<=q_to)):
		return None
	else:
		s_pos=s_from
		q_pos=q_from
		temp=0
		while(temp<len(q_aln)):
			s_val=s_aln[temp]
			q_val=q_aln[temp]
			if(st=="query" and a_pos==q_pos):
				if(s_val=="-"):
					return True
				else:
					return False
			if(st=="subject" and a_pos==s_pos):
				if(q_val=="-"):
					return True
				else:
					return False
			if(s_aln[temp]!="-"):
				s_pos+=1
			if(q_aln[temp]!="-"):
				q_pos+=1
			temp+=1
		return aln #Q then S		





#given a length l, make a list of that
#length whose entries are empty strings
def makeEmptyArrayOfStringsOfLen(l):
	empty_str_arr=list()
	#non-negative numbers only please!
	l=max(0,l)
	for i in range(l):
		empty_str_arr.append("")
	return empty_str_arr


#given a length l, make a list of digits of
#length whose entries are the given digit or 0 if not given
def makeEmptyArrayOfDigitsOfLen(l,d=0):
	empty_dig_arr=list()
	#non-negative numbers only please!
	l=max(0,l)
	for i in range(l):
		empty_dig_arr.append(d)
	return empty_dig_arr

#given a hash map, set all values
# to a given value 
def makeAllMapValuesVal(m,v):
	if(m is not None):
		for k in m:
			m[k]=v
	return m
		


#translate a sequence with BIOPYTHON
#using IUPAC standards
#N-padded to achieve length of multiple of 3
#biopython returns X for ambiguous, but the AA for non-ambiguous
def biopythonTranslate(dna):
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	num_extraN=0
	if(len(dna)==0):
		num_extraN=3
	elif(len(dna)%3!=0):
		num_extraN=3-(len(dna)%3)
	extraNs=("N"*num_extraN)
	#print "nen=",num_extraN
	#print "en=",extraNs
	dna=dna+extraNs
	#print "dna=",dna
	coding_dna = Seq(dna, IUPAC.ambiguous_dna)
	#print "coding_dna" ,  coding_dna
	template_dna = coding_dna.reverse_complement()
	#print "template_dna" , template_dna
	messenger_rna = coding_dna.transcribe()
	#print "messenger_rna" , messenger_rna
	back_transcripted=messenger_rna.back_transcribe()
	#print "back_transcripted" , back_transcripted
	translation=messenger_rna.translate()
	#print "translation ", translation
	return translation






def isIntegral(s):
	ire=re.compile(r'^\s*(\d+)\s*$')
	sr=re.search(ire,s)
	if(sr):
		return True
	else:
		return False



def countFastaReads(fastaFile):
	f=open(fastaFile,'r')
	fc=0
	for line in f:
		if(line.startswith(">")):
			fc+=1
	f.close()
	return fc




def fileLineCount(f):
	reader=open(f,'r')
	line_count=0
	for line in reader:
		line_count+=1
	reader.close()
	return line_count


#using a map and key 
#increment the value
#init at 1 if no such key is in the map
def counter_map_inc(counter_map,key):
	if key in counter_map:
		counter_map[key]+=1
	else:
		counter_map[key]=1
	return counter_map



#scan to see if a string contains an apostrophe
def doesStringContainAnApostrophe(s):
	if(re.search("'",s)):
		return True
	else:
		return False


#return TRUE if the first argument (a)
# is a substring of the second (b)
def a_subseq_of_b(a,b):
	try:
		a_index=b.index(a)
		return True
	except ValueError:
		return False




#return true if the two sequences in an alignment are equal
#allow exclusion of head/tail indels as an option (defaults to False)
def isAlignmentFreeOfBothMutationsAndIndels(alignment,TryRemovalOfHeadAndTailDash=False):
	if(alignment[0]==alignment[1]):
		#query and subject are the same, return true
		return True
	else:
		if(TryRemovalOfHeadAndTailDash and (removeHeadAndTailMultiDash(alignment[0])==removeHeadAndTailMultiDash(alignment[1]))):
			#see if they're equal if head/tail indels are removed
			return True
		#nope!
		return False



#read multiple fasta files into a single map
def readMultipleMapFilesIntoSingleMapWithGlob(g):
	map_files=glob.glob(g)
	total_map=dict()
	for map_file in map_files:
		single_map=read_map_from_file(map_file)
		for key in single_map:
			total_map[key]=single_map[key]
	return total_map


#given a subset list 'subset'
#write data from the source (input)
# that is in the subset to the output path
def fastaSubset(inputFastaPath,subset,outputFastaPath):
	#print "**************************************"
	#print "FASTA SUBSET CALLED"
	#print "FASTA=",inputFastaPath
	q_set=read_fasta_file_into_map(inputFastaPath)
	#for q in q_set:
	#	print "q=",q," len =",len(q_set[q])
	#print "len fasta =",len(q_set)
	#print "SUBSET:",

	#print "size of subset : ",len(subset)
	q_set=read_fasta_file_into_map(inputFastaPath)
	subset_writer=open(outputFastaPath,'w')
	for item in subset:
		#print "looking at item ",item," in subset to determine if it's in fasta...."
		if item in q_set :
			#print "it's in fasta...now writing it to ",outputFastaPath
			subset_writer.write(">"+item+"\n")
			subset_writer.write(q_set[item]+"\n")
		else:
			pass
			#print "it's not in fasta!"
		#print "\n\n"
	subset_writer.close()


#using NW alignment data
#map unmapped items to a list
#(this is for comparing IGBLAST and IMGT data)
def nw_map_from_fastas(q_fasta_path,s_fasta_path,logPath,clone_names,alleleList,mappedPath,unmappedPath,resEq=True):
	q_map=read_fasta_file_into_map(q_fasta_path)
	db_map=read_fasta_file_into_map(s_fasta_path)
	mapped_map=find_best_nw_from_maps_KeepPerfects(q_map,db_map,logPath,clone_names,alleleList,resEq)
	mapped_writer=open(mappedPath,'w')
	unmapped_writer=open(unmappedPath,'w')
	#write the mapped file
	write_map_to_file(mapped_map,mappedPath)
	#write the unammped list
	unmapped_list=list()
	for item in q_map:
		if(not(item in mapped_map)):
			unmapped_list.append(item)
	write_list_to_file(unmapped_list,unmappedPath)
	

#from the map/alignment data and clone list
#find the alignments
def find_best_nw_from_maps_KeepPerfects(q_map,db_map,logPath,clone_names_map,alleleList,resEq=True):
	score_map_name=dict()
	score_map_val=dict()
	best_aln=dict()
	#first, get the best scores
	for q_desc in q_map:
		#nw.global_align("CEELECANTH", "PELICAN", matrix="/home/esalina2/PAM250")
		best_score=(-1)
		gap_open=-5
		gap_extend=-2
		pam_path="/home/esalina2/PAM250"
		for d_desc in db_map:
			alignment=nw.global_align(q_map[q_desc],db_map[d_desc],matrix=pam_path)
			score=nw.score_alignment(alignment[0],alignment[1],gap_open,gap_extend,matrix=pam_path)
			#if(q_desc=="IGHD1-7*01"):
				#print "qkey=",q_desc
				#print "qval=",q_map[q_desc]
				#print "dkey=",d_desc
				#print "d_val=",db_map[d_desc]
				#print "Now attempting a global alignment with :",q_map[q_desc]+", and "+db_map[d_desc]
				#print "alignment=",alignment
				#print "score=",score
			if q_desc not in score_map_val:
				score_map_val[q_desc]=score
				score_map_name[q_desc]=d_desc
				best_aln[q_desc]=alignment
			else:
				existing_score=score_map_val[q_desc]
				if(existing_score<score):
					best_aln[q_desc]=alignment
					score_map_val[q_desc]=score
					score_map_name[q_desc]=d_desc
	good_map=dict()
	for name in score_map_name:
		names_eq=(score_map_name[name]==name)
		seqs_eq=(q_map[name]==db_map[score_map_name[name]])
		perfect_aln=isAlignmentFreeOfBothMutationsAndIndels(best_aln[name])
		perfect_aln_except_for_head_tail_indel=isAlignmentFreeOfBothMutationsAndIndels(best_aln[name],True)
		best_alignment_no_HTIndels=map(removeHeadAndTailMultiDash,best_aln[name])
		aln_frac=float(len(best_alignment_no_HTIndels[0]))/float(min(len(q_map[name]),len(db_map[score_map_name[name]])))
		try:
			log=open(logPath,'a')
			log.write("query("+name+")="+q_map[name]+"\n")
			log.write("ref("+score_map_name[name]+")="+db_map[score_map_name[name]]+"\n")
			log.write("clone info : ")
			if(name in clone_names_map):
				log.write(clone_names_map[name])
			else:
				log.write("no clone name!")
			log.write("\n")
			log.write(best_aln[name][0]+"\n"+best_aln[name][1]+"\n")
			log.write("score="+str(score_map_val[name])+"\n")
			if(seqs_eq or perfect_aln):
				log.write("The sequences are equal or the alignment is perfect!")
			else:
				log.write("Either the sequences are unequal and the alignment is imperfect!")
			inAllelListStatus=None
			if(extractIMGTNameFromKey(score_map_name[name]) in alleleList):
				inAllelListStatus=True
			else:
				inAllelListStatus=False
			log.write("IMGT hit in allele list status : "+str(inAllelListStatus))
			log.write("\n\n\n")
			log.close()
		except:
			print "ERROR WRITING ALGINMENT LOG FILE ",logPath
			print "Exception in user code:"
			print '-'*60
			traceback.print_exc(file=sys.stdout)
			print '-'*60
		if((seqs_eq or perfect_aln) and (extractIMGTNameFromKey(score_map_name[name]) in alleleList) ):
			good_map[name]=extractIMGTNameFromKey(score_map_name[name])
	return good_map




#act as a wrapper getting a value
#from a map with 0 as the default value
#if the key is not in the map
def counter_map_get(counter_map,key):
	if key in counter_map:
		return counter_map[key]
	else:
		return 0



#wrapper to touch a file
#see http://stackoverflow.com/questions/1158076/implement-touch-using-python
def touch(fname, times=None):
	try:
		with file(fname, 'a'):
			os.utime(fname, times)
	except Exception:
		print "Error in touching file",fname


#given a file path attempt to find out if it does NOT exist
#would be touchable and would be deletable
def test_nonexistent_touchable_and_del(f):
	if(os.path.exists(f)):
		return False
	touch(f)
	if(not os.path.exists(f)):
		return False	
	os.remove(f)
	if(os.path.exists(f)):
		return False
	return True



#get the IMGT name from the RDS descriptor
def getIMGTNameFromRefDirSetDescriptor(desc):
	pieces=desc.split('|')
	imgt_name=pieces[1]
	return imgt_name


#given a list of IMGT descriptors
#return the list back, but just the IMGT allele names
def getIMGTNameListFromFastaMap(fm):
	imgtList=list()
	for k in fm:
		imgt_name=getIMGTNameFromRefDirSetDescriptor(k)
		imgtList.append(imgt_name)
	return imgtList




#given a tree object, recursively explore it and
#aggregate a list of all the alleles (ABC*01)
#and return that list
def get_list_of_alleles_appearing_in_tree(t):
	listOfAlleles=list()
	for k in t:
		name=str(k)
		if(looksLikeAlleleStr(name)):
			listOfAlleles.append(k)
		recList=get_list_of_alleles_appearing_in_tree(t[k])
		listOfAlleles=listOfAlleles+recList
	return listOfAlleles



def get_segment_list():
	return ["V","D","J"]




def blastFormatFNAInDir(d,blastformtexecpath):
	fnas=glob.glob(d+"/*.fna")
	for fna in fnas:
		if(fna.find("_D")!=(-1)):
			#needing parse_seqids seems to be needed for D segment
			format_cmd=blastformtexecpath+" -parse_seqids -dbtype nucl -in "+fna
		else:
			# parse_seqids seems to be causing issues for V segment
			format_cmd=blastformtexecpath+" -dbtype nucl -in "+fna
		if(True):
			format_cmd=blastformtexecpath+" -dbtype nucl -in "+fna
		print "BLASTDB Formatting ",fna," with ",format_cmd
		out_path=fna+".blastfmt.out"
		err_path=fna+".blastfmt.err"
		out_handle=open(out_path,'w')
		err_handle=open(err_path,'w')
		proc=subprocess.Popen(format_cmd,stdout=out_handle,stderr=err_handle,shell="/bin/bash")
		ret_val=proc.wait()
		out_handle.close()
		err_handle.close()
		



def blastFormatFLEX(fna,blastformtexecpath,parse=False,recLevel=0):
		format_cmd=blastformtexecpath+" -dbtype nucl -in "+fna
		if(parse):
			format_cmd+=" -parse_seqids"
		print "BLASTDB Formatting ",fna," with ",format_cmd
		out_path=fna+".blastfmt.out"
		err_path=fna+".blastfmt.err"
		if(recLevel>1):
			print "WARNING, too much recursion in readFileIntoString, check ",err_path," for errors!"
			sys.exit(0)
		out_handle=open(out_path,'w')
		err_handle=open(err_path,'w')
		proc=subprocess.Popen(format_cmd,stdout=out_handle,stderr=err_handle,shell="/bin/bash")
		ret_val=proc.wait()
		out_handle.close()
		err_handle.close()
		errLineCount=fileLineCount(err_path)
		if(errLineCount>0):
			print "\nWARNING, got error during formatting ",fna,"\n"
			errTxt=readFileIntoString(err_path)
			print errTxt,"\n\n"
			
		#if(errLineCount>0 and parse):
		#	print "WARNING encountered BLAST error when -parse_seqids :"
		#	errTxt=readFileIntoString(err_path)
		#	print errTxt
		#	print "now to retry formatting without '-parse_seqids'"
		#	blastFormatFLEX(fna,blastformtexecpath,False,recLevel+1)
		#elif(errLineCount>0 and not(parse)):
		#	print "WARNING encountered BLAST error when -parse_seqids is OFF :"
		#	errTxt=readFileIntoString(err_path)
		#	print errTxt
		#	print "now to retry formatting wit '-parse_seqids'"
		#	blastFormatFLEX(fna,blastformtexecpath,True,recLevel+1)			
			
			

def runCmdButGetNoStdOutOrStdErr(cmd):
	out_handle=open("/dev/null",'w')
	err_handle=open("/dev/null",'w')
	proc=subprocess.Popen(cmd,stdout=out_handle,stderr=err_handle,shell="/bin/bash")
	ret_val=proc.wait()
	out_handle.close()
	err_handle.close()	



#for all strings in a list, force them to be "allelic"
#(if they're not already allelic) by adding *01 at the end
def allelifyList(l):
	#print "NOW ALLEFYING A LIST, THE INPUT LIST IS :"
	#printList(l)
	goodlist=list()
	for i in l:
		#print "Examining ",i," in the list..."
		#print "GOOD LIST IS CURRENTLY:"
		#printList(goodlist)
		alleleRegex=re.compile(r'\*\d+$')
		if(not (looksLikeAlleleStr(i))):
			#print "NEED TO ALLELIFY IT!"
			goodlist.append(i+"*01")
		else:
			#print "it's already allelified!"
			goodlist.append(i)
	#print "RETURNING A 'GOOD' LIST:"
	#printList(goodlist)
	return goodlist


#remove the *\d+ at the end of a string if it is present
def deAllelifyName(allele):
	alleleRe=re.compile('\*\d+$')
	gene=re.sub(alleleRe,'',allele)
	return gene



#read a fasta file into a map (descriptors as keys, sequences as values)
def read_fasta_file_into_map(fasta_path,alwaysSeqToUpper=True):
	fasta_data=read_fasta(fasta_path,alwaysSeqToUpper)
	fasta_map=read_fasta_into_map(fasta_data)
	return fasta_map

#read a fasta string (from a fasta file) into an array of data
#elements 0,2,... become descriptors
#elements 1,3,... become sequences
def read_fasta_string(fasta_string,alwaysSeqToUpper=True,removeNonIUPAC=True):
	lines=fasta_string.split("\n")
	data=list()
	for line in lines:
		temp_line=line.strip()
		if(temp_line.startswith(">")):
			data.append(temp_line[1:])
			data.append("")
		else:
			if(removeNonIUPAC):
				#SEE http://search.cpan.org/~cjfields/BioPerl-1.6.923/Bio/Tools/IUPAC.pm
				iupac_re=re.compile('[A-Z\*]',re.IGNORECASE)
				temp_line=re.sub(iupac_re,'',temp_line)
			if(alwaysSeqToUpper):
				data[len(data)-1]+=temp_line.upper()
			else:
				data[len(data)-1]+=temp_line
				
	return data





#read a single fasta file into an array of data
#elements 0,2,... become descriptors
#elements 1,3,... become sequences
def read_fasta(fasta_path,alwaysSeqToUpper=True):
	data=list()
	temp_seq=""
	infile=open(fasta_path,'r')
	for line in infile:
		temp_line=line.strip()
		if(temp_line.startswith(">")):
			data.append(temp_line[1:])
			data.append("")
		else:
			if(alwaysSeqToUpper):
				data[len(data)-1]+=temp_line.upper()
			else:
				data[len(data)-1]+=temp_line
	infile.close()			
	return data


#repeat a string a given number of times
def repStr(s,n):
	if(n<=0):
		return ""
	else:
		return s+repStr(s,n-1)



#return a nicely-formatted date and time string
def formatNiceDateTimeStamp():
	#	>>> datetime.now().strftime("%A, %d. %B %Y %I:%M%p")
	#	'Thursday, 09. January 2014 03:49PM'
	#>>> #Thu Jan  9 15:49:55 CST 2014
	#... #
	#... #
	#... '
	#  File "<stdin>", line 4
	#    #'#
	#    ^
	#SyntaxError: EOL while scanning string literal
	return datetime.now().strftime("%A, %B %d, %H:%M %z, %Y")




#remove any "-" (0 or more) from the head and tail of a string
def removeHeadAndTailMultiDash(s):
	s=re.sub(r'\-+$',"",s)
	s=re.sub(r'^\-+',"",s)
	return s



#little utility to return the modes
def get_domain_modes():
	domain_list=["imgt","kabat"]
	#domain_list=["kabat"]
	return domain_list




#remove any ; found at the end of a string
#remove 1 or more of them .... as many are as found
def removeTerminatingSemicolonIfItExists(s):
	s=re.sub(r';+$',"",s)
	return s
	



#wrapper to call makeblastdb to format a blast DB
def format_blast_db(fastaPath,dbtype="nucl",formatdbbin="/usr/local/bin/makeblastdb"):
	format_db_cmd=formatdbbin+" -dbtype "+dbtype+" -in "+fastaPath
	print "The format cmd is ",format_db_cmd
	try:
		call(format_db_cmd.split(" "))
	except Exception, e:
		print "BLAST_FORMAT_CMD error({0})"+format( e.strerror)
	return format_db_cmd


#call NCBI blast
def blast_db(query_path,db_arg,blast_output_path,xmlMode=True):
	blast_cmd="/usr/local/bin/blastn -query "+query_path#+" -word_size 4"
	if(xmlMode):
		blast_cmd+=" -outfmt 5"
	blast_cmd+=" -num_alignments 10 -num_descriptions 10 -max_target_seqs 10 -db "+db_arg+" -out "+blast_output_path
	print "the blast cmd is ",blast_cmd
	try:
		call(blast_cmd.split(" "))
	except Exception, e:
		print "BLAST error({0})"+format( e.strerror)
	return blast_cmd


#read fasta data into a map
#with descriptors as keys
#and squences as values
#input has elements 0,2,4, etc. as descriptors, 1,3,5,...etc at sequences
def read_fasta_into_map(fasta_data):
	i=0
	fasta_map=dict()
	while (i<len(fasta_data)):
		#print "i=",i
		fasta_map[fasta_data[i]]=fasta_data[i+1]
		i+=2
	return fasta_map


#given a 'plain' dict() (python primitive)
#convert it to a default dict
# see http://docs.python.org/2/library/collections.html#collections.defaultdict
# see http://stackoverflow.com/questions/5900578/how-collections-defaultdict-work
def convertToDefaultStringDict(m,s):
	dd=defaultdict(s)
	for k in m:
		dd[k]=m[k]
	return dd



#given a directory hierarchy, unGZ files under it
#(NOTE not functional)...
def uncompressFilesUnderDir(d):
	uncompRE=re.compile(r'\.(Z|gz)$')
	for root, dirs, files in os.walk(d, topdown=False):
		for f in files:
			print "a file=",f
			compressed_path=root+"/"+f
			search_result=re.search(uncompRE,f)
			if(search_result):
				#print "it matches!"
				compEx=search_result.group(1)
				#print "matched on ",compEx
				uncompressed_name=f[0:(len(f)-len(compEx)-1)]
				print "the uncompressed name is ",uncompressed_name
				uncompressed_path=root+"/"+uncompressed_name
				print "the uncompressed path is",uncompressed_path
			else:
				#print "no match"
				pass



#given two maps of keys and values
#return a third that is their union
#if a key is shared mb takes precedence
#and its value overwrites whatever ma may have
def merge_maps(ma,mb):
	mc=dict()
	for a in ma:
		mc[a]=ma[a]
	for b in mb:
		mc[b]=mb[b]
	return mc



#write a list to a file with line numbers (0-based)
def write_list_to_file(m,path,doSort=True):
	writer=open(path,'w')
	if(doSort):
		b=list(m)
		b.sort()
	else:
		b=m
	for i in b:
		writer.write(i+"\n")
	writer.close()


#wrapper to exectute a BASH script
def execute_bash_script(scriptPath,outPath="/dev/stdout",errPath="/dev/stderr"):
	cmd="/bin/bash "+scriptPath.strip()
	#print "popen called with cmd=/bin/bash "+scriptPath+" stdout=",outPath," stderr=",errPath," and shell=/bin/bash"
	out_handle=open(outPath,'w')
	err_handle=open(errPath,'w')
	proc=subprocess.Popen("/bin/bash "+scriptPath,stdout=out_handle,stderr=err_handle,shell="/bin/bash")
	ret_val=proc.wait()
	out_handle.close()
	err_handle.close()
	if(ret_val!=0):
		print "Error in script "+scriptPath
	return ret_val





#wrapper to write a bash script
def write_temp_bash_script(cmd_string,scriptFilePath=None):
	if(scriptFilePath==None):
		scriptFilePath="/tmp/vdj_server_rand_"+binascii.hexlify(os.urandom(16))+".sh"
	try:
		script_writer=open(scriptFilePath,'w')
		script_writer.write("#!/bin/bash\n"+cmd_string+"\n")
		script_writer.close()
	except:
		print "ERROR IN WRITING SCRIPT FILE ",scriptFilePath
	return scriptFilePath	
	

#get the number of lines in a file
def getNumberLinesInFile(p):
	reader=open(p,'r')
	lines=0
	for line in reader:
		lines+=1
	return lines

#alias for getNumberLinesInFile 
#looking at the file as a list
def getLengthOfListInFile(p):
	mylist=read_list_from_file(p)
	return len(mylist)


#return the number of base pairs (throwing out "-")
def getNumberBpInAlnStr(a_str):
	if(a_str.find("-")==(-1)):
		return len(a_str)
	else:
		a_str=re.sub(r'\-','',a_str)
		return len(a_str)


#write a map to a file (optionally using a sorted list of keys)
# and prefix line/key numbers to the output
def write_map_to_file(m,path,doSort=True):
	writer=open(path,'w')
	if(doSort):
		key_set=list(m.keys())
		key_set.sort()
	else:
		key_set=m.keys()
	for k in key_set:
		writer.write(k+"\t"+m[k]+"\n")
	writer.close()


#read a list of lines from a file into a list
def read_list_from_file(path):
	reader=open(path,'r')
	file_list=list()
	for line in reader:
		line=line.strip()
		file_list.append(line)
	reader.close()	
	return file_list

#read a map from a file
#formatted as tab-separated values
#ignore lines that aren't 2-columns (key and value)
def read_map_from_file(path):
	reader=open(path,'r')
	m=dict()
	for line in reader:
		line=line.strip()
		pieces=line.split('\t')
		if(len(pieces)==2):
			m[pieces[0]]=pieces[1]
	reader.close()
	return m



#print a map to the screen (stdout)
#optionally sorting the keys before doing so
def printMap(m,sortKeys=False):
	i=0
	if(not(sortKeys)):
		for k in m:
			print "#"+str(i+1),k,"->",m[k]
			i+=1
	else:
		keys=m.keys()
		#print "pre sort ",keys
		keys.sort()
		#print "pst sort ",keys
		for k in keys:
			print "#"+str(i+1),k,"->",m[k]
			i+=1

#get the time from the DATE command
def returnTimeStamp():
	date_cmd="date"
	date_cmd_args="+%m_%d_%Y_%H_%M_%S.%N"
	date_result=subprocess.check_output([date_cmd,date_cmd_args])


#print a LIST to STDOUT
#using indices
def printList(l):
	for i in range(len(l)):
		print i,":",l[i]



#given two lists and two titles
#print a "side-by-side" style "diff -sy"-like
#output showing what's in one list but not the
#other!
def extendedSortedSetDiff(a,b,a_title,b_title):
	#coveredInA=list()
	blank=""
	print a_title,"\tD\t",b_title
	for a_i in a:
		#coveredInA.append(a_i)
		if(a_i in b):
			print a_i,"\t\t\t",a_i
		else:
			print a_i,"\t<\t\t",blank
	for b_i in b:
		if(not(b_i in a)):
			print blank,"\t>\t",b_i


#given two lists and two titles
#print  messages show differences
def briefSetDiff(a,b,a_title,b_title):
	if(len(a)==len(b) and set(a)==set(b)):
		print "Two sets are the same!"
	else:
		for a_i in a:
			if(not(a_i in b)):
				print "item",a_i," is in ",a_title,"but not in",b_title
		for b_i in b:
			if(not(b_i in a)):
				print "item",b_i," is in ",b_title,"but not in",a_title				


#read a file into a string
def readFileIntoString(path):
	f=open(path,'r')
	data=f.read()
	f.close()
	return data



def readFileURL(url):
	url=url.strip()
	file_part="file://"
	url=url[len(file_part):]
	reader=open(url,'r')
	data=reader.read()
	return data

def isFileURL(url):
	url=url.strip()
	if(url.startswith("file://")):
		return True
	else:
		return False


#read a URL using wget or local
def readAURL(url):
	#f = urllib2.urlopen(url,None,60)
	#html=f.read()
	#from urlgrabber import urlopen
	#swith to urlgrabber!
	#fo = urlopen(url)
	#data = fo.read()
	#fo.close()
	if(isFileURL(url)):
		return readFileURL(url)
	url=url.strip()
	cmd_and_args=["wget","-O","/dev/stdout",url]
	data = Popen( cmd_and_args  , stdout=subprocess.PIPE).communicate()[0]
	return data

#assuming a directory exists write a string to a file whose path is in
#that directory
def writeStringToFilePathAssumingDirectoryExists(string,file_path):
	outfile=open(file_path,'w')
	outfile.write(string)
	outfile.close()



#download a URL and write it to a local file assuming the 
#directory where the file will be written to exists
def downloadURLToLocalFileAssumingDirectoryExists(url,local):
	print "\n\n******************************"
	print "Trying to read URL ",url
	print "and save to",local
	print "******************************\n\n"	
	url_content=None
	url_content=readAURL(url)
	writeStringToFilePathAssumingDirectoryExists(url_content,local)





#return true if a string ends with *XY with X and Y digits
def looksLikeAlleleStr(a):
	are=re.compile(r'\*\d+$')
	res=re.search(are,a)
	if(res):
		return True
	else:
		return False





#uncompress a .Z file using the uncompress program
#uncompress in place or write to the file
#without the .Z extension
#if no .Z extension, then return
#if the output already exists, then return
def uncompressZFile(path,uncompressInPlace=False):
	if(uncompressInPlace):
		cmd="uncompress "+path
		p1 = subprocess.Popen(cmd.split(' '))
		p1.wait()
	else:
		noZ=re.sub(r'\.Z$','',path)
		if(noZ==path):
			#wasn't a Z file?!
			return
		if(os.path.exists(noZ)):
			print "uncompressed file "+noZ+" already exists!"
			return
		cmd="uncompress -c "+path
		#http://stackoverflow.com/questions/19020557/redirecting-output-of-pipe-to-a-file-in-python
		p1 = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE)
		uncomp_writer=open(noZ,'w')
		out_line = p1.stdout.readline()
		uncomp_writer.write(out_line+"\n")
		while out_line:
			out_line=p1.stdout.readline()
			uncomp_writer.write(out_line+"\n")
		p1.stdout.close()
		uncomp_writer.close()
		p1.wait()
		
		



def test():
	tb=">dsfsdf|1111|dslfkjsdf|dklsfjdsf"
	imgt=extractIMGTNameFromKey(tb)
	print "the extraction is ",imgt
	print "Running test...."
	timestamp=returnTimeStamp()
	print "A timestamp is ",timestamp
	a_str="art"
	b_str="dart"
	if(a_subseq_of_b(a_str,b_str)):
		print a_str,"is a substring of",b_str
	else:
		print a_str,"is NOT a substring of",b_str
	if(a_subseq_of_b(b_str,a_str)):
		print b_str,"is a substring of",a_str
	else:
		print b_str,"is NOT a substring of",a_str
	zpath="/tmp/imgt_down/www.imgt.org/download/LIGM-DB/imgt.dat.Z.copy.Z"
	uncompressZFile(zpath,True)





if (__name__=="__main__"):
	import sys
	test()
