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







#given a length l, make a list of digits of
#length whose entries are the given digit or 0 if not given
def makeEmptyArrayOfDigitsOfLen(l,d=0):
	empty_dig_arr=list()
	#non-negative numbers only please!
	l=max(0,l)
	for i in range(l):
		empty_dig_arr.append(d)
	return empty_dig_arr




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
	dna=dna+extraNs
	coding_dna = Seq(dna, IUPAC.ambiguous_dna)
	messenger_rna = coding_dna.transcribe()
	translation=messenger_rna.translate()
	return translation






def isIntegral(s):
	ire=re.compile(r'^\s*(\d+)\s*$')
	sr=re.search(ire,s)
	if(sr):
		return True
	else:
		return False




def fileLineCount(f):
	reader=open(f,'r')
	line_count=0
	for line in reader:
		line_count+=1
	reader.close()
	return line_count







#return TRUE if the first argument (a)
# is a substring of the second (b)
def a_subseq_of_b(a,b):
	try:
		a_index=b.index(a)
		return True
	except ValueError:
		return False








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
	return date_result


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
		if(a_i in b):
			#this is the intersection
			print a_i,"\t\t\t",a_i
		else:
			#this is what is in a but NOT in b
			print a_i,"\t<\t\t",blank
	for b_i in b:
		if(not(b_i in a)):
			#this is what is in b, but NOT in a
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
