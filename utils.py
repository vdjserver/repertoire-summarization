#!/usr/bin/env python
from datetime import datetime, date, time
import pprint, pickle
import os
import subprocess
import urllib2
import re
from subprocess import call
import nwalign as nw
import pickle





def getRevCompInterval(i,seq_len_in_bp):
	#form an interval in one strand
	#find the interval in the reverse complement strand
	i_s=i[0]
	i_e=i[1]
	int_len=abs(i_s-i_e)
	ne_s=seq_len_in_bp-i_e+1
	ne_e=ne_s+int_len
	return [ne_s,ne_e]



#given a seq record, an interval into it, and an inverted flag (telling wether interval is in the opposite strand or no
#return the subsequece (inclusive) indicated by the interval
def extractSubSeq(interval,seq_rec,is_inverted):
	if(is_inverted):
		rc_inteval=getRevCompInterval(interval,len(str(seq_rec.seq)))
		return extractSubSeq(interval,seq_rec,False)
	else:
		s=str(seq_rec.seq)
		s=s[interval[0]-1,interval[1]
		






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

def makeEmptyArrayOfStringsOfLen(l):
	empty_str_arr=list()
	#non-negative numbers only please!
	l=max(0,l)
	for i in range(l):
		empty_str_arr.append("")
	return empty_str_arr





def makeAllMapValuesVal(m,v):
	if(m is not None):
		for k in m:
			m[k]=v
	return m
		



def rev_dna(dna):
	reversed_dna=""
	for ub in range(len(dna)):
		b=len(dna)-ub-1
		reversed_dna+=dna[b]
	return reversed_dna
		

def comp_dna(dna):
	comped=""
	comp=dict()
	comp['A']='T'
	comp['T']='A'
	comp['G']='C'
	comp['C']='G'
	for b in range(len(dna)):
		if(dna[b] in comp):
			comped+=comp[dna[b]]
		else:
			comped+=dna[b]
	return comped
			

def rev_comp_dna(dna):
	return comp_dna(rev_dna(dna))






def biopythonTranslate(dna):
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC
	#coding_dna = Seq(dna, IUPAC.unambiguous_dna)
	if(len(dna)>0 and len(dna)<3):
		return "X"
	elif(len(dna)==0):
		return ""
	mod_3=((len(dna))%3)
	#print "mod3=",mod_3
	if(mod_3!=0):
		dna_good=dna[0:(len(dna)-mod_3)]
		dna=dna_good
		return biopythonTranslate(dna)
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









def counter_map_inc(counter_map,key):
	if key in counter_map:
		counter_map[key]+=1
	else:
		counter_map[key]=1
	return counter_map


def doesStringContainAnApostrophe(s):
	if(re.search("'",s)):
		return True
	else:
		return False



def a_subseq_of_b(a,b):
	try:
		a_index=b.index(a)
		return True
	except ValueError:
		return False


def isAlignmentFreeOfBothMutationsAndIndels(alignment,TryRemovalOfHeadAndTailDash=False):
	if(alignment[0]==alignment[1]):
		return True
	else:
		if(TryRemovalOfHeadAndTailDash and (removeHeadAndTailMultiDash(alignment[0])==removeHeadAndTailMultiDash(alignment[1]))):
			return True
		return False


def readMultipleMapFilesIntoSingleMapWithGlob(g):
	map_files=glob.glob(g)
	total_map=dict()
	for map_file in map_files:
		single_map=read_map_from_file(map_file)
		for key in single_map:
			total_map[key]=single_map[key]
	return total_map



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


def counter_map_get(counter_map,key):
	if key in counter_map:
		return counter_map[key]
	else:
		return 0




def touch(fname, times=None):
	try:
		with file(fname, 'a'):
			os.utime(fname, times)
	except Exception:
		print "Error in touching file",fname


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



def read_fasta_file_into_map(fasta_path,alwaysSeqToUpper=True):
	fasta_data=read_fasta(fasta_path,alwaysSeqToUpper)
	fasta_map=read_fasta_into_map(fasta_data)
	return fasta_map


def read_fasta_string(fasta_string,alwaysSeqToUpper=True):
	lines=fasta_string.split("\n")
	data=list()
	for line in lines:
		temp_line=line.strip()
		if(temp_line.startswith(">")):
			data.append(temp_line[1:])
			data.append("")
		else:
			if(alwaysSeqToUpper):
				data[len(data)-1]+=temp_line.upper()
			else:
				data[len(data)-1]+=temp_line	
	return data






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





def removeHeadAndTailMultiDash(s):
	s=re.sub(r'\-+$',"",s)
	s=re.sub(r'^\-+',"",s)
	return s



#little utility to return the modes
def get_domain_modes():
	domain_list=["imgt","kabat"]
	#domain_list=["kabat"]
	return domain_list





def removeTerminatingSemicolonIfItExists(s):
	s=re.sub(r';+$',"",s)
	return s
	




def format_blast_db(fastaPath):
	format_db_cmd="/usr/local/bin/makeblastdb -dbtype nucl -in "+fastaPath
	print "The format cmd is ",format_db_cmd
	try:
		call(format_db_cmd.split(" "))
	except Exception, e:
		print "BLAST_FORMAT_CMD error({0})"+format( e.strerror)
	return format_db_cmd



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



def read_fasta_into_map(fasta_data):
	i=0
	fasta_map=dict()
	while (i<len(fasta_data)):
		#print "i=",i
		fasta_map[fasta_data[i]]=fasta_data[i+1]
		i+=2
	return fasta_map



def convertToDefaultStringDict(m,s):
	dd=defaultdict(s)
	for k in m:
		dd[k]=m[k]
	return dd




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




def merge_maps(ma,mb):
	mc=dict()
	for a in ma:
		mc[a]=ma[a]
	for b in mb:
		mc[b]=mb[b]
	return mc




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
	


def getNumberLinesInFile(p):
	reader=open(p,'r')
	lines=0
	for line in reader:
		lines+=1
	return lines


def getLengthOfListInFile(p):
	mylist=read_list_from_file(p)
	return len(mylist)



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



def read_list_from_file(path):
	reader=open(path,'r')
	file_list=list()
	for line in reader:
		line=line.strip()
		file_list.append(line)
	reader.close()	
	return file_list


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


def returnTimeStamp():
	date_cmd="date"
	date_cmd_args="+%m_%d_%Y_%H_%M_%S.%N"
	date_result=subprocess.check_output([date_cmd,date_cmd_args])



def printList(l):
	for i in range(len(l)):
		print i,":",l[i]




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


def readFileIntoString(path):
	f=open(path,'r')
	data=f.read()
	f.close()
	return data



def readAURL(url):
	f = urllib2.urlopen(url)
	html=f.read()
	return html


def writeStringToFilePathAssumingDirectoryExists(string,file_path):
	outfile=open(file_path,'w')
	outfile.write(string)
	outfile.close()




def downloadURLToLocalFileAssumingDirectoryExists(url,local):
	url_content=readAURL(url)
	writeStringToFilePathAssumingDirectoryExists(url_content,local)



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
