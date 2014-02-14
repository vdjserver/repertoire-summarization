#!/usr/bin/env python

import os
import re
import glob
from os.path import basename
from segment_utils import getFastaListOfDescs,getQueryIndexGivenSubjectIndexAndAlignment,getAdjustedCDR3StartFromRefDirSetAllele,getADJCDR3EndFromJAllele
from imgt_utils import imgt_db
from utils import read_fasta_file_into_map,biopythonTranslate,rev_comp_dna


def getCDR3Length(VLine,JLine):
	return 1


#test

#define inputs for run
vfasta="/home/esalina2/round1_imgt/human_IG_V.fna"
dfasta="/home/esalina2/round1_imgt/human_IG_D.fna"
jfasta="/home/esalina2/round1_imgt/human_IG_J.fna"
fasta_file_list=[vfasta,dfasta,jfasta]
blast_out="/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna.igblastn.imgt.out"

#define inputs for database
base_dir="/home/data/DATABASE/01_22_2014/"
organism="human"

#create DB OBJ
imgtdb_obj=imgt_db(base_dir)
#getAdjustedCDR3StartFromRefDirSetAllele(allele,imgtdb_obj,organism="human"):

vlist=getFastaListOfDescs(fasta_file_list[0])
dlist=getFastaListOfDescs(fasta_file_list[1])
jlist=getFastaListOfDescs(fasta_file_list[2])

blast_data=open(blast_out,'r')
line_num=0
counts_map=dict()
capture_flag=False
is_heavy=True
VHitLine=None
JHitLine=None
currentV=None
currentJ=None
currentD=None

trans=0
tot=0
nona=0
max_len=0
cdr3_hist=dict()
query_seq_map=read_fasta_file_into_map("/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna")
domain_list=["imgt","kabat"]
rsmap=dict()
remap=dict()
for d in domain_list:
	rsmap[d]=dict()
	remap[d]=dict()
	cdr3_hist[d]=dict()
	for i in range(0,100):
		cdr3_hist[d][i]=0
v_alleles_review=dict()




vqmap=dict()


files=glob.glob("/home/esalina2/Downloads/HSEQS/IMGT_HighV-QUEST_individual_files_folder/*")
files=[]
for f in files:
	reader=open(f,'r')
	num=0
	fbase=basename(f)
	fbpieces=fbase.split("_")
	readname=str(fbpieces[0].strip())
	vqmap[readname]=dict()
	for line in reader:
		segs=['V','D','J']
		for seg in segs:
			if(not(seg in vqmap[readname])):
				vqmap[readname][seg]=None
			if(line.startswith(seg+"-GENE")):
				pieces=line.split(";")
				if(len(pieces)>=2):
					info=pieces[1]
					segre=re.compile(r'homsap\s+(\S+)\s+',re.IGNORECASE)
					mr=re.search(segre,info)
					if(mr):
						seg_res=mr.group(1)
						#print "readname is ",readname
						vqmap[readname][seg]=seg_res
						print "\n\n\n*************************"
						print "readname="+readname
						print "seg="+seg
						print "storage="+seg_res
						print "*******************************\n\n\n"
							
					else:
						print "FAIL re MATCH ON ",readname," : ",info," for seg=",seg
		if(num>=40):
			break;
		num+=1
	reader.close()


def getNumDashInStr(s):
	orig_len=len(s)
	s_no_dash=re.sub(r'\-','',s)
	nodash_len=len(s_no_dash)
	return orig_len-nodash_len


def repStr(s,n):
	if(n<=0):
		return ""
	else:
		return s+repStr(s,n-1)



def annotatedCDR3(s_aln,q_aln,s_c,q_c,s_s,q_s):
	#get CDR3 starts relative to the alignment starts
	aln_rel_s=s_c-s_s
	aln_rel_q=q_c-q_s
	print "subject from is ",s_s
	print "subject cdr3 is ",s_c
	print "aln rel sub is ",aln_rel_s
	print "query from is ",q_s
	print "query CDR3 is ",q_c
	print "aln rel q is ",aln_rel_q
	#sys.exit(0)


	#find position including gaps in the alignment for the SUBJECT STAR
	#all this code is 1-based indices
	#except for actually doing string comparison and string-printing which is 0-based
	s_temp=0
	if(s_aln[0]=="-"):
		s_pos=0
	else:
		s_pos=1
	while(s_pos<=aln_rel_s):
		if(s_aln[s_temp]!="-"):
			s_pos+=1
		s_temp+=1
	s_actual=s_temp
	#0-based indices!
	q_temp=0
	if(q_aln[0]=="-"):
		q_pos=0
	else:
		q_pos=1
	while(q_pos<=aln_rel_q):
		if(q_aln[q_temp]!="-"):
			q_pos+=1
		q_temp+=1
	q_actual=q_temp 
	#0-based indices
	#if(q_actual!=s_actual):
	#	print "***WARNING***"
	annLines=["","","",""]
	annLines[0]=repStr(" ",s_actual)
	annLines[0]+="X"
	annLines[1]=s_aln
	annLines[2]=q_aln
	annLines[3]=repStr(" ",q_actual)
	annLines[3]+="X"
	annotated=""
	for i in range(len(annLines)):
		#print annLines[i]
		annotated+=str(i)+" : "+annLines[i]
		if(i<len(annLines)-1):
			annotated+="\n"
	return annotated
	

	





def CDR3ANAL(currentV,currentJ,vHitLine,jHitLine,currentD,currentQueryName):
	global transbiopythonTranslate
	global v_alleles_review
	global trans
	global nona
	reason=None
	currentQueryName=str(currentQueryName.strip())
	vpieces=VHitLine.split("\t")
	jpieces=JHitLine.split("\t")
	if(vpieces[6]==currentV and jpieces[6]==currentJ and (not(currentJ==None)) and (not(currentV==None))   ):
		#print "WE'RE IN BUSINESS!"
		if(not(currentV in rsmap)):
			V_CDR3_START=getAdjustedCDR3StartFromRefDirSetAllele(currentV,imgtdb_obj,"human")
			rsmap[currentV]=V_CDR3_START
		else:
			V_CDR3_START=rsmap[currentV]
		print "The adjusted CDR3 start is ",V_CDR3_START
		if(V_CDR3_START!=(-1)):
			print "c:v,d,j =",currentV,",",currentD,",",currentJ
			#print VHitLine
			if(not(currentJ==None)):
				#got both
				#vdat=imgtdb_obj.getIMGTDatGivenAllele(currentV,"human")
				#print vdat+"\n\n"
				#jdat=imgtdb_obj.getIMGTDatGivenAllele(currentJ,"human")
				#print jdat+"\n\n"
				#ddat=imgtdb_obj.getIMGTDatGivenAllele(currentD,"human")
				#print ddat+"\n\n\n"
				vq_aln=vpieces[18]
				vs_aln=vpieces[19]
				vq_f=int(vpieces[14])
				vq_t=int(vpieces[15])
				#print "qf,qt=",vq_f,",",vq_t
				vs_f=int(vpieces[16])
				vs_t=int(vpieces[17])
				#print "sf,st=",vs_f,",",vs_t
				query_cdr3_start=getQueryIndexGivenSubjectIndexAndAlignment(vq_aln,vs_aln,vq_f,vq_t,vs_f,vs_t,V_CDR3_START)
				print "QUERY cdr3 start is ",query_cdr3_start
				#sys.exit(0)
				if(not(query_cdr3_start==(-1))):
					annotated=annotatedCDR3(vs_aln,vq_aln,V_CDR3_START,query_cdr3_start,vs_f,vq_f)
					#annotated=annotatedCDR3(s_aln ,q_aln,s_c,q_c,s_s,q_s):
					print "s_to=",vs_t
					print "q_to=",vq_t
					print "The annotated CDR3 area (queryname="+currentQueryName+") is :\n"
					print annotated+"\n"
					jq_aln=jpieces[18]
					js_aln=jpieces[19]
					jq_f=int(jpieces[14])
					jq_t=int(jpieces[15])
					js_f=int(jpieces[16])
					js_t=int(jpieces[17])
					if(not(currentJ in remap)):
						J_CDR3_END=getADJCDR3EndFromJAllele(currentJ,imgtdb_obj,"human")
						remap[currentJ]=J_CDR3_END
					else:
						J_CDR3_END=remap[currentJ]
					print "ADJ. CDR3 END IN J : ",J_CDR3_END
					if(J_CDR3_END==(-1)):
						reason="REASON: J RECORD HAS NO APPARENT CDR3_END in imgt.dat"
						print reason
					else:
						query_cdr3_end=getQueryIndexGivenSubjectIndexAndAlignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t,J_CDR3_END)
						print "QUERY cdr3 end is ",query_cdr3_end
						if(query_cdr3_end==(-1) and js_f>J_CDR3_END):
							reason="REASON: JCDR3 END IS BEFORE THE ALIGNMENT STARTS"
							print reason
						if(query_cdr3_end==(-1) and js_t<J_CDR3_END):
							reason="REASON: JCDR3 END IS AFTER THE ALIGNMENT ENDS"
							print reason
						if(query_cdr3_end!=(-1) and query_cdr3_start!=(-1)):
							if(query_cdr3_end<=query_cdr3_start):
								raise Exception("why is start>=end?????")
							print "I want to translate!"
							query_coding_seq=query_seq_map[currentQueryName]
							if(vHitLine.find("reversed|"+currentQueryName)!=(-1)):
								query_coding_seq=rev_comp_dna(query_coding_seq)
							coding_seq=query_coding_seq[(query_cdr3_start-1):(query_cdr3_end-1)]
							if(coding_seq.find("N")==(-1)):
								translation=biopythonTranslate(coding_seq)
								print "the coding seq  is : ",coding_seq
								print "The translation is : ",translation
							trans+=1
							#annotatedCDR3(s_aln,q_aln,s_c,q_c,s_s,q_s):
							jAnnotated=annotatedCDR3(js_aln,jq_aln,J_CDR3_END,query_cdr3_end,js_f,jq_f)
							print "THE_J_ANN IS"
							print jAnnotated
							cdr3_len=int((query_cdr3_end-query_cdr3_start)/3.0)
							print "CDR3_LEN="+str(cdr3_len)
							cdr3_hist["imgt"][cdr3_len]+=1
							#if(currentQueryName in vqmap):
							#	if(vqmap[currentQueryName]['V']==currentV and vqmap[currentQueryName]['J']==currentJ and currentD==vqmap[currentQueryName]['D']):
							#		print "THIS IS A SUPER_MATCH!"
						else:
							reason="REASON: J END NOT MAPPABLE VIA ALIGNMENT ; USE JUNCTION END"
							print reason
							#query_cdr3_end=jq_f
							#take CDR3 to be the junction end
							#query_coding_seq=query_seq_map[currentQueryName]
							#query_coding_seq=query_seq_map[currentQueryName]
							#if(vHitLine.find("reversed|"+currentQueryName)!=(-1)):
							#	query_coding_seq=rev_comp_dna(query_coding_seq)
							#if(query_cdr3_end-3-query_cdr3_start>=3):
							#	coding_seq=query_coding_seq[(query_cdr3_start-1):(query_cdr3_end-1)]
							#	trans+=1
							#	cdr3_len=int(len(coding_seq)/3)
							#	cdr3_hist[cdr3_len]+=1
							#	if(coding_seq.find("N")==(-1)):
							#		translation=biopythonTranslate(coding_seq)
							#		print "The JUNCTION coding seq is :",coding_seq
							#		print "The JUNCTION translation is : ",translation
				else:
					reason="REASON: CDR3 START NOT MAPPABLE VIA ALIGNMENT"
					print reason
					if(V_CDR3_START>vs_t):
						reason="REASON: THE CDR3 START IS AFTER THE ALIGNMENT ENDS"
						print reason
						reason="REASON: THE CDR3 START IS AFTER THE ALIGNMENT ENDS and D,J are :"+str([currentD,currentJ])
						print reason
		else:
			#V_CDR3_START is (-1)
			print "V_CDR3_START is negative 1"
			print "c:v,d,j =",currentV,",",currentD,",",currentJ
			reason="REASON: V RECORD "+currentV+" HAS NO APPARENT CDR3_START in imgt.dat"
			print reason
			v_alleles_review[currentV]=imgtdb_obj.getIMGTDatGivenAllele(currentV)
	else:
		reason="REASON: currentJ, currentV, vline or jline None!"
		nona+=1
		print reason
	if(currentQueryName in vqmap):
		if(vqmap[currentQueryName]['V']==currentV and vqmap[currentQueryName]['J']==currentJ and currentD==vqmap[currentQueryName]['D']):
			print "FULL VQUEST MATCH ON ",currentQueryName," V=",currentV," J=",currentJ
		elif(vqmap[currentQueryName]['V']==currentV):
			print "VQUEST MATCH ONLY ON ",currentQueryName," for V=",currentV
		elif(vqmap[currentQueryName]['J']==currentJ):
			print "VQUEST MATCH ONLY ON ",currentQueryName," for J=",currentJ
		else:
			print "VQUEST mismatch on ",currentQueryName," V=",currentV," VQV=",vqmap[currentQueryName]['V']," J=",currentJ," VQJ=",vqmap[currentQueryName]['J']
	if(reason):
		print">"+currentQueryName
		print query_seq_map[currentQueryName]

currentQuery=None
for line in blast_data:
	if(line.startswith("# Query")):
		tot+=1
		#print "At a query : ",line
		if(VHitLine!=None and JHitLine!=None and currentV!=None and currentJ!=None):
			CDR3ANAL(currentV,currentJ,VHitLine,JHitLine,currentD,currentQuery)
		else:
			nona+=1
		#reset values
		VHitLine=None
		JHitLine=None
		currentV=None
		currentJ=None
		currentD=None
		# Query: G7Q0IVE01ADLSM
		qr=re.compile(':\s+([A-Z0-9]+)\s*$')
		mr=re.search(qr,line)
		if(mr):
			currentQuery=mr.group(1)
		else:
			print "fail qre!"
			sys.exit(0)
		print "\n\n\n"
	elif(re.match(r'^[VJ]\s+',line)):
		line_pieces=line.split("\t")
		if(len(line_pieces)>=34):
			#got data!
			if(VHitLine==None and line.startswith("V\t")):
				VHitLine=line
			elif(JHitLine==None and line.startswith("J\t")):
				JHitLine=line
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
				#segment=getSegmentName(pieces[i],list_to_use)
				segment=pieces[i]
				#print "GOT A SEGMENT "+segment+" heavy_status="+str(is_heavy)
				if(list_to_use==vlist):
					currentV=segment
					#print "saved as V",currentV
				elif(list_to_use==jlist):
					currentJ=segment
					#print "saved as J",currentJ
				elif(list_to_use==dlist):
					currentD=segment
		capture_flag=False
	if(line.startswith("# V(D)J rearrangement")):
		#print line
		if(line.find('Top D')!=(-1)):
			is_heavy=True
		else:
			is_heavy=False
		capture_flag=True
	line_num+=1
	#print "at line="+str(line_num)




print "tot ",tot
print "trx ",trans
diff=(tot-trans)
print "dif ",diff
print "nna ",nona
print "tr/t",float(trans)/float(tot)

for d in domain_list:
	print "HISTOGRAM",d
	for i in range(0,100):
		print "LEN=",i,"count=",cdr3_hist[d][i]


#print "***********************************************************"
#print "REVIEWING ALLELELS "
#for a in v_alleles_review:
#	print "********************************************"
#	print "REVIEW ",a
#	print "REVIEW DESCRIPTOR ",imgtdb_obj.extractDescriptorLine(a)
#	print "********************************************"
#	print v_alleles_review[a]
#	print "\n\n\n\n"		

