#!/usr/bin/env python

import os
from subprocess import call
import binascii
from Bio.Blast import NCBIXML
from Bio import SeqIO
import shutil
import re
from igblast_utils import *
from utils import *

#SOME ROUTINES USED IN IMGT/IGBLAST mapping


def isIMGTNameInAllKeyList(imgt_n,s_map):
	imgt_list=list()
	for key in s_map:
		imgt_key=extractIMGTNameFromKey(key)
		imgt_list.append(imgt_key)
	if(imgt_n in imgt_list):
		return True
	else:
		return False



def calc_pct_id(q,s):
	nid=0
	for i in range(min(len(q),len(s))):
		c_q=q[i]
		c_s=s[i]
		if(c_q==c_s):
			nid+=1
	pct_id=float(nid)/float(min(len(q),len(s)))
	return pct_id



def extractIMGTNameFromKey(k):
	pieces=k.split("|")
	#print "the pieces are :",pieces
	#sys.exit(0)
	return pieces[1]





def make_q_conformAsSubsequences(q,s):
	if(len(q)>len(s)):
		return [q,0,0]
		#new_sub=make_q_conformAsSubsequences(s,q)
		#new_q=make_q_conformAsSubsequences(q,new_sub[0])
		#return new_q
	if(len(q)==len(s)):
		return [q,0,0]
	try:
		i=s.index(q)
		numHead=0
		numTail=0
		if(i>0):
			q=s[0:i]+q
			numHead=i
		if(len(q)<len(s)):
			q=q+s[len(q):]
			numTail=len(s[len(q):])
		return [q,numHead,numTail]
	except ValueError:
		return [q,0,0]
	



def make_q_HeadOrTailNMatchSubjectOnSameLengthData(q,s):
	numNHead=0
	numNTail=0
	if(q.startswith("N")):
		if(re.compile("^(N+)[^N]").match(q)):
			rem=re.search("^(N+)[^N]",q)
			if(rem):
				NText=rem.group(1)
				numNHead=len(NText)
				#print "q before:",q
				q=q[len(NText):]
				q=s[0:len(NText)]+q
				#print "q  after:",q		
	if(q.endswith("N")):
		#print "tests as ENDED WITH N"
		num_sub=0
		#print "q before:",q
		while(q.endswith("N")):
			q=q[0:len(q)-1]
			num_sub+=1
		numNTail=num_sub
		q=q+s[len(s)-num_sub:]
		#print "q after :",q		
	return [q,numNHead,q,numNTail]

def get_key_from_blast_title(title):
	space_index=title.find(' ')
	space_index+=1
	#print "from title=",title," the space index is ",space_index
	key=title[space_index:]
	return str(key)
	
	
def hold():
	return

def getAnalysisLine(q_fasta_seq,q_fasta_name,blast_out_path,ref_fasta_path,ref_fasta_map):
	result_handle = open(blast_out_path)
	query_length=len(q_fasta_seq)
	blast_records = NCBIXML.parse(result_handle)
	analysis_line=list()
	descs=list()
	for blast_record in blast_records:
		if(len(blast_record.alignments)!=1):
			print "TOO MANY ALIGNMENTS",blast_out_path
			sys.exit(1)
		alignment=blast_record.alignments[0]
		#print "An alignment length is ",alignment.length
		if(len(alignment.hsps)!=1):
			print "TOO MANY HSPs",blast_out_path
			sys.exit(1)
		hsp=alignment.hsps[0]
		#print "alignment title:",alignment.title
		s_map_key=get_key_from_blast_title(alignment.title)
		#print "retrieved key:",s_map_key
		#print "query name   :",q_fasta_name
		subject_seq=ref_fasta_map[s_map_key].upper()
		#print "from map:",subject_seq
		#print "QUERY=",q_fasta_seq
		#print "SBJCT=",subject_seq
		imgt_name=extractIMGTNameFromKey(s_map_key)
		analysis_line.append(q_fasta_name);descs.append("q_name")
		analysis_line.append(isIMGTNameInAllKeyList(q_fasta_name,ref_fasta_map));descs.append("q_name_in_ref_db")
		analysis_line.append(imgt_name==q_fasta_name);descs.append("names_eq?")
		other_list=list()
		for key in ref_fasta_map:
			key_seq=ref_fasta_map[key]
			if(a_subseq_of_b(subject_seq,key_seq)):
				other_list.append(extractIMGTNameFromKey(key))
		sep_col=" : "
		analysis_line.append(sep_col.join(other_list));descs.append("subject_is_a_subseq_of_these")
		other_list=list()
		for key in ref_fasta_map:
			key_seq=ref_fasta_map[key]
			if(a_subseq_of_b(key_seq,subject_seq)):
				other_list.append(extractIMGTNameFromKey(key))
		analysis_line.append(sep_col.join(other_list));descs.append("these_are_subseqs_of_sbjct")
		other_list=list()
		for key in ref_fasta_map:
			key_seq=ref_fasta_map[key]
			if(subject_seq==key_seq):
				other_list.append(extractIMGTNameFromKey(key))
		analysis_line.append(sep_col.join(other_list));descs.append("subject_is_SAME_as_these")
		analysis_line.append(len(q_fasta_seq));descs.append("q_len")
		analysis_line.append(imgt_name);descs.append("imgt_subject_name_extraction")
		analysis_line.append(s_map_key);descs.append("full_ref_name")
		q_stat=""
		if(len(q_fasta_seq)<len(subject_seq)):
			q_stat="shorter"
		elif(len(q_fasta_seq)==len(subject_seq)):
			q_stat="same"
		else:
			q_stat="longer"
		analysis_line.append("query is "+q_stat+" compared to subject");descs.append("len_comp_note")
		uncond_match=False
		cond_n_match=False
		cond_n_sub_match=False
		note=""
		if(((subject_seq==q_fasta_seq) or (1.0==calc_pct_id(q_fasta_seq,subject_seq) )) or (imgt_name==q_fasta_name)  ):
			#print "PERFECT UNCONDITIONAL MATCH"
			note="it matches unconditionally"
			analysis_line.append(calc_pct_id(q_fasta_seq,subject_seq));descs.append("pct_id")
			uncond_match=True
			analysis_line.append("N/A");descs.append("pct_id_nMod_qsadd")
			hold()
		else:
			query_mod=make_q_HeadOrTailNMatchSubjectOnSameLengthData(q_fasta_seq,subject_seq)
			analysis_line.append(calc_pct_id(query_mod[0],subject_seq));descs.append("pct_id_Nmod")
			if(query_mod[0]==subject_seq):
				#print "PERFECT MATCH AFTER N-REPLACE ON HEAD AND TAIL"
				cond_n_match=True
				note+="matches after replacing query head/tail N with subject bases when query is a substring of the data"
			else:
				query_mod_mod=make_q_conformAsSubsequences(query_mod[0],subject_seq)
				analysis_line.append(calc_pct_id(query_mod_mod[0],subject_seq));descs.append("pct_id_nMod_qsadd")
				if(query_mod_mod[0]==subject_seq):
					#print "PERFECT MATCH AFTER N-REPLACE AND SUBSEQING"
					cond_n_sub_match=True
					note+="matches after both replacing head/tail Ns and adding sequence on head and tail if the query is a subsequence"
				else:
					#print "\n\n"
					note+="DOESN'T match despite 'rescue' efforts of N-beheading/betailing and subject_data reconstruction"
					#if(q_fasta_name=="IGHV2-5*07"):
					#	print "plain query : ",q_fasta_seq
					#	print "subject     : ",subject_seq
					#	print "first mod   : ",query_mod[0]
					#	print "subject     : ",subject_seq
					#	print "second mod  : ",query_mod_mod[0]
					#	print "subject     : ",subject_seq
					#	system.exit(0)
					#print "retrieved key:",s_map_key
					#print "query name   :",q_fasta_name
					#print "query length :",len(q_fasta_seq)
					#print "Query name in search DB ?:",isIMGTNameInAllKeyList(q_fasta_name,ref_fasta_map)				
					#print "IMPERFECT MATCH AFTER N-REPLACE AND SUBSEQING"
					#print hsp.query
					#print hsp.match
					#print hsp.sbjct
					#print "pct id orig:",calc_pct_id(q_fasta_seq,subject_seq)
					#print "pct id mod:",calc_pct_id(query_mod[0],subject_seq)
					#print "pct id mod2:",calc_pct_id(query_mod_mod[0],subject_seq)
		analysis_line.append(note);descs.append("note")
		
		#descs=["q_name","exists_in_ref_DB","names_of_q_and_s_same?","q_len","imgt_extraction","sub_full_name","Len_comp_status","note"]
		#print "Analysis line :",analysis_line
		return [descs,analysis_line]




			
q_fasta="/usr/local/igblast_from_lonestar/database/human_gl_V.fna"
imgt_fasta="/home/esalina2/Downloads/imgt.2/www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP"

q_fasta_data=read_fasta(q_fasta)
s_fasta_data=read_fasta(imgt_fasta)
s_fasta_map=read_fasta_into_map(s_fasta_data)



rand_dir="/tmp/tmp_imgt_"+binascii.hexlify(os.urandom(16))
print "rand dir : ",rand_dir
os.mkdir(rand_dir)
i=0
while (i<len(q_fasta_data)):
	q_name=q_fasta_data[i]
	q_seq=q_fasta_data[i+1]
	o_fasta_file=rand_dir+"/"+str(i)+".fna"
	o_fasta=open(o_fasta_file,'w')
	o_fasta.write(str(">")+str(q_name)+"\n"+str(q_seq)+"\n")
	o_fasta.close()
	tmp_imgt=rand_dir+"/imgt.fna"
	if(i==0):
		o_fasta=open(tmp_imgt,'w')
		j=0
		while j<len(s_fasta_data):
			o_fasta.write(">"+s_fasta_data[j]+"\n")
			o_fasta.write(s_fasta_data[j+1]+"\n")
			j+=2
		o_fasta.close()
		format_db_cmd="/usr/local/bin/makeblastdb -dbtype nucl -in "+tmp_imgt
		print "The format cmd is ",format_db_cmd
		call(format_db_cmd.split(" "))
	blast_output_path=o_fasta_file+".blast.out"
	blast_cmd="/usr/local/bin/blastn -query "+o_fasta_file+" -outfmt 5 -num_alignments 1 -num_descriptions 1 -max_target_seqs 1 -db "+tmp_imgt+" -out "+blast_output_path
	call(blast_cmd.split(" "))
	ret_data=getAnalysisLine(q_seq.upper(),q_name,blast_output_path,tmp_imgt,s_fasta_map)
	al=ret_data[1]
	#print "len al is ",len(al)
	#print "al is ",al
	descs=ret_data[0]
	#print "len descs is ",len(descs)
	sep="\t"
	if(i==0):
		print sep.join(descs)
	print sep.join(map(str,al))
	#print "i=",i
	i+=2




