#!/usr/bin/env python

#code perhaps to later be used to combine igblast and repchar

def makeIGBlastCmd(igdata,fasta,output):
	cmd=igdata+"/bin/igblastn "
	cmd+=" -num_threads 6 "
	#todo allow domain change
	cmd+=" -domain_system imgt "
	cmd+=" -query "+str(fasta)+" "
	#todo allow change database
	cmd+=" -germline_db_V /home/data/DATABASE/01_22_2014/human/ReferenceDirectorySet/human_IG_V.fna -germline_db_D /home/data/DATABASE/01_22_2014/human/ReferenceDirectorySet/human_IG_D.fna -germline_db_J /home/data/DATABASE/01_22_2014/human/ReferenceDirectorySet/human_IG_J.fna "
	#todo allow seq_type change
	cmd+=" -ig_seqtype Ig "
	cmd+=" -auxiliary_data "+igdata+"/optional_file/human_gl.aux "
	cmd+=" -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop' "
	cmd+=" -o "+output+" "
	return cmd

	

def make_rep_char_cmd(igbout,fasta):
	cmd="python rep_char.py "
	cmd+=" -json_out "+igbout+".json "
	cmd+=" -cdr3_hist_out "+igbout+".cdr3_hist.tsv "
	cmd+=" "+igbout+" "
	cmd+=" "+igbout+".vdjml "
	cmd+=" "+igbout+".rc_out.tsv " 
	cmd+=" /home/data/DATABASE/01_22_2014 "
	cmd+=" "+fasta
	return cmd




#program to merge CDR3 histograms for kabat and imgt modes
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='igblast/rep_char pipeline')
	args=parser.parse_args()
	igdata="/usr/local/igblast-1.3.0/"
	fasta=

