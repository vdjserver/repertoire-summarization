"""
Repertoire-summarization
"""

import sys
import vdjml
import utils
from .version import __version__
import argparse
import json
import defaults
from imgt_utils import imgt_db
from Bio import SeqIO
import summarize
from cdr3_hist import histoMapClass
from segment_utils import IncrementMapWrapper, recombFreqManager

def makeSummaryParserArgs():
	"""Command line arguments for repsum"""
	parser = argparse.ArgumentParser();
	parser.description='Summarize VDJ assignment in VDJML for immune repertoire sequencing data. VERSION '+__version__
	parser.add_argument('vdjml_in',type=str,nargs=1,help="the path to the input VDJML file")
	parser.add_argument('char_out',type=str,nargs=1,help="the path to the output TSV file of read-level repertoire characterization data")
	parser.add_argument('vdj_db_root',type=str,nargs=1,help="path to the VDJ directory root REQUIRED")
	parser.add_argument('qry_fasta',type=str,nargs=1,help="path to the input fasta file of query (Rep-Seq) data input to IgBLAST")
	parser.add_argument('db_organism',type=str,nargs=1,default="human",help="the organism IgBLASTed against;must exist under vdj_db_root",choices=["human","Mus_musculus"])
	return parser

def main():
	"""Generate a repertoire summarization"""
	parser = makeSummaryParserArgs()
	args = parser.parse_args()

	if (not args):
		args.print_help()
		sys.exit()

	organism = utils.extractAsItemOrFirstFromList(args.db_organism)
	vdjml_file = utils.extractAsItemOrFirstFromList(args.vdjml_in)
	vr = vdjml.Vdjml_reader(vdjml_file)
	meta = vr.meta()
	imgtdb_obj = imgt_db(utils.extractAsItemOrFirstFromList(args.vdj_db_root))
	#print "a fact is ",fact
	#print "the file is ",args.igblast_in[0]
	query_fasta = utils.extractAsItemOrFirstFromList(args.qry_fasta)
	tot_reads=None
	fasta_reader=SeqIO.parse(open(query_fasta, "r"), "fasta")
	modes=['kabat','imgt']
	my_cdr3_map=histoMapClass(modes)
	segment_counter=IncrementMapWrapper()
	combo_counter=recombFreqManager()
	read_num=1
	rep_char_out=utils.extractAsItemOrFirstFromList(args.char_out)
	fHandle=open(rep_char_out,'w')
	if(False):
			logFile=rep_char_out+".log"
			logHandle=open(logFile,'w')
	else:
			logHandle=None
	every_read=100

	# step through VDJML file
	while vr.has_result():
			read_result_obj = vr.result()

			if(read_num>1 and read_num%every_read==0):
					print "Processing read",read_num,"..."
			elif(read_num==1):
					print "Processing reads..."
			query_record=fasta_reader.next()

			#analyze a read's results
			read_analysis_results=summarize.rep_char_read(read_result_obj,meta,organism,imgtdb_obj,query_record,False)

			#handle cdr3 length/histogram
			if('cdr3_length_results' in read_analysis_results):
					cdr3_res=read_analysis_results['cdr3_length_results']
					for mode in modes:
							my_cdr3_map.inc(mode,cdr3_res[mode])

			#handle segment counting
			if('VDJ' in read_analysis_results):
					segments=read_analysis_results['VDJ']
					for s in segments:
							actual=segments[s]
							if(actual is not None):
									segment_counter.increment(actual)
					vseg=segments['V']
					dseg=segments['D']
					jseg=segments['J']
					combo_counter.addVDJRecombination(vseg,dseg,jseg,False)

			#write rep-char
			summarize.appendAnnToFileWithMap(fHandle,read_analysis_results['ann_map'],read_num,query_record.id,None,"None",logHandle)

			#increment the read number
			read_num+=1
			vr.next()

	print "Processed total of",str(read_num-1),"reads!"

	#close rep_char out
	fHandle.close()

	#close log
	if(logHandle is not None):
			logHandle.close()

