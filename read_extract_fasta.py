#!/usr/bin/env python

from Bio import SeqIO
from utils import extractAsItemOrFirstFromList,readFileIntoString
import argparse

def fasta_filter(fasta_file,name_list_path):
	name_list_str=readFileIntoString(name_list_path)
	name_list_arr=name_list_str.split('\n')
	#print "name list len =",len(name_list_arr)
	#print "name list content ",name_list_arr
	fasta_reader=SeqIO.parse(open(fasta_file, "r"), "fasta")
	#query_record=fasta_reader.next()
	for query_record in fasta_reader:
		#print "got a record : ",query_record.id
		if(query_record.id in name_list_arr):
			#print "in list"
			print ">"+str(query_record.id)+"\n"+query_record.seq
		else:
			#print "not in list"
			pass


if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Extract sequences from a fasta file by name ; write result to stdout')
	parser.add_argument('fasta',type=str,nargs=1,help="path to fasta")
	parser.add_argument('name_list',type=str,nargs=1,help="path to file of list of names (one per line)")
	args=parser.parse_args()
	if(args):
		fasta_file=extractAsItemOrFirstFromList(args.fasta)
		name_list_file=extractAsItemOrFirstFromList(args.name_list)
		fasta_filter(fasta_file,name_list_file)
	else:
		parser.print_help()




