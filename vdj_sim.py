#!/usr/bin/env python
import argparse

#global variable for V(D)J data!
vMap=None
dMap=None
jMap=None

#read_fasta_file_into_map(fasta_path,alwaysSeqToUpper=True,removeNonIUPAC=True):
sim_light(vMap,jMap):
	vIndex=sample()




if (__name__=="__main__"):
	parser = argparse.ArgumentParser(description='Given V,D,J fastA files, simulate V(D)J recombination with somatic hypermutation')
	parser.add_argument('vfasta',type=str,nargs=1,help="path to the V fasta file")
	parser.add_argument('jfasta',type=str,nargs=1,help="path to the J fasta file")
	parser.add_argument('-dfasta',type=str,nargs=1,help="path to the D fasta file for heavy chains")
	parser.add_argument('-no_light',action='store_true',dest='no_light',help="exclude simulation of 'light' chains")
	parser.set_defaults(no_light=False)
	parser.add_argument('-no_heavy',action='store_true',dest='no_heavy',help="exclude simulation of 'heavy' chains")
	parser.set_defaults(no_heavy=False)	
	args = parser.parse_args()
	if(args.no_light and args.no_heavy):
		print "ERROR! Must allow simulation of at least either light or heavy chains!"
		parser.print_help()
	


