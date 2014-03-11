#!/usr/bin/env python
import argparse
import random
from utils import read_fasta_file_into_map

#global variable for V(D)J data!
vMap=None
dMap=None
jMap=None

def rand_jun(min_len,max_len):
	bps=["A","C","G","T"]
	rand_len=random.choice(range(min_len,max_len+1))
	junc=""
	for i in range(rand_len):
		junc+=random.choice(bps)
	return junc



def sim_light_recomb(vMap,jMap):
	vKey=random.choice(vMap.keys())
	jKey=random.choice(jMap.keys())
	vj_junc=rand_jun(0,10)
	light_seq=vMap[vKey]+vj_junc+jMap[jKey]
	return [vKey,jKey,vj_junc,light_seq]
	

def sim_heavy_recomb(vMap,dMap,jMap):
	vKey=random.choice(vMap.keys())
	dKey=random.choice(dMap.keys())
	jKey=random.choice(jMap.keys())
	vd_junc=rand_jun(0,10)
	dj_junc=rand_jun(0,10)
	heavy_seq=vMap[vKey]+vd_junc+dMap[dKey]+dj_junc+jMap[jKey]
	return {
		'vKey':vKey,
		'dKey':dKey,
		'jKey':jKey,
		'vd_junc':vd_junc,
		'dj_junc':dj_junc,
		'seq':heavy_seq
		}
	#return [vKey,dKey,jKey,vd_junc,dj_junc,heavy_seq]
	

def vdj_sim(vFasta,dFasta,jFasta,no_light,no_heavy,max_sim=float("inf")):
	vMap=read_fasta_file_into_map(vFasta)
	jMap=read_fasta_file_into_map(jFasta)
	dMap=None
	if(no_heavy):
		chain_list=["light"]
	elif(no_light):
		chain_list=["heavy"]
	else:
		chain_list=["ligth","heavy"]
	if("heavy" in chain_list):
		dMap=read_fasta_file_into_map(dFasta)
	num_sim=0
	while(num_sim<max_sim):
		chain_selection=random.choice(chain_list)
		print "The chain selection is ",chain_selection
		if(chain_selection=="heavy"):
			heavy_recomb=sim_heavy_recomb(vMap,dMap,jMap)
			print heavy_recomb
		num_sim+=1
	



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
	elif(not(args.no_heavy) and not(args.dfasta)):
		print "ERROR! Must supply D fasta if not excluding heavy chain simulation!"
		parser.print_help
	else:
		vFasta=args.vfasta[0]
		dFasta=None
		if(args.dfasta):
			dFasta=args.dfasta[0]
		jFasta=args.jfasta[0]
		print "vdj are ",vFasta," and ",dFasta," and ",jFasta
		vdj_sim(vFasta,dFasta,jFasta,args.no_light,args.no_heavy,max_sim=float("inf"))

	


