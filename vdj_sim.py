#!/usr/bin/env python
import argparse
import random
from utils import read_fasta_file_into_map
import numpy

#global variable for V(D)J data!
vMap=None
dMap=None
jMap=None
bps=["A","C","G","T"]

#make a random junction string
def rand_jun(min_len,max_len):
	global bps	
	rand_len=random.choice(range(min_len,max_len+1))
	junc=""
	for i in range(rand_len):
		junc+=random.choice(bps)
	return junc


#pick V and J and recombine with a junction
#note, can currently allow IGHV to combine with IGKJ
def sim_light_recomb(vMap,jMap):
	vKey=random.choice(vMap.keys())
	jKey=random.choice(jMap.keys())
	vj_junc=rand_jun(0,10)
	light_seq=vMap[vKey]+vj_junc+jMap[jKey]
	return {
		'vKey':vKey,
		'dKey':"None",
		'jKey':jKey,
		'vj_junc':vj_junc,
		'seq':light_seq
		}
	#return [vKey,jKey,vj_junc,light_seq]
	
#pick V,D, and J and recombine with junctions
#note, can currently allow IGHV to combine with IGKJ
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
	



#dummy SHM
#TODO : region-aware SHM
def dummySHM(seq,mut_param_lambda):
	num_muts=numpy.random.poisson(mut_param_lambda,1)[0]
	if(num_muts==0):
		return seq
	mut_types=["ins","del","bsb"]
	shm_seq=list(seq)
	join_char=""
	global bps
	for mut_id in range(num_muts):
		bp_id=random.choice(range(len(shm_seq)))+1
		mut_type=random.choice(mut_types)
		if(mut_type=="bsb"):
			shm_seq[bp_id-1]=random.choice(bps)
		elif(mut_type=="ins"):
			new_nuc=random.choice(bps)
			shm_seq.insert(bp_id-1,new_nuc)
		else:
			#del #ts_shm=ts[0:1]+ts[2:]
			shm_seq=shm_seq[0:bp_id-1]+shm_seq[bp_id:]
	return join_char.join(shm_seq)




#from a full name extract the class
def getClassName(n):
	class_name=n[0:4]
	return class_name


#given a single map (not partitioned)
#partition into submaps
#for light/heavy loci respectively
#return a list of named maps
#Classificiation is based on name string comparison/extraction of the first 4 characters
def partitionIntoClassMaps(mainMap):
	#first gather the class names
	classList=list()
	for k in mainMap:
		class_name=getClassName(k)
		#print "From ",k," extracted ",class_name
		if(not(class_name in classList)):
			classList.append(class_name)	
	#initialize a map of maps based on class names
	mom=dict()
	for c in classList:
		mom[c]=dict()
	#now map into that map-of-maps the items from k
	for k in mainMap:
		this_class_name=getClassName(k)
		mom[this_class_name][k]=mainMap[k]
		#mom[c][k]="\n"
	return mom



#ae two classes recombination-compatible?
#if the first 3 characters are the same, then yes!
#otherwise no!
def compatibleForRecombination(c1,c2):
	a=c1[0:3]
	b=c2[0:3]
	if(a==b):
		return True
	else:
		return False




#using the input fastas and flags (and max sim)
#writer simulated VDJs to STDOUT
def vdj_sim(vFasta,dFasta,jFasta,no_light,no_heavy,max_sim=float("inf")):
	vMap=read_fasta_file_into_map(vFasta)
	vMom=partitionIntoClassMaps(vMap)
	print "VMOM : ",vMom
	sys.exit(0)
	jMap=read_fasta_file_into_map(jFasta)
	dMap=None
	sorted_h_keys=['vKey','dKey','jKey','vd_junc','dj_junc']
	sorted_l_keys=['vKey','jKey','vj_junc']
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
		descriptor=">"+str(num_sim+1)
		if(chain_selection=="heavy"):
			recomb=sim_heavy_recomb(vMap,dMap,jMap)
			for sk in sorted_h_keys:
				descriptor+="|"+sk+"="+recomb[sk]
			descriptor+="|chain_type=heavy"
		else:
			recomb=sim_light_recomb(vMap,jMap)
			for sk in sorted_l_keys:
				descriptor+="|"+sk+"="+recomb[sk]
			descriptor+="|chain_type=light"
		recomb['shm_seq']=dummySHM(recomb['seq'],0.5)
		print descriptor+"\n"+recomb['shm_seq']
		num_sim+=1
	


#get VDJ fasta, heavy/light flags and simulate
#with 'dumb' SHM (includes insertions, deletions, and base substitutions)
if (__name__=="__main__"):
	parser = argparse.ArgumentParser(description='Given V,D,J fastA files, simulate V(D)J recombination with somatic hypermutation.  Note SHM includes insertions, deletions, and base substitutions is not "region"-aware (e.g. of CDR1, CDR2, etc).  Simulated sequences are written to STDOUT')
	parser.add_argument('vfasta',type=str,nargs=1,help="path to the V fasta file")
	parser.add_argument('jfasta',type=str,nargs=1,help="path to the J fasta file")
	parser.add_argument('-dfasta',type=str,nargs=1,help="path to the D fasta file for heavy chains")
	parser.add_argument('-no_light',action='store_true',dest='no_light',help="exclude simulation of 'light' chains")
	parser.set_defaults(no_light=False)
	parser.add_argument('-no_heavy',action='store_true',dest='no_heavy',help="exclude simulation of 'heavy' chains")
	parser.set_defaults(no_heavy=False)	
	parser.add_argument('-num_seqs',type=int,default=float("inf"),nargs=1,help="the number of sequences to simulate")
	args = parser.parse_args()
	if(args.num_seqs!=float("inf")):
		max_sim=int(args.num_seqs[0])
	else:
		max_sim=float("inf")
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
		vdj_sim(vFasta,dFasta,jFasta,args.no_light,args.no_heavy,max_sim)


	


