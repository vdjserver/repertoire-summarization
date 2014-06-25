#!/usr/bin/env python
import argparse
import random
from utils import read_fasta_file_into_map
from imgt_utils import get_loci_list,get_heavy_loci
import numpy
import sys


#global variable for V(D)J data!
vMap=None
dMap=None
jMap=None
bps=["A","C","G","T"]
seg_counts=dict()


#given a map of (integral) counts return a frequence map (for normalization)
def get_freqs_from_counts(counts):
	total=0
	for s in counts:
		total+=counts[s]
	freqs=dict()
	for s in counts:
		freqs[s]=float(counts[s])/float(total)
	return freqs


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
def getLocusClassName(n):
	class_name=n[0:3]
	return class_name


#given a single map (not partitioned)
#partition into submaps
#for light/heavy loci respectively
#return a list of named maps
#Classificiation is based on name string comparison/extraction of the first 4 characters
def partitionIntoClassMaps(mainMap):
	#first gather the class names
	classList=list()
	allowed_classes=getDefaultClasses()
	for k in mainMap:
		class_name=getLocusClassName(k)
		#print "From ",k," extracted ",class_name
		if(not(class_name in classList)):
			classList.append(class_name)
		if(not(class_name in allowed_classes)):
			raise Exception("Error, found sequence named ",k," with unknown class!")
	#initialize a map of maps based on class names
	mom=dict()
	for c in classList:
		mom[c]=dict()
	#now map into that map-of-maps the items from k
	for k in mainMap:
		this_class_name=getLocusClassName(k)
		mom[this_class_name][k]=mainMap[k]
		#mom[c][k]="\n"
	return mom



#are two classes recombination-compatible?
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
def vdj_sim(vFasta,dFasta,jFasta,selected_loci_classes,mut_lam,max_sim=float("inf")):
	global seg_counts
	vMap=read_fasta_file_into_map(vFasta)
	vMom=partitionIntoClassMaps(vMap)
	num_miss_V=0
	for v_class in vMom.keys():
		if(not(v_class in selected_loci_classes)):
			num_miss_V+=1
	if(num_miss_V==len(vMom.keys())):
		print "ERROR, loci classes from selected loci classes is ",selected_loci_classes," but classes read from V FASTA are ",vMom.keys(),". No data found for selected classes ! Abort!"
		import sys
		sys.exit(0)
	jMap=read_fasta_file_into_map(jFasta)
	jMom=partitionIntoClassMaps(jMap)
	num_miss_J=0
	for j_class in jMom.keys():
		if(not(j_class in selected_loci_classes)):
			num_miss_J+=1
	if(num_miss_J==len(jMom.keys())):
		print "ERROR, loci classes from selected loci classes is ",selected_loci_classes," but classes read from J FASTA are ",jMom.keys(),". No data found for selected classes ! Abort!"
		import sys
		sys.exit(0)
	dMap=None
	sorted_h_keys=['vKey','dKey','jKey','vd_junc','dj_junc']
	sorted_l_keys=['vKey','jKey','vj_junc']
	if(not(dFasta is None) and isAtLeastOneSelectedLocusClassHeavy(selected_loci_classes)):
		dMap=read_fasta_file_into_map(dFasta)
		dMom=partitionIntoClassMaps(dMap)
		num_miss_D=0
		for d_class in dMom.keys():
			if(not(d_class in selected_loci_classes)):
				num_miss_D+=1
		if(num_miss_D==len(dMom.keys())):
			print "ERROR, loci classes from selected loci classes is ",selected_loci_classes," but classes read from D FASTA are ",dMom.keys(),". No data found for selected classes ! Abort!"
			import sys
			sys.exit(0)
	num_sim=0
	fasta_available_loci_classes=list()
	for m in vMom:
		if(not(m in fasta_available_loci_classes)):
			fasta_available_loci_classes.append(m)
	for m in jMom:
		if(not(m in fasta_available_loci_classes)):
			fasta_available_loci_classes.append(m)
	selected_and_available_loci_classes=list()
	for a in fasta_available_loci_classes:
		if(not(a in selected_loci_classes)):
			#if something is available, but not selected, then give no warning
			pass
		else:
			selected_and_available_loci_classes.append(a)
	list_of_selected_but_unavailble=list()
	for s in selected_loci_classes:
		if(not(s in fasta_available_loci_classes)):
			list_of_selected_but_unavailble.append(s)
	if(len(list_of_selected_but_unavailble)>0):
		import sys
		warn_msg="WARNING : The following loci are selected, but unavailable from the supplied FASTA files : "+str(list_of_selected_but_unavailble)
		sys.stderr.write(warn_msg+"\n")
	if(len(selected_and_available_loci_classes)==0):
		warn_msg="ERROR : availble loci (from FASTA) are : "+str(fasta_available_loci_classes)+" but selected loci are "+str(selected_loci_classes)+"!"
		import sys
		sys.stderr.write(warn_msg+"\n")
		sys.exit(2)
	while(num_sim<max_sim):
		#print (num_sim+1)
		#make a random class choice from the list of selected classes
		#MOM='map of maps'
		selected_random_class=random.choice(selected_and_available_loci_classes)
		vPool=vMom[selected_random_class]
		jPool=jMom[selected_random_class]
		dPool=None
		if(locusClassIsHeavy(selected_random_class)):
			dPool=dMom[selected_random_class]
		#print "SELECTED RANDOM CLASS IS ",selected_random_class
		#print "V POOL ",vPool
		#print "D POOL ",dPool
		#print "J POOL ",jPool
		descriptor=">"+str(num_sim+1)
		if(dPool is not None):
			#heavy recombination
			recomb=sim_heavy_recomb(vPool,dPool,jPool)
			for sk in sorted_h_keys:
				descriptor+="|"+sk+"="+recomb[sk]
			descriptor+="|chain_type=heavy"
		else:
			#light recombination
			recomb=sim_light_recomb(vPool,jPool)
			for sk in sorted_l_keys:
				descriptor+="|"+sk+"="+recomb[sk]
			descriptor+="|chain_type=light"
		recomb['shm_seq']=dummySHM(recomb['seq'],0.5)
		print descriptor+"\n"+recomb['shm_seq']
		seg_keys=['vKey','dKey','jKey']
		for s in seg_keys:
			if(s in recomb):
				segment=recomb[s]
				if(segment!="None"):
					if(segment in seg_counts):
						seg_counts[segment]+=1
					else:
						seg_counts[segment]=1						
		num_sim+=1
	

#given a locus, get the class
def getLocusClass(s):
	return s[0:3]



def isAtLeastOneSelectedLocusClassHeavy(slc):
	for s in slc:
		if(locusClassIsHeavy(s)):
			return True
	return False

#get default classes from IMGT loci
def getDefaultClasses():
	loci=get_loci_list()
	classes=list()
	for locus in loci:
		class_name=getLocusClass(locus)
		if(not(class_name in classes)):
			classes.append(class_name)
	return classes


#determine if a locus class is heavy
def locusClassIsHeavy(lc):
	heavy_loci=get_heavy_loci()
	for heavy_locus in heavy_loci:
		heavy_locus_class=getLocusClass(heavy_locus)
		if(heavy_locus_class==lc):
			return True
	return False






#extract from a list the items identified as "heavy"
def returnHeavyClassItems(l):
	hi=list()
	heavy_loci=get_heavy_loci()
	heavy_classes=list()
	for hl in heavy_loci:
		heavy_classes.append(getLocusClass(hl))
	for i in l:
		if(i in heavy_classes):
			hi.append(i)
	return hi



#get VDJ fasta, heavy/light flags and simulate
#with 'dumb' SHM (includes insertions, deletions, and base substitutions)
if (__name__=="__main__"):
	parser = argparse.ArgumentParser(description='Given V,D,J fastA files, simulate V(D)J recombination with somatic hypermutation.  Note SHM includes insertions, deletions, and base substitutions is not "region"-aware (e.g. of CDR1, CDR2, etc).  Since SHM is not "region"-aware there is *no* bias for the positions of any of the mutations (base substitions, insertions, deletions). Simulated sequences are written to STDOUT')
	parser.add_argument('vfasta',type=str,nargs=1,help="path to the V fasta file")
	parser.add_argument('jfasta',type=str,nargs=1,help="path to the J fasta file")
	parser.add_argument('-dfasta',type=str,nargs=1,help="path to the D fasta file for heavy chains")
	parser.add_argument('-loci',type=str,nargs=1,help="comma-separated list of allowed loci classes (default these 7 : \"IGH,IGK,IGV,TRA,TRB,TRD,TRG\")")
	parser.add_argument('-num_seqs',type=int,default=float("inf"),nargs=1,help="the number of sequences to simulate (default : infinity!) ")
	parser.add_argument('-mut_lam',type=float,default=float(0.5),nargs=1,help="the lambda parameter for SHM")
	args = parser.parse_args()
	mut_lam=None
	if(type(args.mut_lam)==list):
		mut_lam=args.mut_lam[0]
	else:
		#print "args.mut_lam is ",args.mut_lam
		mut_lam=args.mut_lam
	if(mut_lam<0):
		print "ERROR, selected mut_lam is negative! Abort !"
		import sys
		sys.exit(1)
	if(not(args.loci is None)):
		selected_loci_classes=args.loci
		selected_loci_classes=selected_loci_classes[0]
		selected_loci_classes=selected_loci_classes.split(',')
		allowed_loci_classes=getDefaultClasses()
		for s in selected_loci_classes:
			if(not(s in allowed_loci_classes)):
				print "ERROR, selected loci ",s," not in list of allowed loci : ",allowed_loci_classes
				import sys
				sys.exit(1)
	else:
		selected_loci_classes=getDefaultClasses()
	if(args.num_seqs!=float("inf")):
		max_sim=int(args.num_seqs[0])
	else:
		max_sim=float("inf")
	heavyClassItems=returnHeavyClassItems(selected_loci_classes)
	if(len(heavyClassItems)>0):
		#ensure the dfasta is set
		if(args.dfasta is None):
			print "The following selected loci are heavy ("+str(heavyClassItems)+") but no DFASTA has been specified ! Heavy loci require D fasta! Abort !"
			import sys
			sys.exit(1)
	vFasta=args.vfasta[0]
	dFasta=None
	if(args.dfasta):
		dFasta=args.dfasta[0]
	jFasta=args.jfasta[0]
	vdj_sim(vFasta,dFasta,jFasta,selected_loci_classes,mut_lam,max_sim)
	freqs=get_freqs_from_counts(seg_counts)
	seg_keys=seg_counts.keys()
	seg_keys.sort()
	sys.stderr.write("Segment\tFrequency\n")
	for seg in seg_keys:
		sys.stderr.write(seg+"\t"+str(freqs[seg])+"\n")
	

	


