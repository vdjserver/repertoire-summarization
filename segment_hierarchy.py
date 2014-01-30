#!/usr/bin/env python

from segment_utils import *
import os
import glob

#define inputs for run
vfasta="/home/esalina2/round1_imgt/human_IG_V.fna"
dfasta="/home/esalina2/round1_imgt/human_IG_D.fna"
jfasta="/home/esalina2/round1_imgt/human_IG_J.fna"
fasta_list=[vfasta,dfasta,jfasta]
blast_out="/home/esalina2/round1_imgt/all_data.processed.Q35.L200.R1.fna.igblastn.imgt.out"

#define inputs for database
base_dir="/home/data/DATABASE/01_22_2014/"
organism="human"


#do hierarchy structure I/O
pickle_file=base_dir+"/hierarchy_data.pkl"
if(os.path.exists(pickle_file)):
	#cache from pickle file
	#delete the file if patches are added so that the cache is 'rebuilt'
	print "Loading cached database from "+pickle_file+"..."
	hier_data=pickleRead(pickle_file)
	#hierarchy first in list, clone names in second, we want hierarchy in this case
	hierarchy=hier_data[0]
else:
	print "Database 'pickle' file "+pickle_file+" not found, so reading gene tables and patches...."
	#setup for gene tables parsing
	geneTables_dir=base_dir+"/"+organism+"/GeneTables"
	loci=['IGHD','IGHJ','IGHV','IGKJ','IGKV','IGLJ','IGLV','TRAJ','TRAV','TRBD','TRBJ','TRBV','TRDD','TRDJ','TRDV','TRGJ','TRGV']
	#loci=list()
	hierarchy=tree()
	for l in range(len(loci)):
		locus=loci[l]
		print "Now anlayzing for locus =",locus,"and organism =",organism
		locus_glob=geneTables_dir+"/"+locus+"*.html"
		html_files=glob.glob(locus_glob)
		for h in range(len(html_files)):
			html_file=html_files[h]
			print "Loading hierarchy from files :",html_file	
			URL="file://"+html_file
			hier_data=hierarchyTreeFromGenetableURL(URL,locus,None,hierarchy[organism][locus])
			hierarchy[organism][locus]=hier_data[0]
			patchPath=html_file+".patch"
			if(os.path.exists(patchPath)):
				print "Patch file "+patchPath+" found, so patching from it..."
				locusHierarchyData=patchlocusHierarchyData(hierarchy[organism][locus],patchPath)
				hierarchy[organism][locus]=locusHierarchyData
print "Using the hierarchy :"
prettyPrintTree(hierarchy)

#get segment counts
print "Now getting counts from ",blast_out
counts_map=get_segment_count_map_from_blast_output(blast_out,fasta_list)
for item in counts_map:
	print item+"->"+str(counts_map[item])+"\n"




JSON=jsonify_hierarchy(hierarchy['human'],"human",counts_map,"value")
print "\n\nTHIS IS THE JSON:\n\n\n"
print JSON
