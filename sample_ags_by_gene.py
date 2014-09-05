#!/usr/bin/env python


import argparse
from utils import extractAsItemOrFirstFromList
import os
from ags_mgr import ags_manager
import re





if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Make gene-level AGS report using rep_char output')
	parser.add_argument('tsv_path',type=str,nargs=1,help="path to rep_char output")
	parser.add_argument('gene_str',type=str,nargs=1,help="regex for gene test/filter")
	#vnames=["IGHV4\-4","IGHV4\-30","IGHV4\-31","IGHV4\-34","IGHV4\-38","IGHV4\-39","IGHV4\-59","IGHV4\-61"]
	#jnames=["IGHJ1","IGHJ2","IGHJ3","IGHJ4","IGHJ5","IGHJ6"]
	myAgsMgr=ags_manager("TNAME")
	args=parser.parse_args()
	if(args):
		tsv_path=extractAsItemOrFirstFromList(args.tsv_path)
		if(not(os.path.exists(tsv_path))):
			raise Exception("Error, file ",tsv_path," not found!")
		else:
			#print "found ",tsv_path
			gene_filter_name=extractAsItemOrFirstFromList(args.gene_str)
			gene_filter_name_escaped=re.sub(r'\-','\\\-',gene_filter_name)
			#print "the name ",gene_filter_name
			#print "the escaped",gene_filter_name_escaped
			regex_to_use="^"+gene_filter_name_escaped+"[^0-9]"
			reader=open(tsv_path,'r')
			line_num=1
			vIndex=None
			jIndex=None
			mutIndex=None
			okIndex=None
			for line in reader:
				temp_line=line
				temp_line=temp_line.strip()
				pieces=temp_line.split('\t')
				if(line_num==1):
					for p in range(len(pieces)):
						#print str(int(p+1)),"\t'"+pieces[p]+"'"
						if(pieces[p]=="V gene (highest score)"):
							vIndex=p
						if(pieces[p]=="J gene (highest score)"):
							jIndex=p
						if(pieces[p]=="AGS Amino Muts"):
							mutIndex=p
						if(pieces[p]=="AGS Q Note"):
							okIndex=p
				else:
					if(mutIndex==None or jIndex==None or vIndex==None):
						print "ln=",line_num
						print "mut=",mutIndex
						print "j=",jIndex
						print "v=",vIndex
						raise Exception("Fail on obtain index ind in ",tsv_path)
					else:
						status=pieces[okIndex]
						if(status=="OK"):
							#print "line_num=",line_num
							vGene=pieces[vIndex]
							jGene=pieces[jIndex]
							#print vGene,"\t",jGene
							muts=eval(pieces[mutIndex])
							matches_on_j=False
							matches_on_v=False
							if(re.search(regex_to_use,vGene)):
								print vGene," matches with ",regex_to_use
								matches_on_v=True
							else:
								#print "NOT", vGene," matches with ",regex_to_use
								pass
							if(re.search(regex_to_use,jGene)):
								matches_on_j=True
							if(matches_on_v or matches_on_j):
								for mut in muts:
									myAgsMgr.receive_numbered_mut(mut)
				line_num+=1
		print tsv_path+"\t"+gene_filter_name+"\t"+str(myAgsMgr.compute_ags()[0])
	#VH4-4
	#VH4-30 (combine all further subdivisions of 4-30 together)
	#VH4-31
	#VH4-34
	#VH4-38
	#VH4-39
	#VH4-59
	#VH4-61
	#
	#For JHgenes, use:
	#JH1
	#JH2
	#JH3
	#JH4
	#JH5
	#JH6
	



