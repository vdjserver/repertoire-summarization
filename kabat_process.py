#!/usr/bin/env python

import re


def isIntegral(s):
	ire=re.compile(r'^\s*(\d+)\s*$')
	sr=re.search(ire,s)
	if(sr):
		return True
	else:
		return False


def writeKabatRegionsFromIGBLASTKabatResult(k,o):
	reg_map=dict()
	regions=["FWR1","CDR1","FWR2","CDR2","FWR3","CDR3"]
	currentQuery=None
	currentMode=None
	reader=open(k,'r')
	for line in reader:
		line=line.strip()
		if(line.startswith("#")):
			if(line.startswith("# Alignment summary")):
				currentMode="summary"
			elif(line.startswith("# Query:")):
				## Query: IGHV5-a*04
				qre=re.compile(r'^#\s+Query:\s*(.*)$')
				rr=re.search(qre,line)
				if(rr):
					currentQuery=rr.group(1)
					currentQuery=currentQuery.strip()
			else:
				currentMode=None
		if(currentMode=="summary"):
			for region in regions:
				if(line.startswith(region)):
					if(not(currentQuery in reg_map)):
						reg_map[currentQuery]=dict()
					if(not(region in reg_map[currentQuery])):
						reg_map[currentQuery][region]=dict()
					pieces=line.split("\t")
					start=pieces[1]
					stop=pieces[2]
					if(isIntegral(start)):
						reg_map[currentQuery][region]["start"]=int(start)
					else:
						reg_map[currentQuery][region]["start"]=(-1)
					print "got start ",start," for q=",currentQuery," region=",region
					if(not(region=="CDR3")):
						if(isIntegral(stop)):
							reg_map[currentQuery][region]["stop"]=int(stop)
						else:
							reg_map[currentQuery][region]["stop"]=(-1)
						print "got stop ",stop," for q=",currentQuery," region=",region
					else:
						reg_map[currentQuery][region]["stop"]=(-1)
	key_list=reg_map.keys()
	Key_list.sort()
	writer=open(o,'w')
	for i in range(len(key_list)):
		k=key_list[i]
		writer.write(k+"\t")
		for region in regions
			writer.write(reg_map[k][region]["start"]+"\t"+reg_map[k][region]["stop"]
		if(i<len(key_list)-1):
			writer.write("\n")
	writer.close()

			

					
					
	
	



igblast="/home/data/DATABASE/01_22_2014/human/ReferenceDirectorySet/KABAT/igblastn.kabat.out"
writeKabatRegionsFromIGBLASTKabatResult(igblast,"/dev/null")
