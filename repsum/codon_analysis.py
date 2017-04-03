#!/usr/bin/env python

from utils import *
from alignment import *
from char_utils import *
from segment_utils import getVRegionStartAndStopGivenRefData,getTheFrameForThisReferenceAtThisPosition,getTheFrameForThisJReferenceAtThisPosition
from vdjml.igblast_parse import rev_comp_dna
import math

#make an instance of the analyzer
#print "INIT IN CODON_ANALYSIS...."
AGSCodonTranslator=CodonAnalysis()
#print "FINISH INIT IN CODON_ANALYSIS...."

#find correct number map based on gene call




#given a region name and length, return the numbering
def acquireNumberingMap(reg_name,align_len,vInfo,organism,imgtdb_obj):
	pos=list()
	region_numbers=list()
	info=reg_name.split("_")
	
	if(vInfo['subject ids'].startswith("IGHV")):
		if(reg_name=="FR1_kabat"):
			region_numbers.extend(range(1,31))
		elif(reg_name=="CDR1_kabat"):
			region_numbers.extend(range(31,36))
			pos=getVRegionStartAndStopGivenRefData(vInfo['subject ids'],organism,imgtdb_obj,info[0],info[1])
			full_length=int((pos[1]-pos[0]+1)/3)
			if(full_length==7):
				region_numbers.extend(["35A","35B"])
			elif(full_length==6):
				region_numbers.extend("35A")
			else:
				pass
		elif(reg_name=="FR2_kabat"):
			region_numbers.extend(range(36,50))
		elif(reg_name=="CDR2_kabat"):
			region_numbers.extend([50,51,52,"52A","52B","52C"])
			pos=getVRegionStartAndStopGivenRefData(vInfo['subject ids'],organism,imgtdb_obj,info[0],info[1])
			insert_size=int((pos[1]-pos[0]+1)/3)-15
			if(insert_size>0):
				for x in xrange (1,insert_size):
					region_numbers.pop()
			region_numbers.extend(range(53,66))
		elif(reg_name=="FR3_kabat"):
			region_numbers.extend(range(66,82))
			region_numbers.extend(["82A","82B","82C"])
			region_numbers.extend(range(83,95))
		else:
			#invalid/unknown region!
			#print "Error, unknown region ",reg_name,"!"
			pass
	elif(vInfo['subject ids'].startswith("IGKV")):
		pass
	elif(vInfo['subject ids'].startswith("IGLV")):
		pass
	else:
		pass
		#invalid/unknown chain!
		#print "Error, unknown gene ",vInfo['subject ids'],"!"
	
#	VH
#	30 (30)
#	31-35,i35A,i35B (5-7)
#	36-49 (14)
#	50-51,i52,i52a,i52b,i52c,53-65 (15-19)
#	66-82,82A,82B,82C,83-94 (32)
#	check CDR1 and CDR2
#	
#	VK
#	23 (23)
#	24-27,27A-F,28-34 (11-17)
#	35-49 (15)
#	50-56 (7)
#	57-88 (32)
#	check CDR1
#	
#	VL
#	1-9,11-23 (22)
#	24-27,27A-F,28-34 (11-17)
#	35-49 (15)
#	50-56 (7)
#	57-88 (32)
#	check CDR1
#	
#	CDR2: L4,L5,L11 = 11
#	CDR2: L9 = 12
#	FR3: L5,L6,L11 = 34
#	
#	numbering=list()
#	letters=["A","B","C"]
#	l_pos=0
#	if(reg_name=="FR1_kabat"):
#		print "FR1_kabat"
#	elif(reg_name=="CDR1_kabat"):
#		print "CDR1_kabat"
#		if(reg_len==5):
#			for p in range(31,36):
#				numbering.append(str(p))
#		elif(reg_len==6):
#			new_nums=["31","31A","32","33","34","35"]
#			for nn in new_nums:
#				numbering.append(nn)
#		elif(reg_len==7):
#			new_nums=["31","31A","31B","32","33","34","35"]
#			for nn in new_nums:
#				numbering.append(nn)					
#	elif(reg_name=="FR2_kabat"):
#		print "FR2_kabat"
#		for p in range(36,50):
#			numbering.append(str(p))
#	elif(reg_name=="CDR2_kabat"):
#		print "CDR2_kabat"
#		for p in range(50,67):
#			numbering.append(str(p))
#	elif(reg_name=="FR3_kabat"):
#		print "FR3_kabat"
#		for p in range(66,83):
#			numbering.append(str(p))
#		for l in range(len(letters)):
#			numbering.append("82"+str(letters[l]))
#		for p in range(83,93):
#			numbering.append(str(p))
#	else:
#		#invalid/unknown region!
#		print "Error, unknown region ",reg_name,"!"
#	self.numberingMapCache[reg_name+str(reg_len)]=numbering
#	return self.numberingMapCache[reg_name+str(reg_len)]
	return region_numbers



#from the regions, acquire a mutation map
def acquire_mutation_map(reg_names,reg_infos,reg_numbers):
	AA_map=list()
	codon_map=list()
	AA_silent_map=list()
	codon_silent_map=list()
	
	for r in range(len(reg_names)):
		reg_info_map=reg_infos[r].getCharMap()
#		print reg_info_map
#		numbering_list=acquireNumberingMap(reg_names[r],reg_info_map['AA'])
		q_codons=reg_info_map['nucleotide read2']
		s_codons=reg_info_map['subject_read2']
		q_aminos=reg_info_map['AA2']
		s_aminos=reg_info_map['AA_ref2']
		numbered_pos=reg_numbers[r]
		#print q_codons
		#print s_codons
		#print q_aminos
		#print s_aminos
		#print numbered_pos
		
		#Check each codon for any base substitutions
		for ci in range(len(q_codons)):
			if(s_codons[ci]!=q_codons[ci] and q_aminos[ci]!="X"):
				#Save mutation information as "codon/amino FROM" "codon number" "codon/amino TO"
				aaP=s_aminos[ci]+str(numbered_pos[ci])+q_aminos[ci]
				cdP=s_codons[ci]+str(numbered_pos[ci])+q_codons[ci]
				#One or more base substitutions cause a replacement mutation
				if(s_aminos[ci]!=q_aminos[ci]):
					AA_map.append(aaP)
					codon_map.append(cdP)
				#One or more base substitutions result in a silent mutation
				else:
					AA_silent_map.append(aaP)
					codon_silent_map.append(cdP)
			else:
				pass
	
	overall_map=dict()
	overall_map['codons']="None"
	overall_map['aminos']="None"
	overall_map['codons_silent']="None"
	overall_map['aminos_silent']="None"
	
	if(len(codon_map)>0):
		overall_map['codons']=codon_map
	if(len(AA_map)>0):
		overall_map['aminos']=AA_map
	if(len(codon_silent_map)>0):
		overall_map['codons_silent']=codon_silent_map
	if(len(AA_silent_map)>0):
		overall_map['aminos_silent']=AA_silent_map
	
	return overall_map



def annotationMutationMap(vInfo,dInfo,jInfo,alignment_output_queue,num_submitted_jobs,imgtdb_obj,organism,read_rec,cdr3_map,homology_filter_val):
	if('subject ids' in vInfo):
		filterNote=""
		get_res=0
		reg_names=[]
		reg_infos=[]
		reg_numbers=[]
		numbering_map=dict()
		numbering_map['FR1_imgt']=list()
		numbering_map['CDR1_imgt']=list()
		numbering_map['FR2_imgt']=list()
		numbering_map['CDR2_imgt']=list()
		numbering_map['FR3_imgt']=list()
		
		#Find and load IMGT germline database numbers for a specific germline
		lookupPath=imgtdb_obj.getBaseDir()+"/"+organism+"/ReferenceDirectorySet/"+"IGHV.fna"
		if(not(os.path.exists(lookupPath)) or not(os.path.isfile(lookupPath))):
			raise Exception("Error, failed to find germline file "+lookupPath+" for organism = "+organism)
		else:
			germline_reader=open(lookupPath,'r')
			germline_found=False
			for line in germline_reader:
				line=line.strip()
				if(germline_found==True):
					#Start number range at first full codon using germline alignment start position information which ignores gaps
					#so gaps are ignored for the alignment, using "gapless_numbers" but counted to ensure proper numbering.
					alignment_start=int(math.ceil((vInfo['s. start']-1)/3.0))+1
					gapless_numbers=0
					#Iterate number map over whole germline up to last full codon
					for r in range(1,int(math.floor(len(line)/3.0))+1):
						#Parse germline by codon starting from positions (0;1;2)
						index_start=(r-1)*3
						index_end=r*3
						#Remove codons with any gaps ("."), since partial codon alignments are discarded for replacement and silent mutation reporting.
						if(line[index_start:index_end].find(".")==(-1)):
							gapless_numbers+=1
							if(gapless_numbers<alignment_start):
								pass
							else:
								if(r<27):
									numbering_map['FR1_imgt'].append(r)
								elif(r<39):
									numbering_map['CDR1_imgt'].append(r)
								elif(r<56):
									numbering_map['FR2_imgt'].append(r)
								elif(r<66):
									numbering_map['CDR2_imgt'].append(r)
								else:
									numbering_map['FR3_imgt'].append(r)
					break
				#Look for germline in full database until found
				elif(line[0]==">"):
					line_pieces=line.split('|')
					if(line_pieces[1]==vInfo['subject ids']):
						germline_found=True
			germline_reader.close()
		
		#Identify all regions which have alignment information for the query
		while(get_res<num_submitted_jobs):
			region_alignment=alignment_output_queue.get()
			if(region_alignment is not None):
				if(region_alignment.getName().startswith("FR1_imgt")):
					reg_names.append("FR1_imgt")
					reg_infos.append(region_alignment)
					reg_numbers.append(numbering_map['FR1_imgt'])
				elif(region_alignment.getName().startswith("CDR1_imgt")):
					reg_names.append("CDR1_imgt")
					reg_infos.append(region_alignment)
					reg_numbers.append(numbering_map['CDR1_imgt'])
				elif(region_alignment.getName().startswith("FR2_imgt")):
					reg_names.append("FR2_imgt")
					reg_infos.append(region_alignment)
					reg_numbers.append(numbering_map['FR2_imgt'])
				elif(region_alignment.getName().startswith("CDR2_imgt")):
					reg_names.append("CDR2_imgt")
					reg_infos.append(region_alignment)
					reg_numbers.append(numbering_map['CDR2_imgt'])
				elif(region_alignment.getName().startswith("FR3_imgt")):
					reg_names.append("FR3_imgt")
					reg_infos.append(region_alignment)
					reg_numbers.append(numbering_map['FR3_imgt'])
				else:
					pass
			get_res+=1
		
		
		mutation_map=acquire_mutation_map(reg_names,reg_infos,reg_numbers)
		return [filterNote,mutation_map]
	else:
		filterNote="No Vgene"
		return [filterNote,mutation_map]





