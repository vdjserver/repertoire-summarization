#!/usr/bin/env python

from utils import *
from alignment import *
from char_utils import *
from segment_utils import getVRegionStartAndStopGivenRefData,getTheFrameForThisReferenceAtThisPosition,getTheFrameForThisJReferenceAtThisPosition
from vdjml.igblast_parse import rev_comp_dna
from imgt_utils import imgt_db
import math



#From the regions, acquire a mutation map
def acquire_mutation_map(reg_names,reg_infos,reg_numbers):
	AA_map=list()
	codon_map=list()
	AA_silent_map=list()
	codon_silent_map=list()
	
	for r in range(len(reg_names)):
		reg_info_map=reg_infos[r].getCharMap()
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
	mutation_map=dict()
	
	#Check if Vgene call exists first
	if('subject ids' in vInfo):
		germline_numbers=imgtdb_obj.getGermlineNumberMap(organism,vInfo['subject ids'])
		#Make a copy to trim without altering germline reference set
		gene_numbers=dict(germline_numbers)
		
		if(gene_numbers==None):
			mutation_map['codons']="No gene numbers found"
			mutation_map['aminos']="No gene numbers found"
			mutation_map['codons_silent']="No gene numbers found"
			mutation_map['aminos_silent']="No gene numbers found"
			return mutation_map
		else:
			get_res=0
			reg_names=[]
			reg_infos=[]
			reg_numbers=[]
			
			#Find number of positions to trim from first region in alignment
			remove_numbers=int(math.ceil(vInfo['s. start']-1)/3.0)
			germline_length=[len(gene_numbers['FR1_imgt']),len(gene_numbers['CDR1_imgt']),len(gene_numbers['FR2_imgt']),len(gene_numbers['CDR2_imgt']),len(gene_numbers['FR3_imgt'])]
			r=0
			while(r<5):
				if(remove_numbers<germline_length[r]):
					break
				else:
					remove_numbers=remove_numbers-germline_length[r]
					r+=1
			if(remove_numbers>0):
				trim_test=True
			else:
				trim_test=False
			
			#Identify all regions which have alignment information for the query
			#Also trim number list for first region with partial alignment
			while(get_res<num_submitted_jobs):
				region_alignment=alignment_output_queue.get()
				if(region_alignment is not None):
					if(region_alignment.getName().startswith("FR1_imgt")):
						reg_names.append("FR1_imgt")
						reg_infos.append(region_alignment)
						if(trim_test is True):
							reg_numbers.append(gene_numbers['FR1_imgt'][remove_numbers:])
							#print "Trimmed",gene_numbers['FR1_imgt'][remove_numbers:]
							trim_test=False
						else:
							reg_numbers.append(gene_numbers['FR1_imgt'])
							#print "Untrimmed",gene_numbers['FR1_imgt']
					elif(region_alignment.getName().startswith("CDR1_imgt")):
						reg_names.append("CDR1_imgt")
						reg_infos.append(region_alignment)
						if(trim_test is True):
							reg_numbers.append(gene_numbers['CDR1_imgt'][remove_numbers:])
							#print "Trimmed",gene_numbers['CDR1_imgt'][remove_numbers:]
							trim_test=False
						else:
							reg_numbers.append(gene_numbers['CDR1_imgt'])
							#print "Untrimmed",gene_numbers['CDR1_imgt']
					elif(region_alignment.getName().startswith("FR2_imgt")):
						reg_names.append("FR2_imgt")
						reg_infos.append(region_alignment)
						if(trim_test is True):
							reg_numbers.append(gene_numbers['FR2_imgt'][remove_numbers:])
							#print "Trimmed",gene_numbers['FR2_imgt'][remove_numbers:]
							trim_test=False
						else:
							reg_numbers.append(gene_numbers['FR2_imgt'])
							#print "Untrimmed",gene_numbers['FR2_imgt']
					elif(region_alignment.getName().startswith("CDR2_imgt")):
						reg_names.append("CDR2_imgt")
						reg_infos.append(region_alignment)
						if(trim_test is True):
							reg_numbers.append(gene_numbers['CDR2_imgt'][remove_numbers:])
							#print "Trimmed",gene_numbers['CDR2_imgt'][remove_numbers:]
							trim_test=False
						else:
							reg_numbers.append(gene_numbers['CDR2_imgt'])
							#print "Untrimmed",gene_numbers['CDR2_imgt']
					elif(region_alignment.getName().startswith("FR3_imgt")):
						reg_names.append("FR3_imgt")
						reg_infos.append(region_alignment)
						if(trim_test is True):
							reg_numbers.append(gene_numbers['FR3_imgt'][remove_numbers:])
							#print "Trimmed",gene_numbers['FR3_imgt'][remove_numbers:]
							trim_test=False
						else:
							reg_numbers.append(gene_numbers['FR3_imgt'])
							#print "Untrimmed",gene_numbers['FR3_imgt']
					else:
						pass
				get_res+=1
			
			mutation_map=acquire_mutation_map(reg_names,reg_infos,reg_numbers)
			return mutation_map
	else:
		mutation_map['codons']="No Vgene found"
		mutation_map['aminos']="No Vgene found"
		mutation_map['codons_silent']="No Vgene found"
		mutation_map['aminos_silent']="No Vgene found"
		return mutation_map



