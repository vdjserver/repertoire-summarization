#!/usr/bin/env python

from utils import *
from alignment import *


#make an instance of the analyzer
global codonAnalyzer

#class for codon counting
#currently for IGHV4 only
class codonCounter:

	gap_kabat=list()
	kabat=list()
	region_kabat=list()
	chothia=list()
	gap_chothia=list()
	region_chothia=list()
	allowGaps=False

	#init
	def __init__(self,pos_file_path,init_allowGaps=False):
		#currently reads a TSV
		#6 fields
		#gap order	KABAT	REGION_KABAT	CHOTHIA	gap order	REGION_CHOTHIA
		reader=open(pos_file_path,'r')
		line_num=1
		for line in reader:
			#print line
			temp=line
			pieces=temp.split('\t')
			for i in range(len(pieces)):
				pieces[i]=pieces[i].strip()
			if(line_num!=1):
				#ignore header line
				self.gap_kabat.append(pieces[0])
				self.kabat.append(pieces[1])
				self.region_kabat.append(pieces[2])
				self.chothia.append(pieces[3])
				self.gap_chothia.append(pieces[4])
				self.region_chothia.append(pieces[5])
			line_num+=1
		reader.close()
		self.allowGaps=init_allowGaps


	#see if valid on the region ; making sure it's of a valid length
	def validate_region(region_info,num_aa_min,num_amino_max,allowGaps=False):
		#first verify if the region is complete.
		region_complete=True
		if(not(region_complete)):
			return False
		q_aa=cdr1_info['q_aminos']
		if(num_aa_min<=len(q_aa) and len(q_aa)<=num_amino_max):
			if(not(self.allowGaps)):
				return True
			else:
				#add code for gapfull pre-examination
				sys.exit(0)
		else:
			#too short or too long!
			return False


	#given the information on the 3 regions, verify
	#that the alignment is suitable for acquisition of mutation counts
	def validate_regions(cdr1_info,fr2_info,fr3_info):
		return False		


	#from the 3 regions, aquire a mutation map
	#POS->AA_MUT_COUNT
	def acquire_mutation_map(cdr1_info,fr2_info,fr3_info):
		return None
	




if (__name__=="__main__"):
	myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")





