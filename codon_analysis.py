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


	#init
	def __init__(self,pos_file_path):
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
				self.chothia.append(region[3])
				self.gap_chothia.append(region[4])
				self.region_chothia.append(region[5])
			line_num+=1


	#see if valid on the region
	def validate_region(region_info,num_aa_min,num_amino_max,allowGaps=False):
		q_aa=cdr1_info['q_aminos'])
		if(num_aa_min<=len(q_aa) and len(q_aa)<=num_amino_max):
			if(allowGaps):
				return True
			else:
				#add code for gap pre-examination
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





