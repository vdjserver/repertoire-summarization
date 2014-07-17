#!/usr/bin/env python

from utils import *




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

	




if (__name__=="__main__"):
	myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")





