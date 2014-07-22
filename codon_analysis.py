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
	def validate_region(region_info,num_aa_min,num_amino_max):
		#first verify if the region is complete.
		region_complete=True
		if(not(region_complete)):
			return False
		q_aa=cdr1_info['q_aminos']
		if(num_aa_min<=len(q_aa) and len(q_aa)<=num_amino_max):
			if(self.allowGaps):
				#unimplemented
				#add code for gapfull pre-examination
				sys.exit(0)
			else:
				pass
		else:
			#too short or too long!
			return False


	#given a region, get the valid lengths possible in a list format
	#lengths in CODONs
	def getRegionValidLengths(self,reg_name):
		if(reg_name=="CDR1"):
			return [5,6,7]
		elif(reg_name=="FR2" or reg_name=="FWR2"):
			return [14]
		elif(reg_name=="CDR2"):
			return [16]
		elif(reg_name=="FWR3" or reg_name=="FR3"):
			return [32]
		else:
			#invalid region
			print "INVALID REGION PASSED "+reg_name
			sys.exit(0)


	#given the information on the 3 regions, verify
	#that the alignment is suitable for acquisition of mutation counts
	def validate_regions(cdr1_info,fr2_info,fr3_info):
		return False		


	#from the 3 regions, aquire a mutation map
	#POS->AA_MUT_COUNT
	def acquire_mutation_map(cdr1_info,fr2_info,fr3_info):
		return None
	
#given an IGHV4 allele (for human), return a hybrid
#interval with FR3 start in KABAT and FR3 end in IMGT
def getHumanHybridInterval(imgtdb_obj,ighv4allelename):
	if(imgtdb_obj==None or ighv4allelename==None):
		return [-1,-1]
	if(not(ighv4allelename.startswith("IGHV4"))):
		return [-1,-1]
	else:
		fr3_interval_kabat=getVRegionStartAndStopGivenRefData(ighv4allelename,"human",imgtdb_obj,"FR3","kabat")
		fr3_interval_imgt=getVRegionStartAndStopGivenRefData(ighv4allelename,"human",imgtdb_obj,"FR3","imgt")
		hybrid_interval=[fr3_interval_kabat[0],fr3_interval_imgt[1]]
		return hybrid_interval


#get indel count from an info map
#return n>=0 if a btop found
#return -1 if no info or no btop found
def getNumberIndelsFromBTOP(info):
	if(not(info==None)):
		if('btop' in info):
			btop=info['btop']
			indel_count=getNumberIndelsFromBTOP(btop)
			return indel_count
		else:
			return -1
	else:
		return -1



#given v,d,j info maps return True if the seq should be 
#skipped due to indels
def shouldFilterOutByIndels(vInfo,dInfo,jInfo):
	if(getNumberIndelsFromBTOP(vInfo)==0 and getNumberIndelsFromBTOP(jInfo)==0 and getNumberIndelsFromBTOP(dInfo)<=0):
		shouldFilter=False
	else:
		shouldFilter=True
	return shouldFilter





if (__name__=="__main__"):
	myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")





