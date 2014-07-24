#!/usr/bin/env python

from utils import *
from alignment import *
from char_utils import *
from segment_utils import getVRegionStartAndStopGivenRefData

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
	def validate_region(self,region_info,num_aa_min,num_amino_max):
		#first verify if the region is complete.
		#region_complete=True
		#if(not(region_complete)):
		#	return False
		q_aa=region_info.getCharMap()['AA']
		print "got ",q_aa
		if(num_aa_min<=len(q_aa) and len(q_aa)<=num_amino_max):
			if(self.allowGaps):
				#unimplemented
				#add code for gapfull pre-examination
				print "unimp"
				sys.exit(0)
			else:
				#print "TODO!"
				return True
		else:
			#too short or too long!
			print "len=",len(q_aa),"too short or too long"
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
			#return [30,31,32]
			return [30]
		else:
			#invalid region
			print "INVALID REGION PASSED "+reg_name
			sys.exit(0)


	



	#given the information on the 4 regions (CDR1,FR2,CDR2,FR3), verify
	#that the alignment is suitable for acquisition of mutation counts
	def validate_regions_for_completenessLength(self,cdr1_info,fr2_info,cdr2_info,fr3_info):
		region_infos=list()
		region_infos.append(cdr1_info)
		region_infos.append(fr2_info)
		region_infos.append(cdr2_info)
		region_infos.append(fr3_info)
		regions_to_analyze=["CDR1","FR2","CDR2","FR3"]
		valid_flags=list()
		valid_on_all=True
		for ri in range(len(regions_to_analyze)):
			valid_lens=self.getRegionValidLengths(regions_to_analyze[ri])
			valid_flag=self.validate_region(region_infos[ri],min(valid_lens),max(valid_lens))
			valid_flags.append(valid_flag)
			if(not(valid_flag)):
				valid_on_all=False
		return valid_on_all
		

	#given a region name and length, return the numbering
	def acquireNumberingMap(self,reg_name,reg_len):
		numbering=list()
		letters=["A","B","C"]
		l_pos=0
		if(reg_name=="CDR1"):
			for p in range(31,36):
				numbering.append(p)
			for m in range(6,8):
				if(reg_len>=m):
					numbering.append("35"+letters[l_pos])
					l_pos+=1
			return numbering
		elif(reg_name=="FR2" or reg_name=="FWR2"):
			for p in range(36,50):
				numbering.append(p)
			return numbering
		elif(reg_name=="CDR2"):
			for p in range(50,53):
				numbering.append(p)
			for m in range(17,20):
				if(reg_len>=m):
					numbering.append("52"+letters[l_pos])
					l_pos+=1
			for p in range(53,66):
				numbering.append(p)
			return numbering
		elif(reg_name=="FR3" or reg_name=="FWR3"):
			for p in range(66,83):
				numbering.append(p)
			for m in range(27,30):
				if(reg_len>=m):
					numbering.append("82"+letters[l_pos])
					l_pos+=1
			p=83
			while(len(numbering)!=reg_len):
				numbering.append(p)
				p+=1
			return numbering
		else:
			#invalid/unknown region!
			print "Error, unknown region ",reg_name,"!"
			sys.exit(0)
			


	#from the 3 regions, aquire a mutation map
	def acquire_mutation_map(self,cdr1_info,fr2_info,cdr2_info,fr3_info):
		reg_names=["CDR1","FR2","CDR2","FR3"]
		reg_infos=[cdr1_info,fr2_info,cdr2_info,fr3_info]
		AA_map=list()
		codon_map=list()
		for r in range(len(reg_names)):
			reg_info_map=reg_infos[r].getCharMap()
			numbering_list=self.acquireNumberingMap(reg_names[r],len(reg_info_map['AA']))
			q_codons=reg_info_map['nucleotide read']
			s_codons=reg_info_map['subject_read']
			q_aminos=reg_info_map['AA']
			s_aminos=reg_info_map['AA_ref']
			for ci in range(len(q_codons)):
				if(s_codons[ci]!=q_codons[ci]):
					#mark a mutation in the numbering system
					numbered_pos=numbering_list[ci]
					aaP=s_aminos[ci]+str(numbered_pos)+q_aminos[ci]
					cdP=s_codons[ci]+str(numbered_pos)+q_codons[ci]
					AA_map.append(aaP)
					codon_map.append(cdP)
				else:
					pass
		overall_map=dict()
		overall_map['codons']=codon_map
		overall_map['aminos']=AA_map
		return overall_map
			


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



#return a hybrid alignment
def extractHybridAlignment(vInfo,imgtdb_obj):
	v_qry_aln=vInfo['query seq']
	v_sbc_aln=vInfo['subject seq']
	q_from=int(vInfo['q. start'])
	q_to=int(vInfo['q. end'])
	s_from=int(vInfo['s. start'])
	s_to=int(vInfo['s. end'])
	v_aln=alignment(v_qry_aln,v_sbc_aln,q_from,q_to,s_from,s_to)
	hybrid_interval=getHumanHybridInterval(imgtdb_obj,vInfo['subject ids'])
	if(len(hybrid_interval)==2):
		if(hybrid_interval[0]!=(-1) and hybrid_interval[1]!=(-1)):
			hybrid_region_aln=v_aln.getSubAlnInc(hybrid_interval[0],hybrid_interval[1],"subject")
			init_frame=0 #since no gaps assume first bp has frame 0
			hybrid_region_aln.setSFM(init_frame)
			hybrid_region_aln.setName("KABAT.IMGT_hybrid_FR3")
			hybrid_region_aln.characterize()
			return hybrid_region_aln
		else:
			return None
	else:
		return None






#get indel count from an info map
#return n>=0 if a btop found
#return -1 if no info or no btop found
def getNumberIndelsFromBTOPInInfo(info):
	if(not(info==None)):
		if('btop' in info):
			btop=info['btop']
			print "EXTRACTED btop=",btop
			indel_count=getNumberIndelsFromBTOP(btop)
			#print "the count is ",indel_count
			return indel_count
		else:
			"no btop avail"
			return -1
	else:
		print "is none"
		return -1



#given v,d,j info maps return True if the seq should be 
#skipped due to indels
def shouldFilterOutByIndels(vInfo,dInfo,jInfo):
	#print "gap v val=",getNumberIndelsFromBTOPInInfo(vInfo)
	#print "gap d val=",getNumberIndelsFromBTOPInInfo(dInfo)
	#print "gap j val=",getNumberIndelsFromBTOPInInfo(jInfo)
	if(getNumberIndelsFromBTOPInInfo(vInfo)==0 and getNumberIndelsFromBTOPInInfo(jInfo)==0 and getNumberIndelsFromBTOPInInfo(dInfo)<=0):
		shouldFilter=False
	else:
		shouldFilter=True
	return shouldFilter





if (__name__=="__main__"):
	myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")





