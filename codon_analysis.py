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
	numberingMapCache=None
	internalAnalyzer=None
	kabat_chothia_trans=None



	#load kabat->chotia translation
	def initKabatChotiaTrans(self,tableFilePath):
		self.kabat_chothia_trans=dict()
		
		if(not(os.path.exists(tableFilePath))):
			sys.exit("ERROR, FAILED TO FIND FILE",tableFilePath," of codon numbering data!")
		with open(tableFilePath, 'r') as f:
			data=f.read()
		codon_pos_re=re.compile('^\d+[A-Z]$',re.IGNORECASE)
		lines=data.split('\n')
		for line_num in range(len(lines)):
			if(line_num!=0):
				pieces=lines[line_num].split('\t')
				kabat_num=pieces[1]
				chothia_num=pieces[3]
				codon_pos_result_kab=codon_pos_re.search(kabat_num)
				codon_pos_result_cht=codon_pos_re.search(kabat_num)
				if(codon_pos_result_kab and codon_pos_result_cht):
					self.kabat_chothia_trans[kabat_num.upper().strip()]=chothia_num.upper().strip()


	#perform KABAT->CHOTHIA numbering translation
	def kabatToChothia(self,pos):
		proper_key=pos.upper().strip()
		if(proper_key in self.kabat_chothia_trans)
			return self.kabat_chothia_trans[proper_key]
		else:
			return "X"
		
	


	def computePathToCodonTableFile(self):
		#path to script
		ownScriptPath=sys.argv[0]
		#its dir
		containingDir=os.path.dirname(ownScriptPath)
		#get relative path to codon data
		tableFilePath=containingDir+"/codon_data/codon_pos_IGHV4"
		return tableFilePath


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
		self.internalAnalyzer=CodonAnalysis()
		tableFilePath=self.computePathToCodonTableFile()
		self.initKabatChotiaTrans(tableFilePath)



	def getBPFromCodonList(self,codon_list,removeGaps=False):
		rbp=""
		for codon in codon_list:
			for bp in codon:
				if(bp=="-" and removeGaps):
					pass
				else:
					rbp+=bp
		return rbp




	#see if valid on the region ; making sure it's of a valid length
	def validate_region(self,region_info,num_aa_min,num_amino_max):
		if(self.allowGaps):
			#unimp!
			sys.exit(0)
		q_aa=region_info.getCharMap()['AA']
		q_codons=region_info.getCharMap()['nucleotide read']
		s_aa=region_info.getCharMap()['AA_ref']
		s_codons=region_info.getCharMap()['subject_read']
		join_str=""
		q_bp=self.getBPFromCodonList(q_codons)
		s_bp=self.getBPFromCodonList(s_codons)
		actual_q_aa=self.getBPFromCodonList(q_aa)
		actual_s_aa=self.getBPFromCodonList(s_aa)
		if(len(q_bp)!=len(s_bp) and not(self.allowGaps)):
			return False
		if(len(actual_q_aa)!=len(actual_s_aa)):
			return False
		if(len(q_aa)!=len(s_aa)):
			return False
		if(num_aa_min<=len(q_aa) and len(q_aa)<=num_amino_max and len(actual_q_aa)==len(q_aa)):
			return True
		else:
			#too short or too long!
			#print "len=",len(q_aa),"too short or too long"
			return False


	#given a region, get the valid lengths possible in a list format
	#lengths in CODONs
	def getRegionValidLengths(self,reg_name):
		if(reg_name=="CDR1"):
			return [5,6,7]
		elif(reg_name=="FR2" or reg_name=="FWR2"):
			return [14]
		elif(reg_name=="CDR2"):
			return [16]   #all IGHV4 sequences have length=16 for CDR2 (52A,52B,52C are NOT used in IGHV4 ; 50-65 (inclusive) are used)
		elif(reg_name=="FWR3" or reg_name=="FR3"):
			#27 if none of 82A,82B,82C are present ; 28 if 82A is present, 29 if 82A and 82B are present, and 30 if 82A,82B,and 83C are present
			return [27,28,29,30]
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
			valid_lengths=self.getRegionValidLengths(regions_to_analyze[ri])
			valid_flag=self.validate_region(region_infos[ri],min(valid_lengths),max(valid_lengths))
			valid_flags.append(valid_flag)
			if(not(valid_flag)):
				valid_on_all=False
		return valid_on_all
		

	#given a region name and length, return the numbering
	def acquireNumberingMap(self,reg_name,reg_len):
		if(self.numberingMapCache==None):
			self.numberingMapCache=dict()
		if(reg_name+str(reg_len) in self.numberingMapCache):
			return self.numberingMapCache[reg_name+str(reg_len)]
		numbering=list()
		letters=["A","B","C"]
		l_pos=0
		if(reg_name=="CDR1"):
			for p in range(31,36):
				numbering.append(p)
			for m in range(6,8):
				#lengths of 6,7
				if(reg_len>=m):
					numbering.append("35"+letters[l_pos])
					l_pos+=1
		elif(reg_name=="FR2" or reg_name=="FWR2"):
			for p in range(36,50):
				numbering.append(p)
		elif(reg_name=="CDR2"):
			for p in range(50,67):
				numbering.append(p)
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
		else:
			#invalid/unknown region!
			print "Error, unknown region ",reg_name,"!"
			sys.exit(0)
		self.numberingMapCache[reg_name+str(reg_len)]=numbering
		return self.numberingMapCache[reg_name+str(reg_len)]
			


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
			#print "\n\n\n"+reg_names[r]
			#print reg_infos[r].getName()
			#print reg_infos[r].getNiceString()
			#print "Q=",q_codons," TRX=",q_aminos
			#print "q AA len=",len(q_aminos)," q codon len=",len(q_codons)
			#print "S=",s_codons," TRX=",s_aminos
			#print "s AA len=",len(s_aminos)," s codon len=",len(s_codons)
			#print "The numbering : ",numbering_list
			for ci in range(len(q_codons)):
				if(s_codons[ci]!=q_codons[ci]):
					#mark a mutation in the numbering system
					numbered_pos=self.kabatToChothia(numbering_list[ci])
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
			hybrid_region_aln.setName("KABAT.IMGT_hybrid_FR3_"+str(vInfo['query id']))
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
			#print "EXTRACTED btop=",btop
			indel_count=getNumberIndelsFromBTOP(btop)
			#print "the count is ",indel_count
			return indel_count
		else:
			#"no btop avail"
			return -1
	else:
		#print "is none"
		return -1



#given v,d,j info maps return True if the seq should be 
#skipped due to indels
def shouldFilterOutByIndels(vInfo,dInfo,jInfo):
	if(vInfo['query id']=="HZ8R54Q02GDVXN"):
		print "gap v val=",getNumberIndelsFromBTOPInInfo(vInfo)
		print "gap d val=",getNumberIndelsFromBTOPInInfo(dInfo)
		print "gap j val=",getNumberIndelsFromBTOPInInfo(jInfo)
	if(getNumberIndelsFromBTOPInInfo(vInfo)==0 and getNumberIndelsFromBTOPInInfo(jInfo)==0 and getNumberIndelsFromBTOPInInfo(dInfo)<=0):
		shouldFilter=False
	else:
		shouldFilter=True
	return shouldFilter



def annotationMutationMap(vInfo,dInfo,jInfo,alignment_output_queue,num_submitted_jobs,imgtdb_obj,myCodonCounter):
	#perform codon mutation counting for IGHV4
	filterNote=""
	mutation_map=dict()
	mutation_map['aminos']=list()
	mutation_map['codons']=list()
	eligibleForScoring=False
	get_res=0
	kabat_CDR1=None
	kabat_FR2=None
	kabat_CDR2=None
	hybrid_FR3=None
	while(get_res<num_submitted_jobs):
		region_alignment=alignment_output_queue.get()
		if(region_alignment is not None):
			#print "\n"
			#print "THE ALIGNMENT NAME '"+region_alignment.getName()+"'"
			#print read_rec.id
			#print region_alignment.getNiceString()
			if(region_alignment.getName().startswith("CDR1_kabat")):
				kabat_CDR1=region_alignment
			elif(region_alignment.getName().startswith("FR2_kabat")):
				kabat_FR2=region_alignment
			elif(region_alignment.getName().startswith("CDR2_kabat")):
				kabat_CDR2=region_alignment
			else:
				pass
		get_res+=1

	if(vInfo is not None):
		#got V hit
		if('subject ids' in vInfo):
			#found subject in it
			if(vInfo['subject ids'].startswith("IGHV4")):
				#V hit found to be IGHV4
				#IGHV4, test regions for completeness and length
				get_res=0
				shouldFilterByIndel=shouldFilterOutByIndels(vInfo,dInfo,jInfo)
				#print "For read=",vInfo['query id']," the shouldFilterByIndel is ",shouldFilterByIndel
				if(shouldFilterByIndel):
					#print "NOTE "+read_rec.id+" filtered out by indels!"
					filterNote="Had Indels!"
					pass
				else:
					#print "NOTE "+read_rec.id+" needs completeness testing..."
					hybrid_aln=extractHybridAlignment(vInfo,imgtdb_obj)
					if(hybrid_aln==None):
						#print "couldn't get a hybrid!"
						filterNote="Failure in hybrid alignment"
					else:
						completeRegionsFlag=myCodonCounter.validate_regions_for_completenessLength(kabat_CDR1,kabat_FR2,kabat_CDR2,hybrid_aln)
						if(completeRegionsFlag):
							filterNote="OK"
							mutation_map=myCodonCounter.acquire_mutation_map(kabat_CDR1,kabat_FR2,kabat_CDR2,hybrid_aln)
							#print "THE MUTATION MAP for ",vInfo['query id']," IS ",mutation_map
						else:
							#print "INCOMPLETE so no mutation counting!"
							filterNote="Incomplete regions"
				#validate_regions_for_completenessLength(cdr1_info,fr2_info,cdr3_info,fr3_info):
			else:
				#print read_rec.id+" isn't an IGHV4 hit!"
				filterNote="Not an IGHV4 hit"
		else:
			#print read_rec.id+" has not subject ids!"
			filterNote="No name found for hit"
	else:
		#print read_rec.id+" has no vinfo"
		filterNote="NoVHit"
	return [filterNote,mutation_map]




if (__name__=="__main__"):
	myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")





