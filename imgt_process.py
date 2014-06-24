#!/usr/bin/env python


from utils import *
from pprint import pprint
from imgt_utils import imgt_db
from cdr3_hist import CDR3LengthAnalysisVDMLOBJ,histoMapClass
from vdjml_utils import getTopVDJItems,getRegionsObjsFromSegmentCombo,getHitInfo,getProductiveRearrangmentFlag,getVDJServerRegionAlignmentFromLargerVAlignmentPyObj,getAlignmentString,getTypeNameScoreMap,getVDJTieMap
from Bio import SeqIO
from segment_utils import IncrementMapWrapper,getVRegionsList,getRegionSpecifcCharacterization,getEmptyRegCharMap,getAdjustedCDR3StartFromRefDirSetAllele,getTheFrameForThisReferenceAtThisPosition,getCDR3RegionSpecificCharacterizationSubAln,getVRegionStartAndStopGivenRefData,getADJCDR3EndFromJAllele,alleleIsTR
from char_utils import getNumberBaseSubsFromBTOP,getNumberIndelsFromBTOP,getIndelMapFromBTOP,getRegPosFromInvertedPos
from alignment import alignment
from CharacterizationThread import CharacterizationThread
from igblast_parse import rev_comp_dna
import re
import Queue
import threading
import time
import multiprocessing
import math






init_db_base="/home/data/DATABASE/06_02_2014/"
imgtdb_obj=imgt_db(init_db_base)
alleles=imgtdb_obj.getAlleles()
print "THE ALLELES ARE ",alleles
genes=imgtdb_obj.getGenes(alleles)
print "THE GENES ARE : ",genes
gene_alleles_map=dict()
for gene in genes:
	for allele in alleles:
		if(allele.startswith(gene)):
			if(gene in gene_alleles_map):
				gene_alleles_map[gene].append(allele)
				gene_alleles_map[gene].sort()
			else:
				gene_alleles_map[gene]=list()
				gene_alleles_map[gene].append(allele)

regions=getVRegionsList(True)
#print "regions : ",regions
for gene in genes:
	print "\n\n\n***********************************************************************\n***********************************************************************\n***********************************************************************"
	print "Analyzing regions for gene : ",gene
	gene_alleles=gene_alleles_map[gene]
	print "\tTo use alleles : ",gene_alleles
	for region in regions:
		print "\t\tNow analyzing for ",region
		for gene_allele in gene_alleles:
			reg_int=imgtdb_obj.getRegionStartStopFromIMGTDat(gene_allele,"human",region)
			print "For region ",region," of ",gene_allele," got ",reg_int
	print "\n\n"







#/home/data/DATABASE/06_02_2014/human/ReferenceDirectorySet/HUMAN_REF/IMGT_HighV-QUEST_individual_files_folder/IGLV9-49_01_461




