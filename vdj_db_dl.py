#!/usr/bin/env python

import argparse
from utils import extractAsItemOrFirstFromList
from imgt_utils import imgt_db
from segment_utils import analyze_download_dir_forVDJserver
import sys
import os

def download_ReferenceDirDataHMAndGeneTables(imgtdb_obj):
	#this will download the FASTAs
	#and gene tables
	imgtdb_obj.download_imgt_RefDirSeqs_AndGeneTables_HumanAndMouse(True)


def downloadIMGTGENE_DB_and_LIGM_DB_and_index(imgtdb_obj):
	#imgtdb_obj.buildAndExecuteWGETDownloadScript()
	imgtdb_obj.indexIMGTDatFile()
	imgtdb_obj.cacheIndex()


#the main wrapper for downloading and preparing the database (except for parts with necessary human intervention)
def downloadAndPrep(imgtdb_obj):
	#download_ReferenceDirDataHMAndGeneTables(imgtdb_obj):
	#analyze_download_dir_forVDJserver(imgtdb_obj.getBaseDir())
	#downloadIMGTGENE_DB_and_LIGM_DB_and_index(imgtdb_obj)


#kabat PROCESSING 
#run igblastn on V data against the IGBLAST-provided database 
#run the kabat routines on the IGBLAST output to make the KABAT region lookup tables
#run blastx with the J nucleotides as the input/query data using IMGT reference data as the 




#get VDJ fasta, heavy/light flags and simulate
#with 'dumb' SHM (includes insertions, deletions, and base substitutions)
if (__name__=="__main__"):
	parser = argparse.ArgumentParser(description='Download IMGT reference data for VDJ Server. ')
	parser.add_argument('imgt_db_base',type=str,nargs=1,help="path to a NON-existent directory where downloading will take place")
	parser.add_argument('igblast_bin',type=str,nargs=1,help="*full* path to the igblastn executable")
	parser.add_argument('igblast_db_
	args=args = parser.parse_args()
	if(args):
		imgt_db_base=extractAsItemOrFirstFromList(args.imgt_db_base)
		print "using base ",imgt_db_base
		if(os.path.isdir(imgt_db_base)):
			print "Error, directory",imgt_db_base," must not exist! Abort!"
			#sys.exit(1)
		elif(os.path.exists(imgt_db_base)):
			print "Error, directory",imgt_db_base," must be a non-existent directory!"
			#sys.exit(1)
		else:
			os.makedirs(imgt_db_base)
		imgtdb_obj=imgt_db(imgt_db_base)
		downloadAndPrep(imgtdb_obj)
