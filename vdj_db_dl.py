#!/usr/bin/env python

import argparse
from utils import *
from imgt_utils import imgt_db
from segment_utils import analyze_download_dir_forVDJserver
import sys
import glob
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
def downloadAndPrep(imgtdb_obj,makeblastdbbin,igblastnbin):
	#download_ReferenceDirDataHMAndGeneTables(imgtdb_obj):
	#analyze_download_dir_forVDJserver(imgtdb_obj.getBaseDir())
	#imgtdb_obj.prepareFASTAForBLASTFormatting()
	#imgtdb_obj.blastFormatFNAInRefDirSetDirs(makeblastdbbin)
	#downloadIMGTGENE_DB_and_LIGM_DB_and_index(imgtdb_obj)
	kabat_process(imgtdb_obj,igblastnbin,None,None,None,None)



#def setupKabatDir(new_dir,V_query,J_query):
	


#kabat PROCESSING 
#run igblastn on V data against the IGBLAST-provided database 
#run the kabat routines on the IGBLAST output to make the KABAT region lookup tables
#run blastx with the J nucleotides as the input/query data using IMGT reference data as the 
def kabat_process(imgtdb_obj,igblastnbin,vdb,ddb,jdb,blastxbin):
	#find  /home/data/DATABASE/01_22_2014/|grep -Pi 'ReferenceDirectorySet/.*_[A-Z]{2}_V\.fna$'
	imgtdb_base=imgtdb_obj.getBaseDir()
	organisms=imgtdb_obj.getOrganismList()
	for organism in organisms:
		print "got ",organism
		IG_V_glob=imgtdb_base+"/"+organism+"/ReferenceDirectorySet/*_IG_V.fna"
		glob_res_v=glob.glob(IG_V_glob)
		TR_V_glob=imgtdb_base+"/"+organism+"/ReferenceDirectorySet/*_TR_V.fna"
		glob_res_tr_v=glob.glob(TR_V_glob)
		if(len(glob_res_v)==1 and len(glob_res_tr_v)==1):
			IG_V=glob_res_v[0]
			TR_V=glob_res_tr_v[0]
			if(os.path.exists(IG_V)):
				print "Found fild ",IG_V
				os.mkdir(imgtdb_base+"/"+organism+"/ReferenceDirectorySet/KABAT/")
				
			else:
				print "file ",IG_V," not found!"
		else:
			print "Error in GLOBBING for files for kabat process!"
			print "GLOB",IG_V_glob," resulted in ",len(glob_res_v)," files found!"
			print "GLOB",TR_V_glob," resulted in ",len(glob_res_tr_v)," files found !"
			print "KABAT V lookups not to be available!!!!"
	



#get VDJ fasta, heavy/light flags and simulate
#with 'dumb' SHM (includes insertions, deletions, and base substitutions)
if (__name__=="__main__"):
	parser = argparse.ArgumentParser(description='Download IMGT reference data for VDJ Server. ')
	parser.add_argument('imgt_db_base',type=str,nargs=1,help="path to a NON-existent directory where downloading will take place")
	parser.add_argument('makeblastdb_bin',type=str,nargs=1,help="*full* path to the makeblastdb binary executable")
	parser.add_argument('igblast_bin',type=str,nargs=1,help="*full* path to the igblastn executable")
	parser.add_argument('blastx_bin',type=str,nargs=1,help="*full path to the blastx binary executable")
	parser.add_argument('map_file',type=str,nargs=1,help="*full* path to the map file (tab-separated key-value pairs")
	args=args = parser.parse_args()
	if(args):
		imgt_db_base=extractAsItemOrFirstFromList(args.imgt_db_base)
		makeblastdb_bin=extractAsItemOrFirstFromList(args.makeblastdb_bin)
		igblast_bin=extractAsItemOrFirstFromList(args.igblast_bin)
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
		downloadAndPrep(imgtdb_obj,makeblastdb_bin,igblast_bin)
