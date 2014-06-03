#!/usr/bin/env python

import argparse
from utils import *
from imgt_utils import imgt_db
from segment_utils import analyze_download_dir_forVDJserver
from kabat_process import *
import sys
import glob
import os



	




def download_ReferenceDirDataHMAndGeneTables(imgtdb_obj):
	#this will download the FASTAs
	#and gene tables
	imgtdb_obj.download_imgt_RefDirSeqs_AndGeneTables_HumanAndMouse(True)

#download FASTA and ANNOTATION data
def downloadIMGTGENE_DB_and_LIGM_DB_and_index(imgtdb_obj):
	imgtdb_obj.buildAndExecuteWGETDownloadScript()
	imgtdb_obj.indexIMGTDatFile()
	imgtdb_obj.cacheIndex()


#the main wrapper for downloading and preparing the database (except for parts with necessary human intervention!)
def downloadAndPrep(imgtdb_obj,makeblastdbbin,igblastnbin,kvMap,blastx_bin):
	#download fasta and annotation data
	downloadIMGTGENE_DB_and_LIGM_DB_and_index(imgtdb_obj)
	#partition the loci (ighv,ighd, etc)
	imgtdb_obj.buildRefDirSetsFromGENEDB()
	imgtdb_obj.prepareFASTAForBLASTFormatting()
	#downloag the gene tables for the IMGT hierarchy
	imgtdb_obj.download_GeneTables()
	#compare sequences in FNA GL files and 
	analyze_download_dir_forVDJserver(imgtdb_obj.getBaseDir())
	imgtdb_obj.blastFormatFNAInRefDirSetDirs(makeblastdbbin)
	#use IgBLAST and BLASTX to get region and CDR3 information
	kabat_process(imgtdb_obj,igblastnbin,blastx_bin,kvMap,makeblastdbbin)




#kabat PROCESSING 
#run igblastn on V data against the IGBLAST-provided database 
#run the kabat routines on the IGBLAST output to make the KABAT region lookup tables
#run blastx with the J nucleotides as the input/query data using IMGT reference data
#IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP as the database
def kabat_process(imgtdb_obj,igblastnbin,blastxbin,kvMap,makeblastdbbin):
	#find  /home/data/DATABASE/01_22_2014/|grep -Pi 'ReferenceDirectorySet/.*_[A-Z]{2}_V\.fna$'
	imgtdb_base=imgtdb_obj.getBaseDir()
	organisms=imgtdb_obj.getOrganismList()
	for organism in organisms:
		#print "got ",organism
		#DO V LOOKUPS HERE : FR1,CDR1,FR2,CDR2,FR3
		stList=["IG"]
		for st in stList:
			kabat_dir=imgtdb_base+"/"+organism+"/ReferenceDirectorySet/KABAT/"
			if(not os.path.exists(kabat_dir)):
				os.mkdir(kabat_dir)
			V_glob=imgtdb_base+"/"+organism+"/ReferenceDirectorySet/*_"+st+"_V.fna"
			glob_res_v=glob.glob(V_glob)
			if(len(glob_res_v)==1):
				queryFile=glob_res_v[0]
				if(os.path.exists(queryFile)):
					print "Found query file ",queryFile
					db_base_key=None
					aux_org=None
					if(organism=="human"):
						db_base_key="HUMAN_GL_IG_"
						aux_org="HUMAN"
					elif(organism=="Mus_musculus"):
						db_base_key="MOUSE_GL_IG_"
						aux_org="MOUSE"
					else:
						print "ERROR!  UNKNOWN ORGANISM '",organism,"' or system/code needs upgrade for this organism!"
						sys.exit(0)
					igblast_cmd=igblastnbin+" -query "+queryFile+" -ig_seqtype Ig "
					for segment in get_segment_list():
						newPart=" -germline_db_"+segment+" "+kvMap[db_base_key+segment]+" "
						igblast_cmd+=newPart
					aux_base=re.sub(r'GL_IG_','AUX',db_base_key)
					igblast_cmd+=" -auxiliary_data "+kvMap[aux_org+"_AUX"]
					igblast_cmd+=" -domain_system kabat "
					igblast_cmd+=" -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'"
					igblast_output=kabat_dir+"igblast.kabat.out"
					igblast_cmd+=" -out "+igblast_output
					print "Placing KABAT lookup data for ",organism," in directory ",kabat_dir,"...."
					print "Running igblast command : ",igblast_cmd					
					runCmdButGetNoStdOutOrStdErr(igblast_cmd)
					errStat=False
					if(not(os.path.exists(igblast_output))):
						print "ERROR, it appears that IGBLAST did not run!?"
						errStat=True
					if(not(errStat)):
						igbLineCount=fileLineCount(igblast_output)
						if(igbLineCount<100):
							print "ERROR, it appears that IGBLAST ran but with errors!"
							errStat=True
						if(not(errStat)):
							Vlookup=kabat_dir+"/Vlookup.tsv"
							print "Now generating V lookup for ",organism," IG from igblast output ",igblast_output," and writing to ",Vlookup
							writeKabatRegionsFromIGBLASTKabatResult(igblast_output,Vlookup)											
				else:
					print "file ",queryFile," not found!"
			else:
				print "Error in GLOBBING for files for kabat process!"
				print "GLOB",V_glob," resulted in ",len(glob_res_v)," files found!"
				print "KABAT V lookup not to be available!!!!"
			J_glob=imgtdb_base+"/"+organism+"/ReferenceDirectorySet/*_"+st+"_J.fna"
			glob_res_j=glob.glob(J_glob)
			if(len(glob_res_j)==1):
				queryFile=glob_res_j[0]
				blastx_fna_db=imgtdb_obj.getBaseDir()+"/www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP"
				if(not(os.path.exists(blastx_fna_db))):
					print "Error, could not find ",blastx_fna_db," for J CDR3 KABAT!"
				else:
					print "Now formatting for blastx :",blastx_fna_db
					format_blast_db(blastx_fna_db,"prot",makeblastdbbin)
					blastx_xml_out=kabat_dir+"/blastx.out.xml"
					blastxcmd=blastxbin+"   -max_target_seqs  1   -query "+queryFile+" -outfmt 5 -db "+blastx_fna_db+" -out "+blastx_xml_out
					print "Now to run blastx with command ",blastxcmd
					runCmdButGetNoStdOutOrStdErr(blastxcmd)
					Jlookup=kabat_dir+"/Jlookup.tsv"
					print "\nWriting Jlookup to ",Jlookup,"\n\n"
					writeKabatJCDR3End(blastx_xml_out,Jlookup)
			else:
				print "Error in GLOBBING for files for kabat process!"
				print "GLOB",J_glob," resulted in ",len(glob_res_j)," files found!"
				print "KABAT J lookup not to be available!!!!"




if (__name__=="__main__"):
	parser = argparse.ArgumentParser(description='Download IMGT reference data for VDJ Server. ')
	parser.add_argument('imgt_db_base',type=str,nargs=1,help="path to a NON-existent directory where downloading will take place")
	parser.add_argument('makeblastdb_bin',type=str,nargs=1,help="*full* path to the makeblastdb binary executable")
	parser.add_argument('igblast_bin',type=str,nargs=1,help="*full* path to the igblastn executable")
	parser.add_argument('blastx_bin',type=str,nargs=1,help="*full* path to the blastx binary executable")
	parser.add_argument('map_file',type=str,nargs=1,help="*full* path to the map file (tab-separated key-value pairs)")
	parser.add_argument('-analyze_only',action='store_true',help="simply analyze the database comparing gene table records with FASTA reference directory records, showing before/after effect of patch files.")
	args=args = parser.parse_args()
	if(args):
		imgt_db_base=extractAsItemOrFirstFromList(args.imgt_db_base)
		analyze_flag=extractAsItemOrFirstFromList(args.analyze_only)
		if(analyze_flag):
			imgtdb_obj=imgt_db(imgt_db_base)
			analyze_download_dir_forVDJserver(imgtdb_obj.getBaseDir())
			sys.exit(0)
		makeblastdb_bin=extractAsItemOrFirstFromList(args.makeblastdb_bin)
		igblast_bin=extractAsItemOrFirstFromList(args.igblast_bin)
		blastx_bin=extractAsItemOrFirstFromList(args.blastx_bin)
		print "using base ",imgt_db_base
		mapPath=extractAsItemOrFirstFromList(args.map_file)
		if(os.path.exists(mapPath)):
			kvMap=read_map_from_file(mapPath)
			printMap(kvMap)
		else:
			print "Error, map file ",mapPath," not found!"
			sys.exit(0)
		if(os.path.isdir(imgt_db_base)):
			print "Error, directory",imgt_db_base," must not exist! Abort!"
			sys.exit(1)
		elif(os.path.exists(imgt_db_base)):
			print "Error, directory",imgt_db_base," must be a non-existent directory!"
			sys.exit(1)
		else:
			os.makedirs(imgt_db_base)
		imgtdb_obj=imgt_db(imgt_db_base)
		downloadAndPrep(imgtdb_obj,makeblastdb_bin,igblast_bin,kvMap,blastx_bin)













