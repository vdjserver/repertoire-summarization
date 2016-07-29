#!/usr/bin/env python

"""
Repertoire-summarization
"""

import vdjml
import utils
from .version import __version__
import argparse
import json
import defaults

def makeParserArgs():
	parser = argparse.ArgumentParser();
	parser.description='Summarization and comparison functions for immune repertoire sequencing data. VERSION '+__version__
	parser.add_argument('input',type=str,nargs=1,help="Input specification file")
	return parser


def main():
	parser = makeParserArgs()
	args = parser.parse_args()

	if (not args):
		args.print_help()
		sys.exit()

	input_file = utils.extractAsItemOrFirstFromList(args.input)
	infile = open(input_file)
	calcSpec = json.load(infile)
	fclose(infile)
	
	print(calcSpec[defaults.metadata_file_key])
	infile = open(calcSpec[defaults.metadata_file_key])
	
"""
		organism=extractAsItemOrFirstFromList(args.db_organism)
		#print "the file is ",args.vdjml_in[0]
                vdjml_file = extractAsItemOrFirstFromList(args.vdjml_in)
                vr = vdjml.Vdjml_reader(vdjml_file)
		#mf=makeVDJMLDefaultMetaAndFactoryFromArgs(args)
		sample_json_path=extractAsItemOrFirstFromList(args.sample_json_out)
		#meta=mf[0]
		#fact=mf[1]
                meta = vr.meta()
		imgtdb_obj=imgt_db(extractAsItemOrFirstFromList(args.vdj_db_root))
		#print "a fact is ",fact
		#print "the file is ",args.igblast_in[0]
		query_fasta=extractAsItemOrFirstFromList(args.qry_fasta)
		if(os.path.exists(query_fasta) and os.path.isfile(query_fasta)):
			#tot_reads=countFastaReads(query_fasta)
			tot_reads=None
			#do we really wanna count the reads?????
		else:
			tot_reads=None
		fasta_reader=SeqIO.parse(open(query_fasta, "r"), "fasta")
		modes=['kabat','imgt']
		my_cdr3_map=histoMapClass(modes)
		cdr3_hist_out_file=args.cdr3_hist_out
		segment_counter=IncrementMapWrapper()
		combo_counter=recombFreqManager()
		read_num=1
		segments_json_out=extractAsItemOrFirstFromList(args.json_out)
		#print "to write vdjml to ",args.vdjml_out
		#rrw = vdjml.Vdjml_writer(extractAsItemOrFirstFromList(args.vdjml_out), meta)
		rep_char_out=extractAsItemOrFirstFromList(args.char_out)
		fHandle=open(rep_char_out,'w')
		if(False):
			logFile=rep_char_out+".log"
			logHandle=open(logFile,'w')
		else:
			logHandle=None
		if(tot_reads==None):
			every_read=100
		else:
			every_read=math.log(tot_reads,10.0)
			if(every_read<=1):
				every_read=1
			else:
				every_read=int(every_read)
				every_read=int(10**every_read)
			#print "every_read at ",every_read

		#test release info
		loadReleaseInfoTagAndHash()
		if(release_info_tag==None or release_info_hash==None):
			print "ERROR IN RELEASE_INFO ?!"


                while vr.has_result():
                        read_result_obj = vr.result()
                        
		#for read_result_obj in scanOutputToVDJML(extractAsItemOrFirstFromList(args.igblast_in),fact,query_fasta):

			#prepare for the iteration and give a possible status message...
			if(read_num>1 and read_num%every_read==0):
				print "Processing read",read_num,"..."
			elif(read_num==1):
				print "Processing reads..."
			query_record=fasta_reader.next()

			#analyze a read's results
			read_analysis_results=rep_char_read(read_result_obj,meta,organism,imgtdb_obj,query_record,args.skip_char)

			#handle cdr3 length/histogram
			if('cdr3_length_results' in read_analysis_results):
				cdr3_res=read_analysis_results['cdr3_length_results']
				for mode in modes:
					my_cdr3_map.inc(mode,cdr3_res[mode])

			#handle segment counting
			if('VDJ' in read_analysis_results):
				segments=read_analysis_results['VDJ']
				for s in segments:
					actual=segments[s]
					if(actual is not None):
						segment_counter.increment(actual)
				vseg=segments['V']
				dseg=segments['D']
				jseg=segments['J']
				combo_counter.addVDJRecombination(vseg,dseg,jseg,False)

			#process for writing
			#rrw(read_result_obj)

			#write rep-char
			appendAnnToFileWithMap(fHandle,read_analysis_results['ann_map'],read_num,query_record.id,None,"None",logHandle)

			#increment the read number
			read_num+=1
                        vr.next()

		print "Processed total of",str(read_num-1),"reads!"

		#close rep_char out
		fHandle.close()

		#close log
		if(logHandle is not None):
			logHandle.close()

		#write the CDR3 hist when non-dev-null
		if(type(cdr3_hist_out_file)==list):
			cdr3_hist_out_file=extractAsItemOrFirstFromList(cdr3_hist_out_file)
		#print "file=",cdr3_hist_out_file
		if(not(cdr3_hist_out_file=="/dev/null")):
			my_cdr3_map.writeToFile(cdr3_hist_out_file)
			print "Wrote CDR3 lengths histogram to ",cdr3_hist_out_file
                print "NOTE : Number IMGT  CDR3 lengtgs not found = "+str(my_cdr3_map.count_map['imgt'].get(-1))
		print "NOTE : Number KABAT CDR3 lengths not found = "+str(my_cdr3_map.count_map['kabat'].get(-1))


		#write AGS/NMO information
		if(organism=="human"):
			print "AGS_SCORE\t"+str(myCodonCounter.computeAGS())+"\tAGS6_TOT_RM\t"+str(myCodonCounter.computeAGS6TotRM())+"\tTOT_RM\t"+str(myCodonCounter.computeSampTotRM())
			print "AGS5_SCRE\t"+str(myCodonCounter.computeAGS5())+"\tAGS5_TOT_RM\t"+str(myCodonCounter.computeAGS5TotRM())+"\tTOT_RM\t"+str(myCodonCounter.computeSampTotRM())
			print "NMO_SCORE\t"+str(myCodonCounter.computeNMO())+"\tNMO_QUERIES_RM\t"+str(myCodonCounter.queriesWithRM)+"\tNMO_SAMP_NUC_TOT\t"+str(myCodonCounter.computeSampNMOTot())
		else:
			print "No AGS/NMO data for non-human analysis!"


		print "Writing JSON segment counts output to ",segments_json_out,"..."
		segment_counter.JSONIFYToFile(
			extractAsItemOrFirstFromList(args.vdj_db_root),
			organism,
			segments_json_out,
			False,
			imgtdb_obj.getPickleFullPath()
			)
		print "Writing JSON segment counts output complete!"
		recomb_out_file=extractAsItemOrFirstFromList(args.combo_out)
		print "Writing JSONS segment combination frequency data to ",recomb_out_file
		combo_counter.writeJSONToFile(recomb_out_file)
		combo_csv_string=combo_counter.writeDataToCSV(False)
		combo_csv_file=extractAsItemOrFirstFromList(args.combo_csv_out)
		print "Writing CSV segment combination frequency data to ",combo_csv_file
		writeStringToFilePathAssumingDirectoryExists(combo_csv_string,combo_csv_file)
		print "Writing sample-level JSON data to ",sample_json_path
		samp_json=generateSampleLevelStatsJSON()
		writeStringToFilePathAssumingDirectoryExists(samp_json,sample_json_path)
"""

