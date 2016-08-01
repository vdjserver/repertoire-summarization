"""
Repertoire calculations and comparison
"""

import sys
import vdjml
import utils
from .version import __version__
import argparse
import json
import defaults
from Bio import SeqIO
import summarize
import metadata

def extract_summary_files(inputDict, metadataDict):
	"""Extract unique set of summary files for given input and metadata"""
	summaryFiles = set()
	groups = inputDict[defaults.groups_key]
	files = inputDict[defaults.files_key]
	for group in groups:
		for sample in groups[group]['samples']:
			file = files[groups[group]['samples'][sample]]
			mfile = metadata.file_with_uuid(metadataDict, file['summary'])
			if (mfile): summaryFiles.add(metadata.filename(mfile))
			print(file['summary'])

	return summaryFiles

def make_parser_args():
	"""Command line arguments for repcalc"""
	parser = argparse.ArgumentParser();
	parser.description='Comparison and calculation functions for immune repertoire sequencing data. VERSION '+__version__
	parser.add_argument('input',type=str,nargs=1,help="Input specification file")
	return parser


def main():
	"""Perform repertoire calculations and comparison"""
	parser = make_parser_args()
	args = parser.parse_args()

	if (not args):
		args.print_help()
		sys.exit()

	# Load input specification
	input_file = utils.extractAsItemOrFirstFromList(args.input)
	try:
			infile = open(input_file)
			inputDict = json.load(infile)
	except:
			print("Could not read input specification file: " + input_file)
			raise
	else:
			infile.close()
	#print(json.dumps(inputDict))

	# Metadata
	#print(inputDict[defaults.metadata_file_key])
	input_file = inputDict[defaults.metadata_file_key]
	try:
			infile = open(input_file)
			metadataDict = json.load(infile)
	except:
			print("Could not read metadata file: " + input_file)
			raise
	else:
			infile.close()
	#print(json.dumps(metadataDict))

	# Extract summary file list
	summaryFiles = extract_summary_files(inputDict, metadataDict)
	print(summaryFiles)

	# Calculations

	# Summary
        
