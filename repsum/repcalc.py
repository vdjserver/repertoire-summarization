"""
Repertoire calculations and comparison
"""

import sys
import argparse
import json
import importlib
import copy

import vdjml
from Bio import SeqIO

# repsum modules
from .version import __version__
import utils
import defaults
import summarize
import metadata
import gldb

def calc_types(inputDict):
        """Return set of calculations types in input specification"""
        calcTypes = set()
        calcs = inputDict[defaults.calculationsKey]
        for calc in calcs: calcTypes.add(calc[defaults.calcTypeKey])
        return calcTypes

def initialize_calculations(inputDict, metadataDict, headerMapping):
        """Initialize each calculation module"""
        calcTypes = calc_types(inputDict)
        for calc in calcTypes:
                m = defaults.calculationModules.get(calc)
                if m is None:
                        print("ERROR: Unknown calculation module type: " + calc)
                        sys.exit()
                try:
                        cmod = importlib.import_module("repsum." + defaults.calculationModules[calc]['filename'])
                except:
                        print("ERROR: Could not load calculation module type: " + calc + "(" + defaults.calculationModules[calc]['filename'] + ".py)")
                        raise
                defaults.calculationModules[calc]['module'] = cmod
                cmod.initialize_calculation_module(inputDict, metadataDict, headerMapping)

def finalize_calculations(inputDict, metadataDict, outputSpec):
        """Finalize and save the calculations"""
        calcs = inputDict[defaults.calculationsKey]
        for calc in calcs:
             cmod = defaults.calculationModules[calc['type']]['module']
             if cmod is not None: cmod.finalize_calculation_module(inputDict, metadataDict, outputSpec, calc)

def extract_summary_files(inputDict, metadataDict):
	"""Extract set of summary files for given input and metadata"""
	summaryFiles = set()
	files = inputDict[defaults.filesKey]
	for file in files:
		#print(file)
		summaryFiles.add(files[file][defaults.summaryKey]['value'])
		#print(files[file][defaults.summaryKey]['value'])
	return summaryFiles

def extract_vdjml_files(inputDict, metadataDict):
	"""Extract set of summary files for given input and metadata"""
	summaryFiles = set()
	files = inputDict[defaults.filesKey]
	for file in files:
		summaryFiles.add(files[file][defaults.vdjmlKey]['value'])
	return summaryFiles

def groups_for_file(inputDict, fileKey, uuid):
        """Return set of applicable groups for a given summary file"""
        groupSet = set()
        files = inputDict[defaults.filesKey]
        groups = inputDict[defaults.groupsKey]
        for group in groups:
		if (groups[group]['type'] == 'file'):
			for f in groups[group]['files']:
				file = files[f]
				if file[fileKey]['value'] == uuid: groupSet.add(group)
		if (groups[group]['type'] == 'sample'):
			for sample in groups[group]['samples']:
                                for f in groups[group]['samples'][sample]:
                                        file = files[f]
                                        if file[fileKey]['value'] == uuid: groupSet.add(group)
		if (groups[group]['type'] == 'sampleGroup'):
			for sample in groups[group]['samples']:
                                for f in groups[group]['samples'][sample]:
                                        file = files[f]
                                        if file[fileKey]['value'] == uuid: groupSet.add(group)
        return groupSet

def make_parser_args():
	"""Command line arguments for repcalc"""
	parser = argparse.ArgumentParser();
	parser.description='Comparison and calculation functions for immune repertoire sequencing data.'
	parser.add_argument('input',type=str,nargs=1,help="Input specification file")
        parser.add_argument('--gldb',type=str,nargs=1,help="Path to germline database")
        parser.add_argument('--output',type=str,help="Output specification file")
	parser.add_argument('-v', '--version', action='version', version='%(prog)s '+__version__)
	return parser

def isRecordFunctional(headerMapping, fields):
    if  fields[headerMapping[defaults.headerNames['NonFuncCondition1']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition2']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition3']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition4']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition5']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition6']]] == 'True':
            return False
    return True

def main():
	"""Perform repertoire calculations and comparison"""
	parser = make_parser_args()
	args = parser.parse_args()

	if (not args):
		args.print_help()
		sys.exit()

        # germline db
        gldb_path = utils.extractAsItemOrFirstFromList(args.gldb)
        gldb.init_germline_db_root(gldb_path)
        #print(gldb.germline_db_root())

	# Load input specification
	input_file = utils.extractAsItemOrFirstFromList(args.input)
	try:
                infile = open(input_file, 'rt')
                inputDict = json.load(infile)
	except:
                print("Could not read input specification file: " + input_file)
                raise
	else:
                infile.close()
        if inputDict.get(defaults.organismKey) is None:
                print("WARNING: No organism defined, assuming human")
                inputDict[defaults.organismKey] = 'human'
	#print(json.dumps(inputDict))

	# Metadata
	#print(inputDict[defaults.metadata_file_key])
	input_file = inputDict[defaults.metadata_file_key]
	try:
                infile = open(input_file, 'rt')
                metadataDict = json.load(infile)
	except:
                print("Could not read metadata file: " + input_file)
                raise
	else:
                infile.close()
	#print(json.dumps(metadataDict))

	# Output specification
        outputSpec = { "files": copy.deepcopy(inputDict[defaults.filesKey]), "groups": copy.deepcopy(inputDict[defaults.groupsKey]) }
        #print(outputSpec)

	# Extract summary file list
	summaryFiles = extract_summary_files(inputDict, metadataDict)
        namesDict = metadata.filenames_from_uuids(metadataDict, summaryFiles)
	#print(summaryFiles)
        #print(namesDict)

	# walk through each file and perform calculations
        first = True
        for sfile in summaryFiles:
                input_file = namesDict[sfile]
                print('Processing file: ' + input_file)
                groupSet = groups_for_file(inputDict, defaults.summaryKey, sfile)
                #print(groupSet)
                try:
                        infile = open(input_file, 'rt')
                        header = infile.readline()
                        headerMapping = summarize.summary_file_header_mappings(header)
                        if first:
                                initialize_calculations(inputDict, metadataDict, headerMapping)
                                first = False
                        calcs = inputDict[defaults.calculationsKey]
                        #print(headerMapping)
                        while True:
                                line = infile.readline()
                                if not line: break
                                fields = line.split('\t')
                                if len(fields) != len(headerMapping):
                                        print("ERROR: Number of fields does not equal number of column headings.")
                                        sys.exit()

                                for calc in calcs:
                                        if defaults.calcFilters in calc and calc[defaults.calcFilters] is not None and 'productive' in calc[defaults.calcFilters] and not isRecordFunctional(headerMapping, fields): continue
                                        cmod = defaults.calculationModules[calc['type']]['module']
                                        cmod.process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields)
                        
                except:
                        print("ERROR: Could not process summary file: " + input_file)
                        raise
                else:
                        infile.close()

        finalize_calculations(inputDict, metadataDict, outputSpec)

	if (args.output):
		with open(args.output, 'w') as output_file:
			json.dump(outputSpec, output_file, indent=2)
