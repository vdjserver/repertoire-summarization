"""
Repertoire calculations and comparison
"""

from __future__ import print_function
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

def isRecordFunctional(headerMapping, fields):
    if  fields[headerMapping[defaults.headerNames['NonFuncCondition1']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition2']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition3']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition4']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition5']]] == 'True' or \
        fields[headerMapping[defaults.headerNames['NonFuncCondition6']]] == 'True':
            return False
    return True

def attach_records(records, name, rel_name, obj_dict, uuid_name):
    for entry in records:
        cp = entry[name]
        if cp.get(uuid_name):
            obj = obj_dict.get(cp[uuid_name]).get('value')
            obj['uuid'] = cp.get(uuid_name)
            entry[rel_name] = obj

def flatten_study_metadata(studyMetadata):
    records = []
    if not studyMetadata.get('nucleicAcidProcessingMetadata'): return records
    sdict = studyMetadata['nucleicAcidProcessingMetadata']
    for r in sdict.keys():
        obj = {'nucleicAcidProcessing': sdict[r].get('value')}
        obj['nucleicAcidProcessing']['uuid'] = r
        records.append(obj)

    attach_records(records, 'nucleicAcidProcessing', 'cellProcessing', studyMetadata.get('cellProcessingMetadata'), 'cell_processing_uuid')
    attach_records(records, 'cellProcessing', 'sample', studyMetadata.get('sampleMetadata'), 'sample_uuid')
    attach_records(records, 'sample', 'subject', studyMetadata.get('subjectMetadata'), 'subject_uuid')

    results = {}
    for r in records:
        results[r['nucleicAcidProcessing']['uuid']] = r

    return results

def restrict_by_logical(sampleList, flatMetadata, sampleGroup):
    newList = []
    for sampleKey in sampleList:
        sample = flatMetadata[sampleKey]
        # get the field value
        logicalSplit = sampleGroup['value']['logical_field'].split('.')
        if sample.get(logicalSplit[0]):
            fieldValue = sample[logicalSplit[0]].get(logicalSplit[1])
        else:
            print('WARNING: sample group (', sampleGroup['uuid'],
                  ') has unknown logical field (', sampleGroup['value']['logical_field'], ')')
            continue
        # do the check
        #print('logical: ', fieldValue)
        result = False
        if sampleGroup['value']['logical_operator'] == '=': result = (fieldValue == sampleGroup['value']['logical_value'])
        if sampleGroup['value']['logical_operator'] == '!=': result = (fieldValue != sampleGroup['value']['logical_value'])
        if sampleGroup['value']['logical_operator'] == '<': result = (fieldValue < sampleGroup['value']['logical_value'])
        if sampleGroup['value']['logical_operator'] == '>': result = (fieldValue > sampleGroup['value']['logical_value'])
        if sampleGroup['value']['logical_operator'] == '>=': result = (fieldValue >= sampleGroup['value']['logical_value'])
        if sampleGroup['value']['logical_operator'] == '<=': result = (fieldValue <= sampleGroup['value']['logical_value'])
        if sampleGroup['value']['logical_operator'] == 'contains': result = (fieldValue.find(sampleGroup['value']['logical_value']) >= 0)
        # passed the test
        if result: newList.append(sampleKey)
    return newList

def group_by_category(sampleList, flatMetadata, sampleGroup):
    categorySplit = sampleGroup['value']['category'].split('.')
    groupBy = {}
    for sampleKey in sampleList:
        sample = flatMetadata[sampleKey]
        if sample.get(categorySplit[0]):
            fieldValue = sample[categorySplit[0]].get(categorySplit[1])
        else:
            print('WARNING: sample group (', sampleGroup['uuid'],
                  ') has unknown category field (', sampleGroup['value']['category'], ')')
            continue
        #print('category: ', fieldValue)
        if not groupBy.get(fieldValue): groupBy[fieldValue] = [ ]
        group = groupBy[fieldValue]
        group.append(sampleKey)
    return groupBy

# generate map of original source files
def init_derivation_map(studyMetadata, jobMetadata):
    derivationMap = {}
    for file in jobMetadata['value']['files']:
        # summary file
        jobFile = jobMetadata['value']['files'][file].get('summary')
        if jobFile and jobFile.get('WasDerivedFrom') and jobFile['WasDerivedFrom'] is not None:
            sourceFile = was_derived_from(studyMetadata, jobMetadata, jobFile['WasDerivedFrom'])
            if sourceFile: derivationMap[jobFile['value']] = sourceFile['uuid']
    #print(derivationMap)
    return derivationMap

def was_derived_from(studyMetadata, processMetadata, fileId):
    fileMetadata = studyMetadata['fileMetadata'].get(fileId)
    if not fileMetadata:
        #print(fileId)
        # could be reference to group
        # TODO: this is specialized for the single use case, file merging in vdjpipe and presto
        if processMetadata['value']['groups'].get(fileId):
            if processMetadata['value']['groups'][fileId].get(processMetadata['value']['process']['appName']):
                file = processMetadata['value']['groups'][fileId][processMetadata['value']['process']['appName']]['files']
                fileEntry = processMetadata['value']['files'].get(file)
                if fileEntry and len(fileEntry.keys()) == 1:
                    fileType = fileEntry.keys()[0]
                    jobFile = processMetadata['value']['files'][file][fileType]
                    if jobFile.get('WasDerivedFrom') and jobFile['WasDerivedFrom'] is not None:
                        return was_derived_from(studyMetadata, processMetadata, jobFile['WasDerivedFrom'])
    elif fileMetadata['name'] == 'projectFile':
        # found original source
        return fileMetadata
    elif fileMetadata['name'] == 'projectJobFile':
        # found job file, continue up WasDerivedFrom links
        jobId = fileMetadata['value']['jobUuid']
        jobMetadata = studyMetadata['processMetadata'][jobId]
        for file in jobMetadata['value']['files']:
            for fileType in jobMetadata['value']['files'][file]:
                jobFile = jobMetadata['value']['files'][file][fileType]
                if jobFile['value'] == fileMetadata['value']['name']:
                    if jobFile.get('WasDerivedFrom') and jobFile['WasDerivedFrom'] is not None:
                        return was_derived_from(studyMetadata, jobMetadata, jobFile['WasDerivedFrom'])
    return None

def file_list_for_sample(config, studyMetadata, derivationMap, projectFile):
    fileList = []
    for file in derivationMap:
        if derivationMap[file] == projectFile:
            for groupFile in config['files']:
                groupSummary = config['files'][groupFile].get('summary')
                if not groupSummary: continue
                summaryMetadata = studyMetadata['fileMetadata'][groupSummary['value']]
                if file == summaryMetadata['value']['name']:
                    fileList.append(groupFile)
    return fileList

def create_config():
    """Generate JSON configuration file for RepCalc"""
    template = { "groups": {}, "files": {}, "calculations": [] }

    parser = argparse.ArgumentParser(description='Generate RepCalc config.')
    parser.add_argument('--init', type=str, nargs=3, help='Create initial config with metadata file', metavar=('metadata file', 'organism', 'seqtype'))
    parser.add_argument('json_file', type=str, help='Output JSON file name')

    parser.add_argument('--geneSegmentLevels', type=str, nargs='*', help='Gene segment levels')
    #parser.add_argument('--geneSegmentSummarizeBy', type=str, nargs='*', help='Gene segment summarize by')
    parser.add_argument('--geneSegmentOperations', type=str, nargs='*', help='Gene segment operations')
    parser.add_argument('--geneSegmentFilters', type=str, nargs='*', help='Gene segment filters')

    parser.add_argument('--cdr3Levels', type=str, nargs='*', help='CDR3 levels')
    parser.add_argument('--cdr3Operations', type=str, nargs='*', help='CDR3 operations')
    parser.add_argument('--cdr3Filters', type=str, nargs='*', help='CDR3 filters')

    parser.add_argument('--cdr3Single', action='store_true', help='Generate CDR3 single group config files')
    parser.add_argument('--cdr3CompareFiles', action='store_true', help='Generate CDR3 pairwise comparison config for files')
    parser.add_argument('--cdr3CompareSamples', action='store_true', help='Generate CDR3 pairwise comparison config for samples')
    parser.add_argument('--cdr3CompareGroups', action='store_true', help='Generate CDR3 pairwise comparison config for sample groups')

    parser.add_argument('--file', type=str, nargs=4, help='File set', metavar=('filekey', 'vdjml', 'summary', 'changeo'))

    parser.add_argument('--samples', type=str, help='Add samples', metavar=('jobId'))

    args = parser.parse_args()
    if args:
        if args.init:
            # save the json
            template['metadata'] = args.init[0];
            template['organism'] = args.init[1];
            template['seqType'] = args.init[2];
            with open(args.json_file, 'w') as json_file:
                json.dump(template, json_file)

        # load json
        with open(args.json_file, 'r') as f:
            config = json.load(f)

        # gene segment usage
        if args.geneSegmentOperations:
            if not args.geneSegmentLevels:
                config['calculations'].append({"type":"gene segment usage",
                                               "levels":[],
                                               "operations":args.geneSegmentOperations,
                                               "filters":args.geneSegmentFilters});
            else:
                config['calculations'].append({"type":"gene segment usage",
                                               "levels":args.geneSegmentLevels,
                                               "operations":args.geneSegmentOperations,
                                               "filters":args.geneSegmentFilters});

        # CDR3
        if args.cdr3Levels:
            if not args.cdr3Operations:
                print ("ERROR: cdr3Operations not specified.", file=sys.stderr);
                sys.exit(1);

            config['calculations'].append({"type":"CDR3",
                                           "levels":args.cdr3Levels,
                                           "operations":args.cdr3Operations,
                                           "filters":args.cdr3Filters});

        # files
        if args.file:
            config['files'][args.file[0]] = { "vdjml": { "value": args.file[1] },
                                              "summary": { "value": args.file[2] },
                                              "changeo": { "value": args.file[3] } }
            config['groups'][args.file[0]] = { "type": "file", "files": [ args.file[0] ] }

        # samples and sample groups
        if args.samples:
            with open(config['metadata'], 'r') as f:
                studyMetadata = json.load(f)

            jobMetadata = studyMetadata['processMetadata'][args.samples]
            derivationMap = init_derivation_map(studyMetadata, jobMetadata)
            flatMetadata = flatten_study_metadata(studyMetadata)
            #print(derivationMap)
            #print(jobMetadata)
            #print(json.dumps(flatMetadata, indent=2))

            # samples
            sampleList = []
            sampleFiles = {}
            sampleCount = 0
            for key in flatMetadata:
                projectFile = flatMetadata[key].get('nucleicAcidProcessing').get('filename_uuid')
                if not projectFile: continue

                # we can only include the samples for which we have data files in the job
                fileList = file_list_for_sample(config, studyMetadata, derivationMap, projectFile)
                #print(fileList)
                if len(fileList) > 0:
                    sampleName = 'sample' + str(sampleCount)
                    sampleList.append(key)
                    sampleFiles[key] = fileList
                    sampleCount += 1
                    config['groups'][sampleName] = { "type": "sample", "samples": { key: fileList } }

            # sample groups
            #print(sampleList)
            sampleSet = set(sampleList)
            groupCount = 0
            for key in studyMetadata['sampleGroups']:
                sampleGroup = studyMetadata['sampleGroups'][key]
                # determine base sample list
                if len(sampleGroup['value']['samples']) > 0: sList = sampleGroup['value']['samples']
                else: sList = sampleList
                # restrict by logical
                if len(sampleGroup['value']['logical_field']) > 0:
                    sList = restrict_by_logical(sList, flatMetadata, sampleGroup)
                # exclude samples not available in input
                groupSet = set(sList)
                iSet = groupSet & sampleSet
                if len(iSet) != len(groupSet):
                    print('WARNING: Not all sample data available for sample group: ', key)
                    sList = list(iSet)
                # group by category
                groupBy = None
                if len(sampleGroup['value']['category']) > 0:
                    groupBy = group_by_category(sList, flatMetadata, sampleGroup)
                    for groupKey in groupBy:
                        groupName = 'sampleGroup' + str(groupCount)
                        #print(groupName)
                        for aKey in groupBy[groupKey]:
                            anEntry = config['groups'].get(groupName)
                            if not anEntry: config['groups'][groupName] = { "type": "sampleGroup",
                                                                            "sampleGroup": sampleGroup['uuid'],
                                                                            "category": groupKey,
                                                                            "samples": { aKey: sampleFiles[aKey] } }
                            else: anEntry['samples'][aKey] = sampleFiles[aKey]
                        groupCount += 1
                else:
                    groupName = 'sampleGroup' + str(groupCount)
                    #print(groupName)
                    for aKey in sList:
                        anEntry = config['groups'].get(groupName)
                        if not anEntry: config['groups'][groupName] = { "type": "sampleGroup", "sampleGroup": sampleGroup['uuid'], "samples": { aKey: sampleFiles[aKey] } }
                        else: anEntry['samples'][aKey] = sampleFiles[aKey]
                    groupCount += 1

        # Separate all the groups into individual files so they can be run in parallel
        if args.cdr3Single:
            for key in config['groups']:
                customConfig = copy.deepcopy(config)
                customConfig['groups'] = {}
                customConfig['groups'][key] = config['groups'][key]
                output_json = key + '_cdr3_' + args.json_file
                with open(output_json, 'w') as json_file:
                    json.dump(customConfig, json_file, indent=2)
                # not for debugging, produce list of scripts for calling program
                print(output_json)

        #
        # Create pairwise comparison
        #
        if args.cdr3CompareFiles:
            # files
            keys = config['groups'].keys()
            for key1 in keys:
                if config['groups'][key1]['type'] != 'file': continue
                idx = keys.index(key1)
                for key2 in keys[idx:]:
                    if key1 == key2: continue
                    if config['groups'][key2]['type'] != 'file': continue
                    customConfig = copy.deepcopy(config)
                    customConfig['groups'] = {}
                    customConfig['groups'][key1] = config['groups'][key1]
                    customConfig['groups'][key2] = config['groups'][key2]
                    output_json = key1 + '_' + key2 + '_cdr3_' + args.json_file
                    with open(output_json, 'w') as json_file:
                        json.dump(customConfig, json_file, indent=2)
                    # not for debugging, produce list of scripts for calling program
                    print(output_json)

        if args.cdr3CompareSamples:
            # samples
            keys = config['groups'].keys()
            for key1 in keys:
                if config['groups'][key1]['type'] != 'sample': continue
                idx = keys.index(key1)
                for key2 in keys[idx:]:
                    if key1 == key2: continue
                    if config['groups'][key2]['type'] != 'sample': continue
                    customConfig = copy.deepcopy(config)
                    customConfig['groups'] = {}
                    customConfig['groups'][key1] = config['groups'][key1]
                    customConfig['groups'][key2] = config['groups'][key2]
                    output_json = key1 + '_' + key2 + '_cdr3_' + args.json_file
                    with open(output_json, 'w') as json_file:
                        json.dump(customConfig, json_file, indent=2)
                    # not for debugging, produce list of scripts for calling program
                    print(output_json)

        if args.cdr3CompareGroups:
            # sample groups
            keys = config['groups'].keys()
            for key1 in keys:
                if config['groups'][key1]['type'] != 'sampleGroup': continue
                idx = keys.index(key1)
                for key2 in keys[idx:]:
                    if key1 == key2: continue
                    if config['groups'][key2]['type'] != 'sampleGroup': continue
                    customConfig = copy.deepcopy(config)
                    customConfig['groups'] = {}
                    customConfig['groups'][key1] = config['groups'][key1]
                    customConfig['groups'][key2] = config['groups'][key2]
                    output_json = key1 + '_' + key2 + '_cdr3_' + args.json_file
                    with open(output_json, 'w') as json_file:
                        json.dump(customConfig, json_file, indent=2)
                    # not for debugging, produce list of scripts for calling program
                    print(output_json)

        # save the json
        with open(args.json_file, 'w') as json_file:
            json.dump(config, json_file, indent=2)

    else:
        # invalid arguments
        parser.print_help()


def make_parser_args():
	"""Command line arguments for repcalc"""
	parser = argparse.ArgumentParser();
	parser.description='Comparison and calculation functions for immune repertoire sequencing data.'
	parser.add_argument('input',type=str,nargs=1,help="Input specification file")
        parser.add_argument('--gldb',type=str,nargs=1,help="Path to germline database")
        parser.add_argument('--output',type=str,help="Output specification file")
	parser.add_argument('-v', '--version', action='version', version='%(prog)s '+__version__)
	return parser

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
                #print('Processing file: ' + input_file)
                groupSet = groups_for_file(inputDict, defaults.summaryKey, sfile)
                #print(groupSet)
                try:
                        infile = open(input_file, 'rt')
                        header = infile.readline().rstrip('\n')
                        headerMapping = summarize.summary_file_header_mappings(header)
                        if first:
                                initialize_calculations(inputDict, metadataDict, headerMapping)
                                first = False
                        calcs = inputDict[defaults.calculationsKey]
                        #print(headerMapping)
                        while True:
                                line = infile.readline().rstrip('\n')
                                if not line: break
                                fields = line.split('\t')
                                if len(fields) != len(headerMapping):
                                        print("ERROR: Number of fields does not equal number of column headings.")
                                        sys.exit()

                                for calc in calcs:
                                        if defaults.calcFilters in calc and calc[defaults.calcFilters] is not None and 'productive' in calc[defaults.calcFilters] and not isRecordFunctional(headerMapping, fields): continue
                                        cmod = defaults.calculationModules[calc['type']]['module']
                                        cmod.process_record(inputDict, metadataDict, sfile, headerMapping, groupSet, calc, fields)
                        
                except:
                        print("ERROR: Could not process summary file: " + input_file)
                        raise
                else:
                        infile.close()

        finalize_calculations(inputDict, metadataDict, outputSpec)

	if (args.output):
		with open(args.output, 'w') as output_file:
			json.dump(outputSpec, output_file, indent=2)
