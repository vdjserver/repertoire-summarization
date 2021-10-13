"""
Repertoire calculations and comparison
"""

#
# repcalc.py
# Main function
#
# VDJServer Analysis Portal
# Repertoire calculations and comparison
# https://vdjserver.org
#
# Copyright (C) 2020 The University of Texas Southwestern Medical Center
#
# Author: Scott Christley <scott.christley@utsouthwestern.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

from __future__ import print_function
import sys
import argparse
import json
import importlib
import copy
import yaml
import yamlordereddictloader

import airr

#import vdjml
from Bio import SeqIO

# repsum modules
from repcalc import __version__
#import utils
import repcalc.defaults as defaults
#import summarize
import repcalc.metadata as metadata
import repcalc.gldb as gldb

def load_data_file(filename):
    ext = filename.split('.')[-1]
    if ext in ('yaml', 'yml'):
        with open(filename, 'r', encoding='utf-8') as handle:
            data = yaml.load(handle, Loader=yamlordereddictloader.Loader)
    elif ext == 'json':
        with open(filename, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
    else:
        if debug:
            sys.stderr.write('Unknown file type: %s. Supported file extensions are "yaml", "yml" or "json"\n' % (ext))
        raise TypeError('Unknown file type: %s. Supported file extensions are "yaml", "yml" or "json"\n' % (ext))
    return data

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
            cmod = importlib.import_module("repcalc." + defaults.calculationModules[calc]['filename'])
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

def make_parser_args():
    """Command line arguments for repcalc"""
    parser = argparse.ArgumentParser();
    parser.description='Comparison and calculation functions for immune repertoire sequencing data.'
    parser.add_argument('input',type=str,help="Input command file")
    #parser.add_argument('--gldb',type=str,nargs=1,help="Path to germline database")
    #parser.add_argument('--output',type=str,help="Output specification file")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+__version__)
    return parser

def main():
    """Perform repertoire calculations and comparison"""
    parser = make_parser_args()
    args = parser.parse_args()

    if (not args):
        args.print_help()
        sys.exit()

    # Load input specification
    inputDict = None
    try:
        with open(args.input, 'r', encoding='utf-8') as handle:
            inputDict = json.load(handle)
    except:
        print("Could not read input specification file: " + args.input)
        raise

    # Metadata
    input_file = inputDict[defaults.metadata_file_key]
    metadataDict = None
    try:
        # Load the repertoires
        data = airr.load_repertoire(input_file)

        # Put repertoires in dictionary keyed by repertoire_id
        metadataDict = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }
    except:
        print("Could not read metadata file: " + input_file)
        raise
    inputDict[defaults.full_metadata_key] = metadataDict

    # check if specific repertoires are requested
    if inputDict.get(defaults.repertoire_key) is not None:
        newDict = {}
        for rep in inputDict[defaults.repertoire_key]:
            if metadataDict.get(rep) is None:
                print("Repertoire", rep, "is missing from metadata file.")
                raise
            else:
                newDict[rep] = metadataDict[rep]
        metadataDict = newDict

    # check if repertoire groups
    if inputDict.get(defaults.groups_file_key) is not None:
        try:
            # Load the repertoire groups
            # TODO: should be able to use AIRR library to load
            #data = airr.load_repertoire(inputDict[defaults.groups_file_key])
            data = load_data_file(inputDict[defaults.groups_file_key])

            # put groups in dictionary keyed by their id
            groupDict = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

            # verify that all repertoires exist
            for group in groupDict:
                for rep in groupDict[group]['repertoires']:
                    if inputDict[defaults.full_metadata_key].get(rep['repertoire_id']) is None:
                        print("Repertoire", rep['repertoire_id'], "for group", group, "is missing from metadata file.")
                        raise

            inputDict[defaults.groups_key] = groupDict
        except:
            print("Could not read repertoire group file: " + inputDict[defaults.groups_file_key])
            raise

    # germline db
    if inputDict.get(defaults.germline_file_key) is not None:
        germline = gldb.loadGermline(inputDict[defaults.germline_file_key])
        if germline is None:
            raise
        else:
            inputDict[defaults.germline_key] = germline

    # Output specification
    #outputSpec = { "files": copy.deepcopy(inputDict[defaults.filesKey]), "groups": copy.deepcopy(inputDict[defaults.groupsKey]) }
    #print(outputSpec)

    # AIRR TSV rearrangement files
    airrFiles = inputDict[defaults.rearrangement_files_key]

    # Extract summary file list
    #summaryFiles = extract_summary_files(inputDict, metadataDict)
    #namesDict = metadata.filenames_from_uuids(metadataDict, summaryFiles)
    #print(summaryFiles)
    #print(namesDict)

    # initialize calculations
    # load and initialize calculation modules
    headerMapping = None
    initialize_calculations(inputDict, metadataDict, headerMapping)

    # iterate through each rearrangement file and perform calculations
    for airrFile in airrFiles:
        print('Processing file: ' + airrFile)
        try:
            # the calculations to be performed
            calcs = inputDict[defaults.calculationsKey]

            # Create an iteratable for rearrangement data
            reader = airr.read_rearrangement(airrFile)
            for row in reader:
                for calc in calcs:
                    cmod = defaults.calculationModules[calc['type']]['module']
                    cmod.process_record(inputDict, metadataDict, airrFile, calc, row)
        except:
            print("ERROR: Could not process rearrangement file: " + airrFile)
            raise

    # finalize calculations
    outputSpec = None
    finalize_calculations(inputDict, metadataDict, outputSpec)

#    if (args.output):
#        with open(args.output, 'w') as output_file:
#            json.dump(outputSpec, output_file, indent=2)
