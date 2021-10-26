"""
Germline database functions
"""

#
# gldb.py
# Germline sets and gene description utility functions
#
# VDJServer Analysis Portal
# Repertoire calculations and comparison
# https://vdjserver.org
#
# Copyright (C) 2021 The University of Texas Southwestern Medical Center
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

import sys
import airr
import json

def loadGermline(filename):
    """Load germline sets and allele descriptions from AIRR file"""
    try:
        with open(filename, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
    except:
        print("Could not read germline file: " + filename)
        raise
    print('Loaded germline database:', filename)

    # organize the data by IDs
    germline = {}
    germline['germline_sets'] = { obj['germline_set_id'] : obj for obj in data['GermlineSet'] }
    germline['allele_descriptions'] = { obj['allele_description_id'] : obj for obj in data['AlleleDescription'] }
    return germline

def getGermlineSet(germline_set_id):
    """Return germline set from its identifier"""
    return None

def getGene(allele_description):
    return allele_description['locus'] + allele_description['sequence_type'] + allele_description['gene_designation']

def getSubgroup(allele_description):
    if allele_description['subgroup_designation'] is None:
        return allele_description['locus'] + allele_description['sequence_type'] + allele_description['gene_designation']
    else:
        return allele_description['locus'] + allele_description['sequence_type'] + allele_description['subgroup_designation']

def getDisplayName(germline, allele_call, level):
    """Return display name for allele call"""
    if germline is None:
        print('Required germline database is missing.')
        sys.exit(1)

    # handle multiple allele calls
    fields = allele_call.split(',')
    distinct_names = []
    for field in fields:
        allele_description = germline['allele_descriptions'].get(field)
        if allele_description is None:
            print('Germline data is missing allele description:', field)
            sys.exit(1)

        segment_name = None
        if level == "allele":
            segment_name = field

        if level == "gene":
            segment_name = getGene(allele_description)

        if level == "subgroup":
            segment_name = getSubgroup(allele_description)

        if segment_name is None:
            continue

        if segment_name not in distinct_names:
            distinct_names.append(segment_name)

    return distinct_names
