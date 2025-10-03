"""
TCR clonal assignment
"""

#
# tcr_clone.py
# TCR clonal assignment
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

# repcalc modules
from repcalc import __version__
import repcalc.defaults as defaults
import repcalc.metadata as metadata
import repcalc.gldb as gldb

import json
import math
import numpy
import csv
import airr

# module operations
tcrKey = "TCR"
module_levels = [ "allele", "gene" ]

cdr3_dict = {}
clone_dict = {}
collapse_dict = {}
file_writers = {}

def intersect_genes(g1, g2):
    ga1 = set(g1.split(','))
    ga2 = set(g2.split(','))
    if len(ga1 & ga2) == 0:
        return False
    else:
        return True

def collapse_clones(clone_list):
    if len(clone_list) == 1:
        return clone_list

    collapse_list = []
    remain_list = clone_list
    while len(remain_list) > 0:
        # find the largest as starting point
        max_idx = 0
        max_count = 0
        idx = 0
        for clone in remain_list:
            if clone['duplicate_count'] > max_count:
                max_count = clone['duplicate_count']
                max_idx = idx
            idx += 1
        max_clone = remain_list[max_idx]
        collapse_list.append(max_clone)
        del remain_list[max_idx]
        #print(max_clone)

        new_list = []
        for clone in remain_list:
            if intersect_genes(max_clone['v_call'], clone['v_call']) and intersect_genes(max_clone['j_call'], clone['j_call']):
                max_clone['clone_count'] += clone['clone_count']
                max_clone['duplicate_count'] += clone['duplicate_count']
            else:
                new_list.append(clone)
        remain_list = new_list
    return collapse_list

def assign_clone(inputDict, rep_id, fields, level):
    # first time we have seen this junction
    if cdr3_dict[level][rep_id].get(fields['junction_aa']) is None:
        cdr3_dict[level][rep_id][fields['junction_aa']] = []

    if level == 'gene':
        germline = inputDict[defaults.germline_key]
        v_gene = gldb.transformToLevel(germline, fields['v_call'], level)
        j_gene = gldb.transformToLevel(germline, fields['j_call'], level)

    # compare with existing clones
    for clone in cdr3_dict[level][rep_id][fields['junction_aa']]:
        if level == 'gene':
            if clone['v_call'] == v_gene:
                if clone['j_call'] == j_gene:
                    clone['clone_count'] += 1
                    clone['duplicate_count'] = defaults.get_duplicate_count(fields)
                    return clone
        else:
            if clone['v_call'] == fields['v_call']:
                if clone['j_call'] == fields['j_call']:
                    clone['clone_count'] += 1
                    clone['duplicate_count'] = defaults.get_duplicate_count(fields)
                    return clone

    # new clone if get here
    clone = {}
    if level == 'gene':
        clone['v_call'] = v_gene
        clone['j_call'] = j_gene
    else:
        clone['v_call'] = fields['v_call']
        clone['j_call'] = fields['j_call']
    clone['clone_id'] = clone_dict[level][rep_id]['clone_id']
    clone['clone_count'] = 1
    clone['duplicate_count'] = defaults.get_duplicate_count(fields)
    clone_dict[level][rep_id]['clone_id'] += 1
    cdr3_dict[level][rep_id][fields['junction_aa']].append(clone)
    return clone

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # TCR clonal assignment
    for level in module_levels:
        cdr3_dict[level] = {}
        clone_dict[level] = {}
        file_writers[level] = {}
        for rep_id in metadataDict:
            cdr3_dict[level][rep_id] = {}
            clone_dict[level][rep_id] = { 'clone_id': 1 }

def process_record(inputDict, metadataDict, currentFile, calc, row):
    """Perform calculation from given fields"""
    # TCR clonal assignment
    if tcrKey in calc['operations']:
        # get repertoire
        # TODO: does the repertoire have to be in the metadata file?
        rep_id = row['repertoire_id']
        if metadataDict.get(rep_id) is None:
            return

        # verify TCR locus
        if 'TR' not in row['locus']:
            return

        # productive
        if not row['productive']:
            return

        # non empty junction
        if len(row['junction_aa']) == 0:
            return

        # assign the rearrangement to a clone
        for level in module_levels:
            clone = assign_clone(inputDict, rep_id, row, level)
    
            # write rearrangement with clone_id
            writer = file_writers[level].get(currentFile)
            if writer is None:
                filename = rep_id
                if inputDict.get(defaults.processing_stage_key) is not None:
                    filename = filename + '.' + inputDict.get(defaults.processing_stage_key)
                filename = filename + '.' + level + '.clone.airr.tsv'
                print('Output file: ' + filename)
                writer = airr.derive_rearrangement(filename, currentFile, fields=['clone_id'])
                file_writers[level][currentFile] = writer
            row['clone_id'] = clone['clone_id']
            writer.write(row)

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    # TCR clonal assignment
    if tcrKey in calc['operations']:
        for level in module_levels:
            for rep_id in metadataDict:
                rep = metadataDict[rep_id]
                filename = f"{rep_id}.{level}.clone.tsv"
                writer = csv.DictWriter(open(filename, 'w'), fieldnames=['clone_id', 'junction_aa', 'v_call', 'j_call', 'clone_count', 'duplicate_count'],
                                        dialect='excel-tab',
                                        extrasaction='ignore', lineterminator='\n')
                writer.writeheader()
                filename = f"{rep_id}.{level}.collapsed_clone.tsv"
                col_writer = csv.DictWriter(open(filename, 'w'), fieldnames=['clone_id', 'junction_aa', 'v_call', 'j_call', 'clone_count', 'duplicate_count'],
                                        dialect='excel-tab',
                                        extrasaction='ignore', lineterminator='\n')
                col_writer.writeheader()
                for cdr3 in cdr3_dict[level][rep_id]:
                    clone_list = cdr3_dict[level][rep_id][cdr3]
                    for clone in clone_list:
                        clone['junction_aa'] = cdr3
                        writer.writerow(clone)
                    collapse_list = collapse_clones(clone_list)
                    for clone in collapse_list:
                        clone['junction_aa'] = cdr3
                        col_writer.writerow(clone)









