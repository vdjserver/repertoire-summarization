"""
Gene segment calculation module
"""

#
# gene_segment.py
# Gene segment usage
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
usageKey = "usage"
comboKey = "combo"
gene_fields = [ "v_call", "d_call", "j_call", "c_call" ]
combo_fields = [ "vj", "vd", "vdj", "dj" ]
module_levels = [ "allele", "gene", "subgroup" ]
calc_modes = [ "unique", "exists", "proportion" ]

segment_counters = {}
segment_counters_productive = {}
combo_counters = {}
combo_counters_productive = {}
summary_combo_counters = {}

def add_segment_count(germline, segment, level, counters, amount):
    """Add segment counts to counters"""
    if segment is None or len(segment) == 0:
        return
    distinct_names = gldb.getDisplayName(germline, segment, level)
    len_float = float(len(distinct_names))
    #print(segment, distinct_names, len(distinct_names))

    # unique mode
    if len(distinct_names) == 1:
        segment_name = distinct_names[0]
        if counters['unique'].get(segment_name) is None:
            counters['unique'][segment_name] = { "sequence_count":0, "sequence_frequency":0.0, "duplicate_count":0, "duplicate_frequency":0.0 }
        counters['unique'][segment_name]['sequence_count'] += 1
        counters['unique'][segment_name]['duplicate_count'] += amount

    # exists mode
    for segment_name in distinct_names:
        if counters['exists'].get(segment_name) is None:
            counters['exists'][segment_name] = { "sequence_count":0, "sequence_frequency":0.0, "duplicate_count":0, "duplicate_frequency":0.0 }
        counters['exists'][segment_name]['sequence_count'] += 1
        counters['exists'][segment_name]['duplicate_count'] += amount

    # proportion
    for segment_name in distinct_names:
        if counters['proportion'].get(segment_name) is None:
            counters['proportion'][segment_name] = { "sequence_count":0.0, "sequence_frequency":0.0, "duplicate_count":0.0, "duplicate_frequency":0.0 }
        counters['proportion'][segment_name]['sequence_count'] += 1.0 / len_float
        counters['proportion'][segment_name]['duplicate_count'] += float(amount) / len_float
    
    #print(counters)
    return

def compute_frequency(counters):
    """Compute frequencies for segment counters"""
    total = 0.0
    dup_total = 0.0
    for value in counters:
        entry = counters[value]
        total = total + float(entry['sequence_count'])
        dup_total = dup_total + float(entry['duplicate_count'])
    for value in counters:
        entry = counters[value]
        entry['sequence_frequency'] = float(entry['sequence_count']) / total
        entry['duplicate_frequency'] = float(entry['duplicate_count']) / dup_total
    return

def add_combo_count(germline, fields, combo, counters, amount):
    """Add segment combo counts to counters"""
    combo_list = []
    if 'v' in combo:
        vcall = fields.get('v_call')
        if not vcall:
            return
        combo_list.append(vcall)
    if 'd' in combo:
        dcall = fields.get('d_call')
        if not dcall:
            # handle empty D call as special
            dcall = 'NONE'
        combo_list.append(dcall)
    if 'j' in combo:
        jcall = fields.get('j_call')
        if not jcall:
            return
        combo_list.append(jcall)

    for level in module_levels:
        level_names = []
        combo_names = []
        for segment in combo_list:
            if segment == 'NONE':
                # special handling of empty D call
                combo_names.append(['NONE'])
            else:
                combo_names.append(gldb.getDisplayName(germline, segment, level))
            level_names.append(level)
        level_name = '|'.join(level_names)
        #print(level_name)

        if counters.get(level_name) is None:
            counters[level_name] = {}
            for mode in calc_modes:
                counters[level_name][mode] = {}

        is_unique = True
        for distinct_names in combo_names:
            if len(distinct_names) != 1:
                is_unique = False
                break

        # unique mode
        if is_unique:
            segment_names = []
            for distinct_names in combo_names:
                segment_names.append(distinct_names[0])
            segment_name = '|'.join(segment_names)
            #print(segment_name)
            if counters[level_name]['unique'].get(segment_name) is None:
                counters[level_name]['unique'][segment_name] = { "sequence_count":0, "sequence_frequency":0.0, "duplicate_count":0, "duplicate_frequency":0.0 }
            counters[level_name]['unique'][segment_name]['sequence_count'] += 1
            counters[level_name]['unique'][segment_name]['duplicate_count'] += amount

        # need to generate all combos, easier with for loops
        if len(combo_names) == 2:
            len_float = float(len(combo_names[0]) * len(combo_names[1]))
            for distinct_name1 in combo_names[0]:
                for distinct_name2 in combo_names[1]:
                    segment_name = distinct_name1 + '|' + distinct_name2
                    #print(segment_name)
                    # exists mode
                    if counters[level_name]['exists'].get(segment_name) is None:
                        counters[level_name]['exists'][segment_name] = { "sequence_count":0, "sequence_frequency":0.0, "duplicate_count":0, "duplicate_frequency":0.0 }
                    counters[level_name]['exists'][segment_name]['sequence_count'] += 1
                    counters[level_name]['exists'][segment_name]['duplicate_count'] += amount
                    # proportion
                    if counters[level_name]['proportion'].get(segment_name) is None:
                        counters[level_name]['proportion'][segment_name] = { "sequence_count":0.0, "sequence_frequency":0.0, "duplicate_count":0.0, "duplicate_frequency":0.0 }
                    counters[level_name]['proportion'][segment_name]['sequence_count'] += 1.0 / len_float
                    counters[level_name]['proportion'][segment_name]['duplicate_count'] += float(amount) / len_float
        else:
            len_float = float(len(combo_names[0]) * len(combo_names[1]) * len(combo_names[2]))
            for distinct_name1 in combo_names[0]:
                for distinct_name2 in combo_names[1]:
                    for distinct_name3 in combo_names[2]:
                        segment_name = distinct_name1 + '|' + distinct_name2 + '|' + distinct_name3
                        #print(segment_name)
                        # exists mode
                        if counters[level_name]['exists'].get(segment_name) is None:
                            counters[level_name]['exists'][segment_name] = { "sequence_count":0, "sequence_frequency":0.0, "duplicate_count":0, "duplicate_frequency":0.0 }
                        counters[level_name]['exists'][segment_name]['sequence_count'] += 1
                        counters[level_name]['exists'][segment_name]['duplicate_count'] += amount
                        # proportion
                        if counters[level_name]['proportion'].get(segment_name) is None:
                            counters[level_name]['proportion'][segment_name] = { "sequence_count":0.0, "sequence_frequency":0.0, "duplicate_count":0.0, "duplicate_frequency":0.0 }
                        counters[level_name]['proportion'][segment_name]['sequence_count'] += 1.0 / len_float
                        counters[level_name]['proportion'][segment_name]['duplicate_count'] += float(amount) / len_float
    return

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # gene and combo segment usage
    for rep_id in metadataDict:
        segment_counters[rep_id] = {}
        segment_counters_productive[rep_id] = {}
        for gene in gene_fields:
            segment_counters[rep_id][gene] = {}
            segment_counters_productive[rep_id][gene] = {}
            for level in module_levels:
                segment_counters[rep_id][gene][level] = {}
                segment_counters_productive[rep_id][gene][level] = {}
                for mode in calc_modes:
                    segment_counters[rep_id][gene][level][mode] = {}
                    segment_counters_productive[rep_id][gene][level][mode] = {}

        combo_counters[rep_id] = {}
        combo_counters_productive[rep_id] = {}
        for combo in combo_fields:
            combo_counters[rep_id][combo] = {}
            combo_counters_productive[rep_id][combo] = {}
    return

def process_record(inputDict, metadataDict, currentFile, calc, fields):
    """Perform calculation from given fields"""

    rep_id = fields.get('repertoire_id')
    if rep_id is None:
        return
    if metadataDict.get(rep_id) is None:
        return
    germline = inputDict.get(defaults.germline_key)
    if germline is None:
        return

    # usage operations
    if usageKey in calc['operations']:
        # TODO: repertoire groups
        for gene in gene_fields:
            segment = fields.get(gene)
            if segment is not None:
                for level in module_levels:
                    add_segment_count(germline, segment, level, segment_counters[rep_id][gene][level], defaults.get_duplicate_count(fields))
                    if fields.get('productive'):
                        add_segment_count(germline, segment, level, segment_counters_productive[rep_id][gene][level], defaults.get_duplicate_count(fields))

    # combo operations
    if comboKey in calc['operations']:
        for combo in combo_fields:
            add_combo_count(germline, fields, combo, combo_counters[rep_id][combo], defaults.get_duplicate_count(fields))
            if fields.get('productive'):
                add_combo_count(germline, fields, combo, combo_counters_productive[rep_id][combo], defaults.get_duplicate_count(fields))


def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""

    # generate usage output
    if usageKey in calc['operations']:
        for rep_id in metadataDict:
            rep = metadataDict[rep_id]
            for gene in gene_fields:
                filename = rep_id + '.' + gene + '.tsv'
                writer = open(filename, 'w')
                writer.write('repertoire_id\tlevel\tmode\tproductive\tgene\tsequence_count\tduplicate_count\tsequence_frequency\tduplicate_frequency\n')
                for level in module_levels:
                    for mode in calc_modes:
                        # compute frequencies
                        compute_frequency(segment_counters[rep_id][gene][level][mode])
                        # output data
                        for value in segment_counters[rep_id][gene][level][mode]:
                            entry = segment_counters[rep_id][gene][level][mode][value]
                            writer.write(rep_id + '\t')
                            writer.write(level + '\t')
                            writer.write(mode + '\t')
                            writer.write('\t')
                            writer.write(value + '\t')
                            writer.write(str(entry['sequence_count']) + '\t')
                            writer.write(str(entry['duplicate_count']) + '\t')
                            writer.write(str(entry['sequence_frequency']) + '\t')
                            writer.write(str(entry['duplicate_frequency']) + '\n')
                        # productive
                        compute_frequency(segment_counters_productive[rep_id][gene][level][mode])
                        # output data
                        for value in segment_counters_productive[rep_id][gene][level][mode]:
                            entry = segment_counters_productive[rep_id][gene][level][mode][value]
                            writer.write(rep_id + '\t')
                            writer.write(level + '\t')
                            writer.write(mode + '\t')
                            writer.write('TRUE\t')
                            writer.write(value + '\t')
                            writer.write(str(entry['sequence_count']) + '\t')
                            writer.write(str(entry['duplicate_count']) + '\t')
                            writer.write(str(entry['sequence_frequency']) + '\t')
                            writer.write(str(entry['duplicate_frequency']) + '\n')
                writer.close()

    # generate combo output
    if comboKey in calc['operations']:
        for rep_id in metadataDict:
            for combo in combo_fields:
                filename = rep_id + '.' + combo + '_combo.tsv'
                writer = open(filename, 'w')
                writer.write('repertoire_id\tlevel\tmode\tproductive\tcombo\tsequence_count\tduplicate_count\tsequence_frequency\tduplicate_frequency\n')
                for level in combo_counters[rep_id][combo]:
                    for mode in calc_modes:
                        # compute frequencies
                        compute_frequency(combo_counters[rep_id][combo][level][mode])
                        # output data
                        for value in combo_counters[rep_id][combo][level][mode]:
                            entry = combo_counters[rep_id][combo][level][mode][value]
                            writer.write(rep_id + '\t')
                            writer.write(level + '\t')
                            writer.write(mode + '\t')
                            writer.write('\t')
                            writer.write(value + '\t')
                            writer.write(str(entry['sequence_count']) + '\t')
                            writer.write(str(entry['duplicate_count']) + '\t')
                            writer.write(str(entry['sequence_frequency']) + '\t')
                            writer.write(str(entry['duplicate_frequency']) + '\n')
                # productive
                for level in combo_counters_productive[rep_id][combo]:
                    for mode in calc_modes:
                        compute_frequency(combo_counters_productive[rep_id][combo][level][mode])
                        # output data
                        for value in combo_counters_productive[rep_id][combo][level][mode]:
                            entry = combo_counters_productive[rep_id][combo][level][mode][value]
                            writer.write(rep_id + '\t')
                            writer.write(level + '\t')
                            writer.write(mode + '\t')
                            writer.write('TRUE\t')
                            writer.write(value + '\t')
                            writer.write(str(entry['sequence_count']) + '\t')
                            writer.write(str(entry['duplicate_count']) + '\t')
                            writer.write(str(entry['sequence_frequency']) + '\t')
                            writer.write(str(entry['duplicate_frequency']) + '\n')
                writer.close()
