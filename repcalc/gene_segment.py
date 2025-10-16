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
import sys

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
group_segment_counters = {}
group_segment_counters_productive = {}
group_combo_counters = {}
group_combo_counters_productive = {}

def add_segment_count(germline, segment, level, counters, amount):
    """Add segment counts to counters"""
    if segment is None or len(segment) == 0:
        return
    distinct_names = gldb.getDisplayName(germline, segment, level)
    len_float = float(len(distinct_names))
    if len_float == 0:
        print('WARNING: skipping segment', segment, file=sys.stderr)
        return
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

def compute_frequency(counters):
    """Compute frequencies"""
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

def compute_usage_frequency(counters):
    """Compute frequencies for segment counters"""
    for gene in gene_fields:
        for level in module_levels:
            for mode in calc_modes:
                compute_frequency(counters[gene][level][mode])

def compute_combo_frequency(counters):
    """Compute frequencies for segment combos"""
    for combo in combo_fields:
        for level in counters[combo]:
            for mode in calc_modes:
                compute_frequency(counters[combo][level][mode])

def compute_std(entry, repertoire_counters, field, N, groupDict, gene, level, mode, value):
    """Compute average and standard deviation"""
    group_total = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            rep_value = 0.0
            if repertoire_counters[rep_id][gene].get(level) is not None:
                if repertoire_counters[rep_id][gene][level].get(mode) is not None:
                    if repertoire_counters[rep_id][gene][level][mode].get(value) is not None:
                        rep_value = float(repertoire_counters[rep_id][gene][level][mode][value][field])
            group_total += rep_value
    if N > 0:
        group_avg = group_total / N
    else:
        group_avg = 0.0
    group_std = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            rep_value = 0.0
            if repertoire_counters[rep_id][gene].get(level) is not None:
                if repertoire_counters[rep_id][gene][level].get(mode) is not None:
                    if repertoire_counters[rep_id][gene][level][mode].get(value) is not None:
                        rep_value = float(repertoire_counters[rep_id][gene][level][mode][value][field])
            group_std += (group_avg - rep_value) * (group_avg - rep_value)
    if N > 1:
        group_std = math.sqrt(group_std / (N - 1.0))
    else:
        group_std = 0.0
    entry['N'] = N
    entry[field + '_avg'] = group_avg
    entry[field + '_std'] = group_std

def compute_group_usage(groupDict, counters, repertoire_counters):
    """Calculate group average and deviation for segment counters"""
    for gene in gene_fields:
        for level in module_levels:
            for mode in calc_modes:
                compute_frequency(counters[gene][level][mode])
    N = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            N += 1.0
    for gene in gene_fields:
        for level in module_levels:
            for mode in calc_modes:
                for value in counters[gene][level][mode]:
                    entry = counters[gene][level][mode][value]
                    compute_std(entry, repertoire_counters, 'sequence_count', N, groupDict, gene, level, mode, value)
                    compute_std(entry, repertoire_counters, 'duplicate_count', N, groupDict, gene, level, mode, value)
                    compute_std(entry, repertoire_counters, 'sequence_frequency', N, groupDict, gene, level, mode, value)
                    compute_std(entry, repertoire_counters, 'duplicate_frequency', N, groupDict, gene, level, mode, value)

def compute_group_combo(groupDict, counters, repertoire_counters):
    """Calculate group average and deviation for segment combos"""
    for combo in combo_fields:
        for level in counters[combo]:
            for mode in calc_modes:
                compute_frequency(counters[combo][level][mode])
    N = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            N += 1.0
    for combo in combo_fields:
        for level in counters[combo]:
            for mode in calc_modes:
                for value in counters[combo][level][mode]:
                    entry = counters[combo][level][mode][value]
                    compute_std(entry, repertoire_counters, 'sequence_count', N, groupDict, combo, level, mode, value)
                    compute_std(entry, repertoire_counters, 'duplicate_count', N, groupDict, combo, level, mode, value)
                    compute_std(entry, repertoire_counters, 'sequence_frequency', N, groupDict, combo, level, mode, value)
                    compute_std(entry, repertoire_counters, 'duplicate_frequency', N, groupDict, combo, level, mode, value)

def write_usage_output(id_name, id_value, stage, group_flag, counters, counters_productive):
    """Output segment counts to TSV files"""
    # DEBUG ***
    with open(f'{id_value}.wuo_counters.items.json', 'w') as f:
        json.dump(counters, f, indent=4)
    with open(f'{id_value}.wuo_counters_productive.items.json', 'w') as f:
        json.dump(counters_productive, f, indent=4)

    for gene in gene_fields:
        if group_flag:
            if stage:
                filename = id_value + '.' + stage + '.group.' + gene + '.tsv'
            else:
                filename = id_value + '.group.' + gene + '.tsv'
            writer = open(filename, 'w')
            writer.write(id_name + '\tlevel\tmode\tproductive\tgene\tsequence_count\tduplicate_count\tsequence_frequency\tduplicate_frequency')
            writer.write('\tN\tsequence_count_avg\tsequence_count_std\tsequence_frequency_avg\tsequence_frequency_std')
            writer.write('\tduplicate_count_avg\tduplicate_count_std\tduplicate_frequency_avg\tduplicate_frequency_std\n')
        else:
            if stage:
                filename = id_value + '.' + stage + '.' + gene + '.tsv'
            else:
                filename = id_value + '.' + gene + '.tsv'
            writer = open(filename, 'w')
            writer.write(id_name + '\tlevel\tmode\tproductive\tgene\tsequence_count\tduplicate_count\tsequence_frequency\tduplicate_frequency\n')
        for level in module_levels:
            for mode in calc_modes:
                # output data
                for value in counters[id_value][gene][level][mode]:
                    entry = counters[id_value][gene][level][mode][value]
                    writer.write(id_value + '\t')
                    writer.write(level + '\t')
                    writer.write(mode + '\t')
                    writer.write('FALSE\t')
                    writer.write(value + '\t')
                    writer.write(str(entry['sequence_count']) + '\t')
                    writer.write(str(entry['duplicate_count']) + '\t')
                    writer.write(str(entry['sequence_frequency']) + '\t')
                    if group_flag:
                        writer.write(str(entry['duplicate_frequency']) + '\t')
                        writer.write(str(entry['N']) + '\t')
                        writer.write(str(entry['sequence_count_avg']) + '\t')
                        writer.write(str(entry['sequence_count_std']) + '\t')
                        writer.write(str(entry['sequence_frequency_avg']) + '\t')
                        writer.write(str(entry['sequence_frequency_std']) + '\t')
                        writer.write(str(entry['duplicate_count_avg']) + '\t')
                        writer.write(str(entry['duplicate_count_std']) + '\t')
                        writer.write(str(entry['duplicate_frequency_avg']) + '\t')
                        writer.write(str(entry['duplicate_frequency_std']) + '\n')
                    else:
                        writer.write(str(entry['duplicate_frequency']) + '\n')
                # output data
                for value in counters_productive[id_value][gene][level][mode]:
                    entry = counters_productive[id_value][gene][level][mode][value]
                    writer.write(id_value + '\t')
                    writer.write(level + '\t')
                    writer.write(mode + '\t')
                    writer.write('TRUE\t')
                    writer.write(value + '\t')
                    writer.write(str(entry['sequence_count']) + '\t')
                    writer.write(str(entry['duplicate_count']) + '\t')
                    writer.write(str(entry['sequence_frequency']) + '\t')
                    if group_flag:
                        writer.write(str(entry['duplicate_frequency']) + '\t')
                        writer.write(str(entry['N']) + '\t')
                        writer.write(str(entry['sequence_count_avg']) + '\t')
                        writer.write(str(entry['sequence_count_std']) + '\t')
                        writer.write(str(entry['sequence_frequency_avg']) + '\t')
                        writer.write(str(entry['sequence_frequency_std']) + '\t')
                        writer.write(str(entry['duplicate_count_avg']) + '\t')
                        writer.write(str(entry['duplicate_count_std']) + '\t')
                        writer.write(str(entry['duplicate_frequency_avg']) + '\t')
                        writer.write(str(entry['duplicate_frequency_std']) + '\n')
                    else:
                        writer.write(str(entry['duplicate_frequency']) + '\n')
        writer.close()

def write_combo_output(id_name, id_value, stage, group_flag, counters, counters_productive):
    """Output segment combo counts to TSV files"""
    for combo in combo_fields:
        if group_flag:
            if stage:
                filename = id_value + '.' + stage + '.group.' + combo + '_combo.tsv'
            else:
                filename = id_value + '.group.' + combo + '_combo.tsv'
        else:
            if stage:
                filename = id_value + '.' + stage + '.' + combo + '_combo.tsv'
            else:
                filename = id_value + '.' + combo + '_combo.tsv'
        writer = open(filename, 'w')
        writer.write(id_name + '\tlevel\tmode\tproductive\tcombo\t')
        if 'v' in combo:
            writer.write('v_level\t')
        if 'd' in combo:
            writer.write('d_level\t')
        if 'j' in combo:
            writer.write('j_level\t')
        writer.write('sequence_count\tduplicate_count\tsequence_frequency\tduplicate_frequency')
        if group_flag:
            writer.write('\tN\tsequence_count_avg\tsequence_count_std\tsequence_frequency_avg\tsequence_frequency_std')
            writer.write('\tduplicate_count_avg\tduplicate_count_std\tduplicate_frequency_avg\tduplicate_frequency_std\n')
        else:
            writer.write('\n')
        for level in counters[id_value][combo]:
            for mode in calc_modes:
                # output data
                for value in counters[id_value][combo][level][mode]:
                    entry = counters[id_value][combo][level][mode][value]
                    writer.write(id_value + '\t')
                    writer.write(level + '\t')
                    writer.write(mode + '\t')
                    writer.write('FALSE\t')
                    writer.write(value + '\t')
                    l = value.split('|')
                    idx = 0
                    if 'v' in combo:
                        writer.write(l[idx] + '\t')
                        idx += 1
                    if 'd' in combo:
                        writer.write(l[idx] + '\t')
                        idx += 1
                    if 'j' in combo:
                        writer.write(l[idx] + '\t')
                        idx += 1
                    writer.write(str(entry['sequence_count']) + '\t')
                    writer.write(str(entry['duplicate_count']) + '\t')
                    writer.write(str(entry['sequence_frequency']) + '\t')
                    if group_flag:
                        writer.write(str(entry['duplicate_frequency']) + '\t')
                        writer.write(str(entry['N']) + '\t')
                        writer.write(str(entry['sequence_count_avg']) + '\t')
                        writer.write(str(entry['sequence_count_std']) + '\t')
                        writer.write(str(entry['sequence_frequency_avg']) + '\t')
                        writer.write(str(entry['sequence_frequency_std']) + '\t')
                        writer.write(str(entry['duplicate_count_avg']) + '\t')
                        writer.write(str(entry['duplicate_count_std']) + '\t')
                        writer.write(str(entry['duplicate_frequency_avg']) + '\t')
                        writer.write(str(entry['duplicate_frequency_std']) + '\n')
                    else:
                        writer.write(str(entry['duplicate_frequency']) + '\n')
        # productive
        for level in counters_productive[id_value][combo]:
            for mode in calc_modes:
                for value in counters_productive[id_value][combo][level][mode]:
                    entry = counters_productive[id_value][combo][level][mode][value]
                    writer.write(id_value + '\t')
                    writer.write(level + '\t')
                    writer.write(mode + '\t')
                    writer.write('TRUE\t')
                    writer.write(value + '\t')
                    l = value.split('|')
                    idx = 0
                    if 'v' in combo:
                        writer.write(l[idx] + '\t')
                        idx += 1
                    if 'd' in combo:
                        writer.write(l[idx] + '\t')
                        idx += 1
                    if 'j' in combo:
                        writer.write(l[idx] + '\t')
                        idx += 1
                    writer.write(str(entry['sequence_count']) + '\t')
                    writer.write(str(entry['duplicate_count']) + '\t')
                    writer.write(str(entry['sequence_frequency']) + '\t')
                    if group_flag:
                        writer.write(str(entry['duplicate_frequency']) + '\t')
                        writer.write(str(entry['N']) + '\t')
                        writer.write(str(entry['sequence_count_avg']) + '\t')
                        writer.write(str(entry['sequence_count_std']) + '\t')
                        writer.write(str(entry['sequence_frequency_avg']) + '\t')
                        writer.write(str(entry['sequence_frequency_std']) + '\t')
                        writer.write(str(entry['duplicate_count_avg']) + '\t')
                        writer.write(str(entry['duplicate_count_std']) + '\t')
                        writer.write(str(entry['duplicate_frequency_avg']) + '\t')
                        writer.write(str(entry['duplicate_frequency_std']) + '\n')
                    else:
                        writer.write(str(entry['duplicate_frequency']) + '\n')
        writer.close()

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # gene and combo segment usage
    if inputDict.get(defaults.groups_key) is not None:
        for group in inputDict[defaults.groups_key]:
            combo_counters[group] = {}
            combo_counters_productive[group] = {}
            segment_counters[group] = {}
            segment_counters_productive[group] = {}
            for rep_id in metadataDict:
                segment_counters[group][rep_id] = {}
                segment_counters_productive[group][rep_id] = {}
                for gene in gene_fields:
                    segment_counters[group][rep_id][gene] = {}
                    segment_counters_productive[group][rep_id][gene] = {}
                    for level in module_levels:
                        segment_counters[group][rep_id][gene][level] = {}
                        segment_counters_productive[group][rep_id][gene][level] = {}
                        for mode in calc_modes:
                            segment_counters[group][rep_id][gene][level][mode] = {}
                            segment_counters_productive[group][rep_id][gene][level][mode] = {}
                
                combo_counters[group][rep_id] = {}
                combo_counters_productive[group][rep_id] = {}
                for combo in combo_fields:
                    combo_counters[group][rep_id][combo] = {}
                    combo_counters_productive[group][rep_id][combo] = {}
    
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

    if inputDict.get(defaults.groups_key) is not None:
        for group in inputDict[defaults.groups_key]:
            group_segment_counters[group] = {}
            group_segment_counters_productive[group] = {}
            for gene in gene_fields:
                group_segment_counters[group][gene] = {}
                group_segment_counters_productive[group][gene] = {}
                for level in module_levels:
                    group_segment_counters[group][gene][level] = {}
                    group_segment_counters_productive[group][gene][level] = {}
                    for mode in calc_modes:
                        group_segment_counters[group][gene][level][mode] = {}
                        group_segment_counters_productive[group][gene][level][mode] = {}

            group_combo_counters[group] = {}
            group_combo_counters_productive[group] = {}
            for combo in combo_fields:
                group_combo_counters[group][combo] = {}
                group_combo_counters_productive[group][combo] = {}

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
                    groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
                    if groupList:
                        for group in groupList:
                            # check and apply filter
                            if defaults.apply_filter(inputDict, group, fields):
                                add_segment_count(germline, segment, level, group_segment_counters[group][gene][level], defaults.get_duplicate_count(fields))
                                if fields.get('productive'):
                                    add_segment_count(germline, segment, level, group_segment_counters_productive[group][gene][level], defaults.get_duplicate_count(fields))
                                add_segment_count(germline, segment, level, segment_counters[group][rep_id][gene][level], defaults.get_duplicate_count(fields))
                                if fields.get('productive'):
                                    add_segment_count(germline, segment, level, segment_counters_productive[group][rep_id][gene][level], defaults.get_duplicate_count(fields))

    # combo operations
    if comboKey in calc['operations']:
        for combo in combo_fields:
            add_combo_count(germline, fields, combo, combo_counters[rep_id][combo], defaults.get_duplicate_count(fields))
            if fields.get('productive'):
                add_combo_count(germline, fields, combo, combo_counters_productive[rep_id][combo], defaults.get_duplicate_count(fields))
            groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
            if groupList:
                for group in groupList:
                    # check and apply filter
                    if defaults.apply_filter(inputDict, group, fields):
                        add_combo_count(germline, fields, combo, group_combo_counters[group][combo], defaults.get_duplicate_count(fields))
                        if fields.get('productive'):
                            add_combo_count(germline, fields, combo, group_combo_counters_productive[group][combo], defaults.get_duplicate_count(fields))
                        add_combo_count(germline, fields, combo, combo_counters[group][rep_id][combo], defaults.get_duplicate_count(fields))
                        if fields.get('productive'):
                            add_combo_count(germline, fields, combo, combo_counters_productive[group][rep_id][combo], defaults.get_duplicate_count(fields))


def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""

    # generate usage output
    if usageKey in calc['operations']:
        for rep_id in metadataDict:
            compute_usage_frequency(segment_counters[rep_id])
            compute_usage_frequency(segment_counters_productive[rep_id])
            write_usage_output('repertoire_id', rep_id, inputDict.get(defaults.processing_stage_key), False, segment_counters, segment_counters_productive)
        # repertoire groups
        if inputDict.get(defaults.groups_key) is not None:
            # compute frequencies for group rearrangement filter counts
            for group in inputDict[defaults.groups_key]:
                for rep in inputDict[defaults.groups_key][group]['repertoires']:
                    rep_id = rep['repertoire_id']
                    compute_usage_frequency(segment_counters[group][rep_id])
                    compute_usage_frequency(segment_counters_productive[group][rep_id])

            for group in inputDict[defaults.groups_key]:
                compute_group_usage(inputDict[defaults.groups_key][group], group_segment_counters[group], segment_counters[group])
                compute_group_usage(inputDict[defaults.groups_key][group], group_segment_counters_productive[group], segment_counters_productive[group])
                write_usage_output('repertoire_group_id', group, inputDict.get(defaults.processing_stage_key), True, group_segment_counters, group_segment_counters_productive)
    # DEBUG ***
    with open('group_segment_counters.items.json', 'w') as f:
        json.dump(group_segment_counters, f, indent=4)
    with open('group_segment_counters_productive.items.json', 'w') as f:
        json.dump(group_segment_counters_productive, f, indent=4)

    # generate combo output
    if comboKey in calc['operations']:
        for rep_id in metadataDict:
            compute_combo_frequency(combo_counters[rep_id])
            compute_combo_frequency(combo_counters_productive[rep_id])
            write_combo_output('repertoire_id', rep_id, inputDict.get(defaults.processing_stage_key), False, combo_counters, combo_counters_productive)
        # repertoire groups
        if inputDict.get(defaults.groups_key) is not None:
            for group in inputDict[defaults.groups_key]:
                compute_group_combo(inputDict[defaults.groups_key][group], group_combo_counters[group], combo_counters[group])
                compute_group_combo(inputDict[defaults.groups_key][group], group_combo_counters_productive[group], combo_counters_productive[group])
                write_combo_output('repertoire_group_id', group, inputDict.get(defaults.processing_stage_key), True, group_combo_counters, group_combo_counters_productive)
