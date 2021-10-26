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
    #print(len(fields), distinct_names, len(distinct_names))

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

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # gene segment usage
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
    return

def process_record(inputDict, metadataDict, currentFile, calc, fields):
    """Perform calculation from given fields"""

    # usage operations
    if usageKey in calc['operations']:
        rep_id = fields.get('repertoire_id')
        germline = inputDict.get(defaults.germline_key)
        if rep_id is None:
            # TODO: error
            return
        if metadataDict.get(rep_id) is None:
            # TODO: error
            return
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
        invert_hierarchy = gldb.getInvertHierarchyBy(inputDict[defaults.organismKey])

        def update_combo(level, combo):
            for group in groupSet:
                if combo_counters[group].get(level) is None: combo_counters[group][level] = {}
                if combo_counters[group][level].get(combo) is None:
                    combo_counters[group][level][combo] = { 'count': 0, 'total_count': 0 }
                combo_entry = combo_counters[group][level][combo]
                combo_entry['count'] = combo_entry['count'] + 1
                combo_entry['total_count'] = combo_entry['total_count'] + defaults.get_dupcount(headerMapping, fields)

        for level in calc['levels']:
            if level == 'vj':
                vcall = fields[headerMapping[defaults.headerNames['V_CALL']]]
                if not vcall or vcall == 'None': continue
                vgene = invert_hierarchy.get(vcall)
                if not vgene: continue
                vgene = vgene.keys()[0]
                if not vgene: continue
                jcall = fields[headerMapping[defaults.headerNames['J_CALL']]]
                if not jcall or jcall == 'None': continue
                jgene = invert_hierarchy.get(jcall)
                if not jgene: continue
                jgene = jgene.keys()[0]
                if not jgene: continue
                combo = vgene + '|' + jgene
                update_combo(level, combo)
            if level == 'vdj':
                vcall = fields[headerMapping[defaults.headerNames['V_CALL']]]
                if not vcall or vcall == 'None': continue
                vgene = invert_hierarchy.get(vcall)
                if not vgene: continue
                vgene = vgene.keys()[0]
                if not vgene: continue
                dcall = fields[headerMapping[defaults.headerNames['D_CALL']]]
                if not dcall or dcall == 'None': continue
                dgene = invert_hierarchy.get(dcall)
                if not dgene: continue
                dgene = dgene.keys()[0]
                if not dgene: continue
                jcall = fields[headerMapping[defaults.headerNames['J_CALL']]]
                if not jcall or jcall == 'None': continue
                jgene = invert_hierarchy.get(jcall)
                if not jgene: continue
                jgene = jgene.keys()[0]
                if not jgene: continue
                combo = vgene + '|' + dgene + '|' + jgene
                update_combo(level, combo)
                    

def make_tree(hierarchy, segment_counters, depth):
    siblings = []
    total_count = 0
    for label, sub_hierarchy in hierarchy.items():
        #print (label)
        if len(sub_hierarchy) > 0:
            children = make_tree(sub_hierarchy, segment_counters, depth + 1)
            count = 0
            for child in children:
                count += child['absolute']
            total_count += count
            siblings.append({
                'label': label,
                'children': children,
                'absolute': count,
            })
        else:
            if label in segment_counters:
                count = segment_counters[label]
                total_count += count
                siblings.append({
                    'label': label,
                    'absolute': count,
                })
            else:
                siblings.append({
                    'label': label,
                    'absolute': 0,
                })
    for sibling in siblings:
        if total_count > 0:
            # don't sum across gene families
            if depth == 1 and sibling['absolute'] > 0: sibling['relative'] = 1.0
            else: sibling['relative'] = float(sibling['absolute'])/float(total_count)
        else:
            sibling['relative'] = 0.0
    return siblings

def prune_tree(tree, key):
    del tree[key]
    if 'children' in tree:
        children = []
        for child in tree['children']:
            children.append(prune_tree(child, key))
        tree['children'] = children
    return tree

def group_tree(inputDict, sampleGroup, samples, hierarchy, sampleTrees):
    N = float(len(samples))

    siblings = []
    for label, sub_hierarchy in hierarchy.items():
        absolute = 0.0
        absolute_std = 0.0
        relative = 0.0
        relative_std = 0.0
        sampleSubTrees = {}
        for sampleName in samples:
            tree = sampleTrees[sampleName]
            sampleNode = None
            for child in tree:
                if child['label'] == label:
                    sampleNode = child
                    break
            if sampleNode:
                absolute += sampleNode['absolute']
                absolute_std += sampleNode['absolute'] * sampleNode['absolute']
                relative += sampleNode['relative']
                relative_std += sampleNode['relative'] * sampleNode['relative']
                if sampleNode.get('children'):
                    sampleSubTrees[sampleName] = sampleNode['children']
        if N > 1:
            absolute_std = math.sqrt((absolute_std - absolute * absolute / N) / (N - 1.0))
            relative_std = math.sqrt((relative_std - relative * relative / N) / (N - 1.0))
        else:
            absolute_std = 0.0
            relative_std = 0.0
        absolute /= N
        relative /= N

        children = None
        nodeEntry = { 'label': label,
                      'absolute': absolute,
                      'absolute_std': absolute_std,
                      'relative': relative,
                      'relative_std': relative_std }
        if len(sub_hierarchy) > 0:
            # internal node
            children = group_tree(inputDict, sampleGroup, samples, sub_hierarchy, sampleSubTrees)
            nodeEntry['children'] = children
        siblings.append(nodeEntry)

    return siblings    

def generate_combo_summary(inputDict, outputSpec, calc):
    groups = inputDict[defaults.groupsKey]
    for group in groups:
        if (groups[group]['type'] == 'sampleGroup'):
            for level in calc['levels']:
                if combo_counters[group].get(level) is None: continue
                if summary_combo_counters.get(group) is None: summary_combo_counters[group] = {}
                if summary_combo_counters[group].get(level) is None: summary_combo_counters[group][level] = {}
                for entry in combo_counters[group][level]:
                    if summary_combo_counters[group][level].get(entry) is None: summary_combo_counters[group][level][entry] = { }
                    combo_entry = summary_combo_counters[group][level][entry]
                    samples = groups[group]['samples']
                    count = 0
                    sampleList = []
                    for sample in samples:
                        sampleName = metadata.sample_for_uuid(inputDict, sample)
                        if combo_counters.get(sampleName) is None: continue
                        if combo_counters[sampleName].get(level) is None: continue
                        if combo_counters[sampleName][level].get(entry) is not None:
                            count += 1
                            sampleList.append(sampleName)
                    combo_entry['count'] = combo_counters[group][level][entry]['count']
                    combo_entry['total_count'] = combo_counters[group][level][entry]['total_count']
                    combo_entry['level'] = count
                    combo_entry['samples'] = sampleList

def generate_combo_detail(inputDict, outputSpec, calc):
    groups = inputDict[defaults.groupsKey]
    for level in calc['levels']:
        for group in groups:
            if combo_counters.get(group) is None: continue
            if combo_counters[group].get(level) is None: continue
            combo_entry = combo_counters[group][level]
            #print(combo_entry)
            filename = group + "_" + level + "_segment_combos.tsv"
            writer = open(filename, 'w')
            if (groups[group]['type'] == 'sampleGroup'):
                summary_combo_entry = summary_combo_counters[group][level]
                writer.write('COMBO\tSHARING_LEVEL\tCOUNT\tTOTAL_COUNT\tSAMPLES\tSAMPLE_COUNT\tSAMPLE_TOTAL_COUNT\n')
                samples = groups[group]['samples']
                for sl in range(1, len(samples)+1, 1):
                    string_array = []
                    count_array = []
                    for entry in summary_combo_entry:
                        sl_entry = summary_combo_entry[entry]
                        if sl_entry['level'] != sl: continue
                        if sl_entry['total_count'] == 0: continue
                        line = entry + '\t' + str(sl)
                        line += '\t' + str(sl_entry['count']) + '\t' + str(sl_entry['total_count'])
                        line += '\t' + ','.join(sl_entry['samples'])
                        for sampleName in sl_entry['samples']:
                            line += '\t' + str(combo_counters[sampleName][level][entry]['count'])
                            line += '\t' + str(combo_counters[sampleName][level][entry]['total_count'])
                        line += '\n'
                        string_array.append(line)
                        count_array.append(sl_entry['total_count'])
                    sort_array = numpy.array(count_array).argsort()
                    for idx in reversed(sort_array):
                        writer.write(string_array[idx])
            else:
                writer.write('COMBO\tCOUNT\tTOTAL_COUNT\n')
                string_array = []
                count_array = []
                for entry in combo_entry:
                    line = entry + '\t' + str(combo_entry[entry]['count']) + '\t' + str(combo_entry[entry]['total_count']) + '\n'
                    string_array.append(line)
                    count_array.append(combo_entry[entry]['total_count'])
                sort_array = numpy.array(count_array).argsort()
                for idx in reversed(sort_array):
                    writer.write(string_array[idx])
            writer.close()
            # output specification for process metadata
            if (not outputSpec['files'].get(group + "_gene_segment_combos")): outputSpec['files'][group + "_gene_segment_combos"] = {}
            outputSpec['groups'][group]['gene_segment_combos'] = { "files": group + "_gene_segment_combos", "type": "output" }
            outputSpec['files'][group + "_gene_segment_combos"][level] = { "value": filename, "description":"Gene Segment Combos", "type":"tsv" }

def generate_combo_comparison(inputDict, outputSpec, calc):
    groups = inputDict[defaults.groupsKey]

    # separate sample groups from singles
    sampleGroups = []
    singleGroups = []
    for group in groups:
        if (groups[group]['type'] == 'sampleGroup'): sampleGroups.append(group)
        else: singleGroups.append(group)

    # single groups
    for level in calc['levels']:
        filename1 = "group_shared_" + level + "_segment_combos.tsv"
        filename2 = "group_diff_" + level + "_segment_combos.tsv"
        writer1 = open(filename1, 'w')
        writer2 = open(filename2, 'w')
        for rowGroup in singleGroups:
            if combo_counters.get(rowGroup) is None: continue
            row_entry = combo_counters[rowGroup].get(level)
            if not row_entry: continue
            A = set(row_entry.keys())
            string_array = []
            count_array = []
            diff_string_array = []
            diff_count_array = []
            for colGroup in singleGroups:
                if metadata.sample_contains_file(inputDict, rowGroup, colGroup): continue
                if combo_counters.get(colGroup) is None: continue
                col_entry = combo_counters[colGroup].get(level)
                if not col_entry: continue
                B = set(col_entry.keys())
                singleShared = A & B
                singleDiff = A - B
                if len(singleShared) > 0:
                    for entry in singleShared:
                        line = entry
                        line += '\t' + rowGroup + '\t' + str(row_entry[entry]['count']) + '\t' + str(row_entry[entry]['total_count'])
                        line += '\t' + colGroup + '\t' + str(col_entry[entry]['count']) + '\t' + str(col_entry[entry]['total_count'])
                        line += '\n'
                        string_array.append(line)
                        count_array.append(row_entry[entry]['total_count'])
                if len(singleDiff) > 0:
                    for entry in singleDiff:
                        line = entry
                        line += '\t' + rowGroup + '\t' + str(row_entry[entry]['count']) + '\t' + str(row_entry[entry]['total_count'])
                        line += '\t' + colGroup
                        line += '\n'
                        diff_string_array.append(line)
                        diff_count_array.append(row_entry[entry]['total_count'])
            writer1.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
            writer1.write('COMBO\tGROUP_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\tCOUNT_B\tTOTAL_COUNT_B\n')
            sort_array = numpy.array(count_array).argsort()
            for idx in reversed(sort_array):
                writer1.write(string_array[idx])
            writer2.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
            writer2.write('COMBO\tGROUP_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\n')
            sort_array = numpy.array(diff_count_array).argsort()
            for idx in reversed(sort_array):
                writer2.write(diff_string_array[idx])
        writer1.close()
        writer2.close()
        # output specification for process metadata
        # TODO: Cheat with app name of RepCalc, we should use the app name in process metadata
        if (not outputSpec['groups'].get("RepCalc")): outputSpec['groups']["RepCalc"] = {}
        outputSpec['groups']['RepCalc']['gene_segment_combos'] = { "files": "RepCalc_gene_segment_combos", "type": "output" }
        if (not outputSpec['files'].get("RepCalc_gene_segment_combos")): outputSpec['files']["RepCalc_gene_segment_combos"] = {}
        outputSpec['files']["RepCalc_gene_segment_combos"]['group_shared_' + level] = { "value": filename1, "description":"Gene Segment Combos", "type":"tsv" }
        outputSpec['files']["RepCalc_gene_segment_combos"]['group_diff_' + level] = { "value": filename2, "description":"Gene Segment Combos", "type":"tsv" }

    # sample groups
    for level in calc['levels']:
        filename1 = "sampleGroup_shared_" + level + "_segment_combos.tsv"
        filename2 = "sampleGroup_diff_" + level + "_segment_combos.tsv"
        writer1 = open(filename1, 'w')
        writer2 = open(filename2, 'w')
        for rowGroup in sampleGroups:
            row_entry = summary_combo_counters[rowGroup][level]
            A = set(row_entry.keys())
            string_array = []
            count_array = []
            diff_string_array = []
            diff_count_array = []
            for colGroup in sampleGroups:
                if rowGroup == colGroup: continue
                col_entry = summary_combo_counters[colGroup][level]
                B = set(col_entry.keys())
                groupShared = A & B
                groupDiff = A - B
                if len(groupShared) > 0:
                    for entry in groupShared:
                        line = entry
                        line += '\t' + rowGroup + '\t' + str(row_entry[entry]['level'])
                        line += '\t' + str(row_entry[entry]['count']) + '\t' + str(row_entry[entry]['total_count'])
                        line += '\t' + colGroup + '\t' + str(col_entry[entry]['level'])
                        line += '\t' + str(col_entry[entry]['count']) + '\t' + str(col_entry[entry]['total_count'])
                        line += '\n'
                        string_array.append(line)
                        count_array.append(row_entry[entry]['total_count'])
                if len(groupDiff) > 0:
                    for entry in groupDiff:
                        line = entry
                        line += '\t' + rowGroup + '\t' + str(row_entry[entry]['level'])
                        line += '\t' + str(row_entry[entry]['count']) + '\t' + str(row_entry[entry]['total_count'])
                        line += '\t' + colGroup
                        line += '\n'
                        diff_string_array.append(line)
                        diff_count_array.append(row_entry[entry]['total_count'])
                writer1.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer1.write('COMBO\tGROUP_A\tSHARE_LEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\tSHARE_LEVEL_B\tCOUNT_B\tTOTAL_COUNT_B\n')
                sort_array = numpy.array(count_array).argsort()
                for idx in reversed(sort_array):
                    writer1.write(string_array[idx])
                writer2.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer2.write('COMBO\tGROUP_A\tSHARE_LEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\n')
                sort_array = numpy.array(diff_count_array).argsort()
                for idx in reversed(sort_array):
                    writer2.write(diff_string_array[idx])
        writer1.close()
        writer2.close()
        # output specification for process metadata
        # TODO: Cheat with app name of RepCalc, we should use the app name in process metadata
        if (not outputSpec['groups'].get("RepCalc")): outputSpec['groups']["RepCalc"] = {}
        outputSpec['groups']['RepCalc']['gene_segment_combos'] = { "files": "RepCalc_gene_segment_combos", "type": "output" }
        if (not outputSpec['files'].get("RepCalc_gene_segment_combos")): outputSpec['files']["RepCalc_gene_segment_combos"] = {}
        outputSpec['files']["RepCalc_gene_segment_combos"]['sampleGroup_shared_' + level] = { "value": filename1, "description":"Gene Segment Combos", "type":"tsv" }
        outputSpec['files']["RepCalc_gene_segment_combos"]['sampleGroup_diff_' + level] = { "value": filename2, "description":"Gene Segment Combos", "type":"tsv" }



def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""

    # usage operations
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
                            writer.write(mode + '\t')
                            writer.write(level + '\t')
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
                            writer.write(mode + '\t')
                            writer.write(level + '\t')
                            writer.write('TRUE\t')
                            writer.write(value + '\t')
                            writer.write(str(entry['sequence_count']) + '\t')
                            writer.write(str(entry['duplicate_count']) + '\t')
                            writer.write(str(entry['sequence_frequency']) + '\t')
                            writer.write(str(entry['duplicate_frequency']) + '\n')
                writer.close()

    if comboKey in calc['operations']:
        generate_combo_summary(inputDict, outputSpec, calc)
        generate_combo_detail(inputDict, outputSpec, calc)
        generate_combo_comparison(inputDict, outputSpec, calc)
                    
