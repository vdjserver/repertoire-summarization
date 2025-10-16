"""
CDR3 calculation module
"""

#
# cdr3.py
# CDR3 junction analysis
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
import sys

# module operations
lengthKey = "length"
distributionKey = "distribution"
sharedKey = "shared"
compareKey = "compare"

#
# junction length functions
#
length_levels = [ 'nucleotide', 'nucleotide_wo_duplicates', 'aa', 'aa_wo_duplicates' ]
rel_length_levels = [ 'rel_' + obj for obj in length_levels ]

# repertoire counts
# first key is repertoire_id
# second key is level like aa, nucleotide, which holds the count array
cdr3_histograms = {}
cdr3_histograms_productive = {}

# group counts
# first key is repertoire_group_id
# second key is repertoire_id (to handle rearrangement filters)
# third key is level like aa, nucleotide, which holds the count array
group_cdr3_histograms = {}
group_cdr3_histograms_productive = {}

# group totals (totals across all repertoires for group)
# first key is repertoire_group_id
# second key is level like aa, nucleotide, which holds the count array
group_cdr3_total_histograms = {}
group_cdr3_total_histograms_productive = {}

def increment_count(anArray, length):
    if length == 0: return
    if length >= len(anArray):
        for i in range(len(anArray), length+1, 1): anArray.append(0.0)
    anArray[length] += 1

def increment_duplicate_count(anArray, length, duplicate_count):
    if length == 0: return
    if length >= len(anArray):
        for i in range(len(anArray), length+1, 1): anArray.append(0.0)
    anArray[length] += duplicate_count;

def group_counts(groupArray, repArray):
    if len(repArray) > len(groupArray):
        for i in range(len(groupArray), len(repArray), 1):
            groupArray.append(0.0)
    for i in range(0, len(repArray)):
        groupArray[i] += repArray[i]

def group_average(inputDict, sampleGroup, level):
    groups = inputDict[defaults.groupsKey]
    samples = groups[sampleGroup]['samples']

    cdr3_histograms[sampleGroup][level + '_std'] = []

    # average/std
    N = float(len(samples))
    for sample in samples:
        #sampleName = metadata.sample_for_uuid(inputDict, sample)
        #print(sampleGroup, sample, sampleName)

        sampleArray = cdr3_histograms[sample][level]
        groupArray = cdr3_histograms[sampleGroup][level]
        stdArray = cdr3_histograms[sampleGroup][level + '_std']
        if len(sampleArray) >= len(groupArray):
            for i in range(len(groupArray), len(sampleArray)+1, 1): groupArray.append(0.0)
            for i in range(len(stdArray), len(sampleArray)+1, 1): stdArray.append(0.0)
        for i in range(0, len(sampleArray), 1): groupArray[i] += sampleArray[i]
        for i in range(0, len(sampleArray), 1): stdArray[i] += sampleArray[i] * sampleArray[i]

    groupArray = cdr3_histograms[sampleGroup][level]
    stdArray = cdr3_histograms[sampleGroup][level + '_std']
    for i in range(0, len(stdArray), 1): 
        if N > 1: stdArray[i] = math.sqrt((stdArray[i] - groupArray[i] * groupArray[i] / N) / (N - 1.0))
        else: stdArray[i] = 0.0
    for i in range(0, len(groupArray), 1): groupArray[i] /= N

def compute_relative(data_dict, group, level):
    lenArray = data_dict[group][level]
    totalCount = 0.0
    for i in range(0, len(lenArray), 1): totalCount += lenArray[i]
    if totalCount == 0: totalCount = 1.0
    relArray = []
    for i in range(0, len(lenArray), 1): relArray.append(lenArray[i] / totalCount)
    data_dict[group]['rel_' + level] = relArray

def output_counts(filename, cntArray, relArray):
    #print(filename)
    #print(cntArray)
    #print(relArray)
    writer = open(filename, 'w')
    writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_RELATIVE\n')
    for i in range(0, len(cntArray), 1):
        writer.write(str(i) + '\t' + str(cntArray[i]) + '\t' + str(relArray[i]) + '\n')
    writer.close()

def output_counts_table(filename, groupCounts, reps, repCounts, level):
    names = [ 'repertoire_id', 'total', 'avg' ]
    for i in range(0, len(groupCounts), 1):
        names.append(str(i))
    writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
    writer.writeheader()
    for rep in reps:
        total = 0
        avg = 0
        entry = { 'repertoire_id': rep['repertoire_id'] }
        cntArray = repCounts[entry['repertoire_id']][level]
        for i in range(0, len(cntArray), 1):
            entry[str(i)] = cntArray[i]
            avg += i * cntArray[i]
            total += cntArray[i]
        entry['total'] = total
        if total == 0:
            avg = 0
        else:
            avg = avg / total
        entry['avg'] = avg
        writer.writerow(entry);

#
# AA/NT distribution operations
#
aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
nt_list = ['A', 'C', 'G', 'T']

# first key is group
# second key is aa/nt
# third key is length
# fourth key is position
dist_counters = {}
dist_counters_productive = {}

group_dist_counters = {}
group_dist_counters_productive = {}

group_dist_counters_total = {}
group_dist_counters_productive_total = {}

def aa_distribution(counter, cdr3):
    l = len(cdr3)
    if l < 1: return

    # skip CDR3 if any invalid residues
    for aa in cdr3:
        if aa not in aa_list:
            return

    # initialize counters if necessary
    cnt_dict = counter.get(l)
    if cnt_dict is None:
        counter[l] = {}
        cnt_dict = counter[l]
        for pos in range(l):
            cnt_dict[pos] = {}
            for aa in aa_list:
                cnt_dict[pos][aa] = 0

    for pos in range(l):
        aa = cdr3[pos]
        cnt_dict[pos][aa] = cnt_dict[pos][aa] + 1

def nt_distribution(counter, cdr3):
    l = len(cdr3)
    if l < 1: return

    # skip CDR3 if any invalid nucleotides
    for nt in cdr3:
        if nt not in nt_list:
            return

    # initialize counters if necessary
    cnt_dict = counter.get(l)
    if cnt_dict is None:
        counter[l] = {}
        cnt_dict = counter[l]
        for pos in range(l):
            cnt_dict[pos] = {}
            for nt in nt_list:
                cnt_dict[pos][nt] = 0

    for pos in range(l):
        nt = cdr3[pos]
        cnt_dict[pos][nt] = cnt_dict[pos][nt] + 1

def output_distribution(filename, cntArray, symbol_list):
    writer = open(filename, 'w')
    writer.write('CDR3_LENGTH\tCDR3_POSITION\t' + '\t'.join(symbol_list) + '\n')
    for l in cntArray:
        for pos in cntArray[l]:
            writer.write(str(l) + '\t' + str(pos))
            for symbol in symbol_list:
                writer.write('\t' + str(cntArray[l][pos][symbol]))
            writer.write('\n')
    writer.close()


def compute_relative_distribution(group, level, symbol_list):
    cntArray = dist_counters[group][level]
    relArray = {}
    for l in cntArray:
        if relArray.get(l) is None: relArray[l] = {}
        for pos in cntArray[l]:
            if relArray[l].get(pos) is None: relArray[l][pos] = {}
            totalCount = 0.0
            for symbol in symbol_list:
                totalCount += cntArray[l][pos][symbol]
            if totalCount == 0: totalCount = 1.0
            for symbol in symbol_list:
                relArray[l][pos][symbol] = cntArray[l][pos][symbol] / totalCount
    dist_counters[group]['rel_' + level] = relArray

def group_total_distribution(group, group_counter, repertoires, rep_counter, level, symbol_list):
    for rep in repertoires:
        rep_id = rep['repertoire_id']
        repArray = rep_counter[rep_id][level]
        groupArray = group_counter[group][level]
        # print("repArray: \n", repArray)
        # print("GroupArray: \n", group_counter)
        for l in repArray:
            if groupArray.get(l) is None: groupArray[l] = {}
            for pos in repArray[l]:
                if groupArray[l].get(pos) is None: groupArray[l][pos] = {}
                for symbol in symbol_list:
                    if groupArray[l][pos].get(symbol) is None: groupArray[l][pos][symbol] = 0
                    groupArray[l][pos][symbol] += repArray[l][pos][symbol]
                        
def group_average_distribution(inputDict, sampleGroup, level, symbol_list):
    groups = inputDict[defaults.groupsKey]
    samples = groups[sampleGroup]['samples']

    dist_counters[sampleGroup][level + '_std'] = {}

    # average/std
    N = float(len(samples))
    for sample in samples:
        sampleArray = dist_counters[sample][level]
        groupArray = dist_counters[sampleGroup][level]
        stdArray = dist_counters[sampleGroup][level + '_std']
        for l in sampleArray:
            if groupArray.get(l) is None:
                groupArray[l] = {}
                stdArray[l] = {}
            for pos in sampleArray[l]:
                if groupArray[l].get(pos) is None:
                    groupArray[l][pos] = {}
                    stdArray[l][pos] = {}
                for symbol in symbol_list:
                    if groupArray[l][pos].get(symbol) is None:
                        groupArray[l][pos][symbol] = 0
                        stdArray[l][pos][symbol] = 0.0
                    groupArray[l][pos][symbol] += sampleArray[l][pos][symbol]
                    stdArray[l][pos][symbol] += sampleArray[l][pos][symbol] * sampleArray[l][pos][symbol]

    groupArray = dist_counters[sampleGroup][level]
    stdArray = dist_counters[sampleGroup][level + '_std']
    for l in groupArray:
        for pos in groupArray[l]:
            for symbol in symbol_list:
                if N > 1: stdArray[l][pos][symbol] = math.sqrt((stdArray[l][pos][symbol] - groupArray[l][pos][symbol] * groupArray[l][pos][symbol] / N) / (N - 1.0))
                else: stdArray[l][pos][symbol] = 0.0
                groupArray[l][pos][symbol] /= N

#
# shared/unique sequence operations
#
cdr3_shared = { "aa": {}, "nucleotide": {} }
# accumulates count of number of unique CDR3s within a metadata group
# if sublevel is defined then counts are per sublevel
# for metadata repertoire groups, counts indicate sharing level between repertoires in the group
share_summary = {}
share_summary_cdr3 = {}

def level_share_summary(inputDict, cdr3_shared_level, share_summary_level, share_summary_cdr3_level):
    # generate summary for each group
    # repertoires are just a count
    # groups produce an array with sharing counts among repertoires
    groups = inputDict.get(defaults.groups_key)
#    if groups is None:
#        return
    # iterate over each CDR3
    #print(len(cdr3_shared_level.keys()))
    for key in cdr3_shared_level.keys():
        group_set = cdr3_shared_level[key].keys()
        #print(group_set)
        # iterate across each repertoire/group
        for group in group_set:
            entry = share_summary_level.get(group)
            #print(entry, group)
#            if not groups.get(group):
            if (groups is None) or (not groups.get(group)):
                # individual repertoire
                if entry is None:
                    share_summary_level[group] = 0
                    share_summary_cdr3_level[group] = {}
                share_summary_level[group] = share_summary_level[group] + 1
                share_summary_cdr3_level[group][key] = cdr3_shared_level[key][group]
            else:
                # repertoire group
                if entry is None:
                    share_summary_level[group] = []
                    share_summary_cdr3_level[group] = {}
                repertoires = groups[group]['repertoires']
                count = 0
                for rep in repertoires:
                    if rep['repertoire_id'] in group_set: count += 1
                if count >= len(share_summary_level[group]): 
                    for i in range(len(share_summary_level[group]), count+1, 1): share_summary_level[group].append(0)
                share_summary_level[group][count] = share_summary_level[group][count] + 1
                if share_summary_cdr3_level[group].get(count) is None: share_summary_cdr3_level[group][count] = {}
                share_summary_cdr3_level[group][count][key] = cdr3_shared_level[key][group]

def generate_share_summary(inputDict, level, sublevel):
    #print(level, sublevel)
    if share_summary.get(level) is None:
        share_summary[level] = {}
        share_summary_cdr3[level] = {}
    if sublevel is None:
        level_share_summary(inputDict, cdr3_shared[level], share_summary[level], share_summary_cdr3[level])
    else:
        if share_summary[level].get(sublevel) is None:
            share_summary[level][sublevel] = {}
            share_summary_cdr3[level][sublevel] = {}
        for sl in cdr3_shared[level][sublevel]:
            if share_summary[level][sublevel].get(sl) is None:
                share_summary[level][sublevel][sl] = {}
                share_summary_cdr3[level][sublevel][sl] = {}
            level_share_summary(inputDict, cdr3_shared[level][sublevel][sl], share_summary[level][sublevel][sl], share_summary_cdr3[level][sublevel][sl])

def write_share_summary(inputDict, metadataDict, outputSpec, level, sublevel):
    cdr3_shared_level = cdr3_shared[level]
    share_summary_level = share_summary[level]
    share_summary_cdr3_level = share_summary_cdr3[level]

    fileTxt = level
    if sublevel is not None:
        cdr3_shared_level = cdr3_shared_level[sublevel]
        share_summary_level = share_summary_level[sublevel]
        share_summary_cdr3_level = share_summary_cdr3_level[sublevel]
        fileTxt = level + '_' + sublevel

    #print(level, sublevel)
    #print(len(cdr3_shared_level.keys()))
    #print(len(share_summary_level.keys()))
    #print(len(share_summary_cdr3_level.keys()))

    cdr3_text = 'junction_aa'
    if level == 'nucleotide': cdr3_text = 'junction'
    if sublevel == 'nucleotide': cdr3_text = 'junction'

    groups = inputDict.get(defaults.groups_key)
    if groups is None:
        groups = []
    filename = "summary_cdr3_" + fileTxt + "_sharing.tsv"
    writer = open(filename, 'w')

    # output specification for process metadata
    # TODO: Cheat with app name of RepCalc, we should use the app name in process metadata
    #if (not outputSpec['groups'].get("RepCalc")): outputSpec['groups']["RepCalc"] = {}
    #outputSpec['groups']['RepCalc']['cdr3_shared'] = { "files": "RepCalc_cdr3_shared", "type": "output" }
    #if (not outputSpec['files'].get("RepCalc_cdr3_shared")): outputSpec['files']["RepCalc_cdr3_shared"] = {}
    #outputSpec['files']["RepCalc_cdr3_shared"]["summary_cdr3_" + fileTxt] = { "value": filename, "description":"CDR3 Summary", "type":"tsv" }

    # inner function for summary group values
    def write_summary(group, summary_level):
        if summary_level.get(group) is None:
            return
        if groups.get(group) is not None:
            for i in range(1, len(summary_level[group]), 1): 
                writer.write('\t' + str(summary_level[group][i]))
        else:
            writer.write('\t' + str(summary_level[group]))
        writer.write('\n')

    # write intra-group summary
    if sublevel is None:
        writer.write('GROUP\tUNIQUE\tSHARING_LEVELS\n')
    else:
        writer.write('GROUP\tLEVEL\tUNIQUE\tSHARING_LEVELS\n')
    for group in groups:
        if sublevel is None:
            writer.write(group)
            write_summary(group, share_summary_level)
        else:
            for sl in share_summary_level:
                if share_summary_level[sl].get(group) is None: continue
                writer.write(group + '\t' + sl)
                write_summary(group, share_summary_level[sl])
    writer.close()

    # write detail for repertoire groups
    for group in groups:
        filename = group + ".group.cdr3_" + fileTxt + "_sharing.tsv"
        writer = open(filename, 'w')

        # output specification for process metadata
        #if (not outputSpec['files'].get(group + "_cdr3_shared")): outputSpec['files'][group + "_cdr3_shared"] = {}
        #outputSpec['groups'][group]['cdr3_shared'] = { "files": group + "_cdr3_shared", "type": "output" }
        #outputSpec['files'][group + "_cdr3_shared"][fileTxt] = { "value": filename, "description":"CDR3 Detail", "type":"tsv" }

        writer.write(cdr3_text)
        if sublevel is not None: writer.write("\tlevel")
        writer.write("\trepertoires\tsequence_count\tduplicate_count\n")

        # inner function for detail group values
        def write_group_detail(sl, summary_level, summary_cdr3_level, cdr3_level):
            if summary_cdr3_level.get(group) is None:
                return
            for i in summary_cdr3_level[group].keys():
                #print(i)
                cdr3_set = summary_cdr3_level[group][i]
                for cdr3 in cdr3_set:
                    writer.write(cdr3)
                    if sl: writer.write('\t' + sl)
                    group_set = cdr3_level[cdr3].keys()
                    repertoires = groups[group]['repertoires']
                    #print(cdr3_level[cdr3])
                    #print(group_set)
                    #print(samples)
                    repSet = []
                    for rep in repertoires:
                        if rep['repertoire_id'] in group_set: repSet.append(rep['repertoire_id'])
                    writer.write('\t' + ','.join(repSet))
                    for rep in repSet:
                        writer.write('\t' + str(cdr3_level[cdr3][rep]['count']))
                        writer.write('\t' + str(cdr3_level[cdr3][rep]['total_count']))
                    writer.write('\n')

        if sublevel is None:
            write_group_detail(None, share_summary_level, share_summary_cdr3_level, cdr3_shared_level)
        else:
            for sl in share_summary_cdr3_level:
                if share_summary_cdr3_level[sl].get(group) is None: continue
                write_group_detail(sl, share_summary_level[sl], share_summary_cdr3_level[sl], cdr3_shared_level[sl])
        writer.close()

    # write detail for repertoires
    for rep_id in metadataDict:
        filename = rep_id + ".repertoire.cdr3_" + fileTxt + "_sharing.tsv"
        writer = open(filename, 'w')
        writer.write(cdr3_text)
        if sublevel is not None: writer.write("\tlevel")
        writer.write("\tsequence_count\tduplicate_count\n")

        # inner function for detail group values
        def write_detail(sl, summary_level, summary_cdr3_level, cdr3_level):
            cdr3_set = summary_cdr3_level[rep_id]
            for cdr3 in cdr3_set:
                writer.write(cdr3)
                if sl: writer.write('\t' + sl)
                writer.write('\t' + str(cdr3_level[cdr3][rep_id]['count']))
                writer.write('\t' + str(cdr3_level[cdr3][rep_id]['total_count']) + '\n')

        if sublevel is None:
            if share_summary_cdr3_level.get(rep_id) is not None:
                write_detail(None, share_summary_level, share_summary_cdr3_level, cdr3_shared_level)
        else:
            for sl in share_summary_cdr3_level:
                if share_summary_cdr3_level[sl].get(rep_id) is not None:
                    write_detail(sl, share_summary_level[sl], share_summary_cdr3_level[sl], cdr3_shared_level[sl])
        writer.close()

def read_share_detail(inputDict, group, level, sublevel):
    groups = inputDict.get(defaults.groups_key)
    if groups is None:
        groups = {}

    fileTxt = level
    if sublevel is not None:
        fileTxt = level + '_' + sublevel

    if groups.get(group) is not None:
        filename = group + ".group.cdr3_" + fileTxt + "_sharing.tsv"
        with open(filename, 'rt') as infile:
            header = infile.readline().rstrip('\n')
            while True:
                line = infile.readline().rstrip('\n')
                if not line: break
                fields = line.split('\t')
                cdr3 = fields[0]
                #print(cdr3)
                if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                if sublevel is not None:
                    samples = fields[2].split(',')
                    countStart = 3
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    sublevel_value = fields[1]
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(sublevel_value) is None: cdr3_shared[level][sublevel][sublevel_value] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][sublevel_value]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = {}
                        cdr3_entry = cdr3_shared_level[cdr3]
                else:
                    samples = fields[1].split(',')
                    countStart = 2
                    cdr3_entry = cdr3_shared[level].get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared[level][cdr3] = {}
                        cdr3_entry = cdr3_shared[level].get(cdr3)
                i = 0
                field_count = 0
                field_total = 0
                for sample in samples:
                    field_count += int(fields[countStart + i])
                    field_total += int(fields[countStart + i + 1])
                    sample_entry = cdr3_entry.get(sample)
                    if sample_entry is None:
                        cdr3_entry[sample] = { 'count': int(fields[countStart + i]), 'total_count': int(fields[countStart + i + 1]) }
                    i += 2
                group_entry = cdr3_entry.get(group)
                if group_entry is None:
                    cdr3_entry[group] = { 'count': int(field_count), 'total_count': int(field_total) }
                #print(cdr3_entry)
                #print(field_count)
                #print(field_total)
    else:
        filename = group + ".repertoire.cdr3_" + fileTxt + "_sharing.tsv"
        with open(filename, 'rt') as infile:
            header = infile.readline().rstrip('\n')
            while True:
                line = infile.readline().rstrip('\n')
                if not line: break
                fields = line.split('\t')
                cdr3 = fields[0]
                if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                if sublevel is not None:
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    sublevel_value = fields[1]
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(sublevel_value) is None: cdr3_shared[level][sublevel][sublevel_value] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][sublevel_value]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = {}
                        cdr3_entry = cdr3_shared_level[cdr3]
                    field_count = fields[2]
                    field_total = fields[3]
                else:
                    cdr3_entry = cdr3_shared[level].get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared[level][cdr3] = {}
                        cdr3_entry = cdr3_shared[level].get(cdr3)
                    field_count = fields[1]
                    field_total = fields[2]

                group_entry = cdr3_entry.get(group)
                if group_entry is None:
                    cdr3_entry[group] = { 'count': int(field_count), 'total_count': int(field_total) }

def repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, level, sublevel):
    singleGroups = []
    for rep_id in metadataDict:
        singleGroups.append(rep_id)

    fileTxt = level
    if sublevel is not None: fileTxt = level + '_' + sublevel

    cdr3_text = 'junction_aa'
    if level == 'nucleotide': cdr3_text = 'junction'
    if sublevel == 'nucleotide': cdr3_text = 'junction'

    # single groups are directly compared to each other
    singleShared = None
    singleDiff = None
    filename1 = "repertoire.shared.matrix.cdr3_" + fileTxt + "_sharing.tsv"
    filename2 = "repertoire.diff.matrix.cdr3_" + fileTxt + "_sharing.tsv"
    writer1 = open(filename1, 'w')
    writer1.write('SHARED')
    writer2 = open(filename2, 'w')
    writer2.write('DIFF')

    numSingle = len(singleGroups)
    for j in range(0, numSingle, 1):
        colGroup = singleGroups[j]
        writer1.write('\t' + colGroup)
        writer2.write('\t' + colGroup)
    writer1.write('\n')
    writer2.write('\n')
    singleShared = numpy.zeros([numSingle, numSingle])
    singleDiff = numpy.zeros([numSingle, numSingle])
    for i in range(0, numSingle, 1):
        rowGroup = singleGroups[i]
        writer1.write(rowGroup)
        writer2.write(rowGroup)
        for j in range(0, numSingle, 1):
            colGroup = singleGroups[j]
            if sublevel is None:
                if share_summary_cdr3[level].get(rowGroup) is None:
                    A = set()
                else:
                    A = set(share_summary_cdr3[level][rowGroup].keys())
                if share_summary_cdr3[level].get(colGroup) is None:
                    B = set()
                else:
                    B = set(share_summary_cdr3[level][colGroup].keys())
            else:
                A = set()
                B = set()
                for gene in share_summary_cdr3[level][sublevel]:
                    if share_summary_cdr3[level][sublevel][gene].get(rowGroup):
                        for cdr3 in share_summary_cdr3[level][sublevel][gene][rowGroup]:
                            A.add(gene + '|' + cdr3)
                for gene in share_summary_cdr3[level][sublevel]:
                    if share_summary_cdr3[level][sublevel][gene].get(colGroup):
                        for cdr3 in share_summary_cdr3[level][sublevel][gene][colGroup]:
                            B.add(gene + '|' + cdr3)
            singleShared[i,j] = len(A & B)
            singleDiff[i,j] = len(A - B)
            #print(singleShared[i,j])
            #print(singleDiff[i,j])
            writer1.write('\t' + str(int(singleShared[i,j])))
            writer2.write('\t' + str(int(singleDiff[i,j])))
        writer1.write('\n')
        writer2.write('\n')
    writer1.close()
    writer2.close()

def repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, level, sublevel):
    singleGroups = []
    for rep_id in metadataDict:
        singleGroups.append(rep_id)

    fileTxt = level
    if sublevel is not None: fileTxt = level + '_' + sublevel

    cdr3_text = 'junction_aa'
    if level == 'nucleotide': cdr3_text = 'junction'
    if sublevel == 'nucleotide': cdr3_text = 'junction'

    # single groups are directly compared to each other
    singleShared = None
    singleDiff = None
    filename1 = "repertoire.shared.detail.cdr3_" + fileTxt + "_sharing.tsv"
    filename2 = "repertoire.diff.detail.cdr3_" + fileTxt + "_sharing.tsv"
    writer1 = open(filename1, 'w')
    writer2 = open(filename2, 'w')
    if sublevel is None:
        writer1.write(cdr3_text + '\tGROUP_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\tCOUNT_B\tTOTAL_COUNT_B\n')
        writer2.write(cdr3_text + '\tGROUP_A\tGROUP_B\tCOUNT_A\tTOTAL_COUNT_A\n')
    else:
        writer1.write(cdr3_text + '\tGROUP_A\tLEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\tLEVEL_B\tCOUNT_B\tTOTAL_COUNT_B\n')
        writer2.write(cdr3_text + '\tGROUP_A\tLEVEL_A\tGROUP_B\tLEVEL_B\tCOUNT_A\tTOTAL_COUNT_A\n')

    numSingle = len(singleGroups)
    for i in range(0, numSingle, 1):
        rowGroup = singleGroups[i]
        if sublevel is None:
            for j in range(0, numSingle, 1):
                if i == j: continue
                colGroup = singleGroups[j]
                if share_summary_cdr3[level].get(rowGroup) is None:
                    A = set()
                else:
                    A = set(share_summary_cdr3[level][rowGroup].keys())
                if share_summary_cdr3[level].get(colGroup) is None:
                    B = set()
                else:
                    B = set(share_summary_cdr3[level][colGroup].keys())
                singleShared = A & B
                singleDiff = A - B
                if len(singleShared) > 0:
                    for cdr3 in singleShared:
                        writer1.write(cdr3)
                        writer1.write('\t' + rowGroup)
                        writer1.write('\t' + str(share_summary_cdr3[level][rowGroup][cdr3]['count']))
                        writer1.write('\t' + str(share_summary_cdr3[level][rowGroup][cdr3]['total_count']))
                        writer1.write('\t' + colGroup)
                        writer1.write('\t' + str(share_summary_cdr3[level][colGroup][cdr3]['count']))
                        writer1.write('\t' + str(share_summary_cdr3[level][colGroup][cdr3]['total_count']))
                        writer1.write('\n')
                if len(singleDiff) > 0:
                    for cdr3 in singleDiff:
                        writer2.write(cdr3)
                        writer2.write('\t' + rowGroup)
                        writer2.write('\t' + colGroup)
                        writer2.write('\t' + str(share_summary_cdr3[level][rowGroup][cdr3]['count']))
                        writer2.write('\t' + str(share_summary_cdr3[level][rowGroup][cdr3]['total_count']))
                        writer2.write('\n')
        else:
            for rowLevel in share_summary_cdr3[level][sublevel]:
                if not share_summary_cdr3[level][sublevel][rowLevel].get(rowGroup): continue
                for j in range(0, numSingle, 1):
                    if i == j: continue
                    colGroup = singleGroups[j]
                    #if metadata.sample_contains_file(inputDict, rowGroup, colGroup): continue
                    if not share_summary_cdr3[level][sublevel][rowLevel].get(colGroup): continue
                    A = set(share_summary_cdr3[level][sublevel][rowLevel][rowGroup].keys())
                    B = set(share_summary_cdr3[level][sublevel][rowLevel][colGroup].keys())
                    singleShared = A & B
                    singleDiff = A - B
                    if len(singleShared) > 0:
                        for cdr3 in singleShared:
                            writer1.write(cdr3)
                            writer1.write('\t' + rowGroup + '\t' + rowLevel)
                            writer1.write('\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][cdr3]['count']))
                            writer1.write('\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][cdr3]['total_count']))
                            writer1.write('\t' + colGroup + '\t' + rowLevel)
                            writer1.write('\t' + str(share_summary_cdr3[level][sublevel][rowLevel][colGroup][cdr3]['count']))
                            writer1.write('\t' + str(share_summary_cdr3[level][sublevel][rowLevel][colGroup][cdr3]['total_count']))
                            writer1.write('\n')
                    if len(singleDiff) > 0:
                        for cdr3 in singleDiff:
                            writer2.write(cdr3)
                            writer2.write('\t' + rowGroup + '\t' + rowLevel)
                            writer2.write('\t' + colGroup + '\t' + rowLevel)
                            writer2.write('\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][cdr3]['count']))
                            writer2.write('\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][cdr3]['total_count']))
                            writer2.write('\n')
    writer1.close()
    writer2.close()


def generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, level, sublevel):
    groups = inputDict.get(defaults.groups_key)
    if groups is None:
        return

    # pairwise group comparison
    sampleGroups = []
    for group in groups:
        sampleGroups.append(group)

    fileTxt = level
    if sublevel is not None: fileTxt = level + '_' + sublevel

    cdr3_text = 'junction_aa'
    if level == 'nucleotide': cdr3_text = 'junction'
    if sublevel == 'nucleotide': cdr3_text = 'junction'

    # output specification for process metadata
    # TODO: Cheat with app name of RepCalc, we should use the app name in process metadata
    #if (not outputSpec['groups'].get("RepCalc")): outputSpec['groups']["RepCalc"] = {}
    #outputSpec['groups']['RepCalc']['cdr3_shared'] = { "files": "RepCalc_cdr3_shared", "type": "output" }
    #if (not outputSpec['files'].get("RepCalc_cdr3_shared")): outputSpec['files']["RepCalc_cdr3_shared"] = {}
    #outputSpec['files']["RepCalc_cdr3_shared"]["group_comparison_cdr3_" + fileTxt] = { "value": filename1, "description":"Shared CDR3 Comparison", "type":"tsv" }
    #outputSpec['files']["RepCalc_cdr3_shared"]["group_diff_cdr3_" + fileTxt] = { "value": filename2, "description":"Unique CDR3 Comparison", "type":"tsv" }

    # sample groups are compared at each level

    def group_comparison():
        #filename1 = filePrefix + "_summary_comparison_cdr3_" + fileTxt + "_sharing.tsv"
        #filename2 = filePrefix + "_diff_cdr3_" + fileTxt + "_sharing.tsv"
        #filename3 = filePrefix + "_comparison_cdr3_" + fileTxt + "_sharing.tsv"
        filename1 = "group.summary_comparison.cdr3_" + fileTxt + "_sharing.tsv"
        filename2 = "group.diff.cdr3_" + fileTxt + "_sharing.tsv"
        filename3 = "group.comparison.cdr3_" + fileTxt + "_sharing.tsv"
        writer1 = open(filename1, 'w')
        writer2 = open(filename2, 'w')
        writer3 = open(filename3, 'w')

        # output specification for process metadata
        # TODO: Cheat with app name of RepCalc, we should use the app name in process metadata
        #outputSpec['files']["RepCalc_cdr3_shared"]["sampleGroup_summary_comparison_cdr3_" + fileTxt] = { "value": filename1, "description":"CDR3 Summary Comparison", "type":"tsv" }
        #outputSpec['files']["RepCalc_cdr3_shared"]["sampleGroup_diff_cdr3_" + fileTxt] = { "value": filename2, "description":"Unique CDR3 Comparison", "type":"tsv" }
        #outputSpec['files']["RepCalc_cdr3_shared"]["sampleGroup_comparison_cdr3_" + fileTxt] = { "value": filename3, "description":"Shared CDR3 Comparison", "type":"tsv" }

        #print(sampleGroups)
        #print(share_summary)
        numGroups = len(sampleGroups)
        for row_i in range(0, numGroups, 1):
            rowGroup = sampleGroups[row_i]
            numRows = len(share_summary[level][rowGroup]) - 1
            for row_j in range(0, numGroups, 1):
                if row_i == row_j: continue
                colGroup = sampleGroups[row_j]
                numCols = len(share_summary[level][colGroup]) - 1
                sharedMatrix = numpy.zeros([numRows, numCols])
                writer1.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer1.write('SHARED')
                writer2.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer2.write(cdr3_text + '\tGROUP_A\tSHARE_LEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\n')
                writer3.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer3.write(cdr3_text + '\tGROUP_A\tSHARE_LEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\tSHARE_LEVEL_B\tCOUNT_B\tTOTAL_COUNT_B\n')
                for j in range(0, numCols, 1):
                    writer1.write('\t' + str(j+1))
                writer1.write('\n')
                count_array = []
                string_array = []
                diff_count_array = []
                diff_string_array = []
                for i in range(0, numRows, 1):
                    if not share_summary_cdr3[level][rowGroup].get(i+1): A = set()
                    else: A = set(share_summary_cdr3[level][rowGroup][i+1].keys())
                    groupDiff = A
                    writer1.write(str(i+1))
                    for j in range(0, numCols, 1):
                        #print(j+1)
                        #print(share_summary_cdr3[level][colGroup][j+1])
                        if not share_summary_cdr3[level][colGroup].get(j+1): B = set()
                        else: B = set(share_summary_cdr3[level][colGroup][j+1].keys())
                        groupShared = A & B
                        groupDiff = groupDiff - B
                        sharedMatrix[i,j] = len(groupShared)
                        writer1.write('\t' + str(int(sharedMatrix[i,j])))
                        if len(groupShared) > 0:
                            for cdr3 in groupShared:
                                line = cdr3
                                line += '\t' + rowGroup + '\t' + str(i+1)
                                line += '\t' + str(share_summary_cdr3[level][rowGroup][i+1][cdr3]['count'])
                                line += '\t' + str(share_summary_cdr3[level][rowGroup][i+1][cdr3]['total_count'])
                                line += '\t' + colGroup + '\t' + str(j+1)
                                line += '\t' + str(share_summary_cdr3[level][colGroup][j+1][cdr3]['count'])
                                line += '\t' + str(share_summary_cdr3[level][colGroup][j+1][cdr3]['total_count'])
                                line += '\n'
                                string_array.append(line)
                                count_array.append(share_summary_cdr3[level][rowGroup][i+1][cdr3]['total_count'])
                    if len(groupDiff) > 0:
                        for cdr3 in groupDiff:
                            line = cdr3
                            line += '\t' + rowGroup + '\t' + str(i+1)
                            line += '\t' + str(share_summary_cdr3[level][rowGroup][i+1][cdr3]['count'])
                            line += '\t' + str(share_summary_cdr3[level][rowGroup][i+1][cdr3]['total_count'])
                            line += '\t' + colGroup + '\n'
                            diff_string_array.append(line)
                            diff_count_array.append(share_summary_cdr3[level][rowGroup][i+1][cdr3]['total_count'])
                    writer1.write('\n')
                sort_array = numpy.array(count_array).argsort()
                for idx in reversed(sort_array):
                    writer3.write(string_array[idx])
                diff_sort_array = numpy.array(diff_count_array).argsort()
                for idx in reversed(diff_sort_array):
                    writer2.write(diff_string_array[idx])
        writer1.close()
        writer2.close()
        writer3.close()

    def group_level_comparison():
        #filename1 = filePrefix + "_comparison_cdr3_" + fileTxt + "_sharing.tsv"
        #filename2 = filePrefix + "_diff_cdr3_" + fileTxt + "_sharing.tsv"
        filename1 = "group.comparison.cdr3_" + fileTxt + "_sharing.tsv"
        filename2 = "group.diff.cdr3_" + fileTxt + "_sharing.tsv"
        writer1 = open(filename1, 'w')
        writer2 = open(filename2, 'w')
 
        # output specification for process metadata
        # TODO: Cheat with app name of RepCalc, we should use the app name in process metadata
        #outputSpec['files']["RepCalc_cdr3_shared"]["sampleGroup_comparison_cdr3_" + fileTxt] = { "value": filename1, "description":"Shared CDR3 Comparison", "type":"tsv" }
        #outputSpec['files']["RepCalc_cdr3_shared"]["sampleGroup_diff_cdr3_" + fileTxt] = { "value": filename2, "description":"Unique CDR3 Comparison", "type":"tsv" }

        numGroups = len(sampleGroups)
        for row_i in range(0, numGroups, 1):
            rowGroup = sampleGroups[row_i]
            for row_j in range(0, numGroups, 1):
                if row_i == row_j: continue
                colGroup = sampleGroups[row_j]
                writer1.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer1.write(cdr3_text + '\tGROUP_A\tLEVEL_A\tSHARE_LEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\tLEVEL_B\tSHARE_LEVEL_B\tCOUNT_B\tTOTAL_COUNT_B\n')
                writer2.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer2.write(cdr3_text + '\tGROUP_A\tLEVEL_A\tSHARE_LEVEL_A\tCOUNT_A\tTOTAL_COUNT_A\tGROUP_B\n')
                count_array = []
                string_array = []
                diff_count_array = []
                diff_string_array = []
                for rowLevel in share_summary_cdr3[level][sublevel]:
                    if not share_summary_cdr3[level][sublevel][rowLevel].get(rowGroup): continue
                    numRows = len(share_summary[level][sublevel][rowLevel][rowGroup]) - 1

                    # sublevel is not in the other group at any share level
                    if not share_summary_cdr3[level][sublevel][rowLevel].get(colGroup):
                        for i in range(0, numRows, 1):
                            if not share_summary_cdr3[level][sublevel][rowLevel][rowGroup].get(i+1): continue
                            A = set(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1].keys())
                            for cdr3 in A:
                                line = cdr3
                                line += '\t' + rowGroup + '\t' + rowLevel + '\t' + str(i+1)
                                line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['count'])
                                line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['total_count'])
                                line += '\t' + colGroup + '\n'
                                diff_string_array.append(line)
                                diff_count_array.append(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['total_count'])
                        continue

                    # some level of sharing at the sublevel
                    numCols = len(share_summary[level][sublevel][rowLevel][colGroup]) - 1
                    for i in range(0, numRows, 1):
                        if not share_summary_cdr3[level][sublevel][rowLevel][rowGroup].get(i+1): continue
                        A = set(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1].keys())
                        groupDiff = A
                        for j in range(0, numCols, 1):
                            if not share_summary_cdr3[level][sublevel][rowLevel][colGroup].get(j+1): continue
                            B = set(share_summary_cdr3[level][sublevel][rowLevel][colGroup][j+1].keys())
                            groupShared = A & B
                            groupDiff = groupDiff - B
                            if len(groupShared) > 0:
                                for cdr3 in groupShared:
                                    line = cdr3
                                    line += '\t' + rowGroup + '\t' + rowLevel + '\t' + str(i+1)
                                    line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['count'])
                                    line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['total_count'])
                                    line += '\t' + colGroup + '\t' + rowLevel + '\t' + str(j+1)
                                    line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][colGroup][j+1][cdr3]['count'])
                                    line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][colGroup][j+1][cdr3]['total_count'])
                                    line += '\n'
                                    string_array.append(line)
                                    count_array.append(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['total_count'])
                        if len(groupDiff) > 0:
                            for cdr3 in groupDiff:
                                line = cdr3
                                line += '\t' + rowGroup + '\t' + rowLevel + '\t' + str(i+1)
                                line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['count'])
                                line += '\t' + str(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['total_count'])
                                line += '\t' + colGroup + '\n'
                                diff_string_array.append(line)
                                diff_count_array.append(share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1][cdr3]['total_count'])
                sort_array = numpy.array(count_array).argsort()
                for idx in reversed(sort_array):
                    writer1.write(string_array[idx])
                diff_sort_array = numpy.array(diff_count_array).argsort()
                for idx in reversed(diff_sort_array):
                    writer2.write(diff_string_array[idx])
        writer1.close()
        writer2.close()

    if len(sampleGroups) != 0:
        if sublevel is None: group_comparison()
        else: group_level_comparison()

#
# Calculation module protocol functions
# initialize_calculation_module()
# process_record()
# finalize_calculation_module()
#

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # length and AA/NT distribution operations
    for rep_id in inputDict[defaults.full_metadata_key]:
        cdr3_histograms[rep_id] = { obj: [] for obj in length_levels }
        cdr3_histograms_productive[rep_id] = { obj: [] for obj in length_levels }
        dist_counters[rep_id] = { obj: {} for obj in length_levels }
        dist_counters_productive[rep_id] = { obj: {} for obj in length_levels }

    # group counts and totals
    if inputDict.get(defaults.groups_key) is not None:
        for group in inputDict[defaults.groups_key]:
            group_cdr3_histograms[group] = {}
            group_cdr3_histograms_productive[group] = {}
            
            group_cdr3_total_histograms[group] = { obj: [] for obj in length_levels }
            group_cdr3_total_histograms_productive[group] = { obj: [] for obj in length_levels }
            
            group_dist_counters[group] = {}
            group_dist_counters_productive[group] = {}
            
            group_dist_counters_total[group] = { obj: {} for obj in length_levels }
            group_dist_counters_productive_total[group] = { obj: {} for obj in length_levels }
            
            rep_list = [rep_id['repertoire_id'] for rep_id in inputDict[defaults.groups_key][group]['repertoires']]
            # print("Repertoire List: ", group, rep_list)
            for rep_id in rep_list:
                group_cdr3_histograms[group][rep_id] = { obj: [] for obj in length_levels }
                group_cdr3_histograms_productive[group][rep_id] = { obj: [] for obj in length_levels }
                # print("Group CDR3 hist: \n", group_cdr3_histograms)  
                group_dist_counters[group][rep_id] = { obj: {} for obj in length_levels }
                group_dist_counters_productive[group][rep_id] = { obj: {} for obj in length_levels }

def process_record(inputDict, metadataDict, currentFile, calc, fields):
    """Perform calculation from given fields"""
    rep_id = fields.get('repertoire_id')
    if rep_id is None:
        print('WARNING: sequence id', fields.get('sequence_id'), 'is missing repertoire_id')
        return
    if inputDict[defaults.full_metadata_key].get(rep_id) is None:
        print('WARNING: sequence id', fields.get('sequence_id'), 'has repertoire_id not in repertoire list')
        return

    # length operations
    if lengthKey in calc['operations']:
        length = fields.get('junction_aa_length')
        if length is None:
            # some old AIRR TSV is missing field, so calculate if necessary
            fields['junction_aa_length'] = len(fields['junction_aa'])
            length = fields.get('junction_aa_length')

        # AA length
        if length is not None:
            # repertoire counts
            increment_count(cdr3_histograms[rep_id]['aa_wo_duplicates'], length)
            increment_duplicate_count(cdr3_histograms[rep_id]['aa'], length, defaults.get_duplicate_count(fields))
            if fields.get('productive'):
                increment_count(cdr3_histograms_productive[rep_id]['aa_wo_duplicates'], length)
                increment_duplicate_count(cdr3_histograms_productive[rep_id]['aa'], length, defaults.get_duplicate_count(fields))
            # group counts
            if inputDict.get(defaults.groups_key) is not None:
                groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
                if groupList:
                    for group in groupList:
                        # we do not bother optimizing if the group has no rearrangement filter
                        # just accumulate repertoire counts for each group, data is duplicated but it is small
                        if defaults.apply_filter(inputDict, group, fields):
                            increment_count(group_cdr3_histograms[group][rep_id]['aa_wo_duplicates'], length)
                            increment_duplicate_count(group_cdr3_histograms[group][rep_id]['aa'], length, defaults.get_duplicate_count(fields))
                            if fields.get('productive'):
                                increment_count(group_cdr3_histograms_productive[group][rep_id]['aa_wo_duplicates'], length)
                                increment_duplicate_count(group_cdr3_histograms_productive[group][rep_id]['aa'], length, defaults.get_duplicate_count(fields))

        # nucleotide length
        length = fields.get('junction_length')
        if length is not None:
            # repertoire counts
            increment_count(cdr3_histograms[rep_id]['nucleotide_wo_duplicates'], length)
            increment_duplicate_count(cdr3_histograms[rep_id]['nucleotide'], length, defaults.get_duplicate_count(fields))
            if fields.get('productive'):
                increment_count(cdr3_histograms_productive[rep_id]['nucleotide_wo_duplicates'], length)
                increment_duplicate_count(cdr3_histograms_productive[rep_id]['nucleotide'], length, defaults.get_duplicate_count(fields))
            # group counts
            if inputDict.get(defaults.groups_key) is not None:
                groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
                if groupList:
                    for group in groupList:
                        # we do not bother optimizing if the group has no rearrangement filter
                        # just accumulate repertoire counts for each group, data is duplicated but it is small
                        if defaults.apply_filter(inputDict, group, fields):
                            increment_count(group_cdr3_histograms[group][rep_id]['nucleotide_wo_duplicates'], length)
                            increment_duplicate_count(group_cdr3_histograms[group][rep_id]['nucleotide'], length, defaults.get_duplicate_count(fields))
                            if fields.get('productive'):
                                increment_count(group_cdr3_histograms_productive[group][rep_id]['nucleotide_wo_duplicates'], length)
                                increment_duplicate_count(group_cdr3_histograms_productive[group][rep_id]['nucleotide'], length, defaults.get_duplicate_count(fields))

    # AA/NT distribution operations
    if distributionKey in calc['operations']:
        cdr3 = fields.get('junction_aa')
        if cdr3 is not None:
            aa_distribution(dist_counters[rep_id]['aa'], cdr3)
            if fields.get('productive'):
                aa_distribution(dist_counters_productive[rep_id]['aa'], cdr3)
            # group aa distribution
            if inputDict.get(defaults.groups_key) is not None:
                groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
                if groupList:
                    for group in groupList:
                        # we do not bother optimizing if the group has no rearrangement filter
                        # just accumulate repertoire counts for each group, data is duplicated but it is small
                        if defaults.apply_filter(inputDict, group, fields):
                            aa_distribution(group_dist_counters[group][rep_id]['aa'], cdr3)
                            if fields.get('productive'):
                                aa_distribution(group_dist_counters_productive[group][rep_id]['aa'], cdr3)
        cdr3 = fields.get('junction')
        if cdr3 is not None:
            aa_distribution(dist_counters[rep_id]['nucleotide'], cdr3)
            if fields.get('productive'):
                aa_distribution(dist_counters_productive[rep_id]['nucleotide'], cdr3)
            # group aa distribution
            if inputDict.get(defaults.groups_key) is not None:
                groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
                if groupList:
                    for group in groupList:
                        # we do not bother optimizing if the group has no rearrangement filter
                        # just accumulate repertoire counts for each group, data is duplicated but it is small
                        if defaults.apply_filter(inputDict, group, fields):
                            aa_distribution(group_dist_counters[group][rep_id]['nucleotide'], cdr3)
                            if fields.get('productive'):
                                aa_distribution(group_dist_counters_productive[group][rep_id]['nucleotide'], cdr3)

    # share/unique sequence operations
    if sharedKey in calc['operations']:
        germline = inputDict.get(defaults.germline_key)
        groups = metadata.groupsWithRepertoire(inputDict, rep_id)

        if "productive" in calc['filters']:
            if not fields.get('productive'):
                return

        for l in calc['levels']:
            cdr3_entry = None
            sublevel = None
            level = None
            # get appropriate entry
            if l == 'aa':
                cdr3 = fields.get('junction_aa')
                if cdr3 is not None and len(cdr3) > 0:
                    cdr3_entry = cdr3_shared[l].get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared[l][cdr3] = {}
                        cdr3_entry = cdr3_shared[l].get(cdr3)

            elif l == 'nucleotide':
                cdr3 = fields.get('junction')
                if cdr3 is not None and len(cdr3) > 0:
                    cdr3_entry = cdr3_shared[l].get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared[l][cdr3] = {}
                        cdr3_entry = cdr3_shared[l].get(cdr3)

            elif l == 'v,aa':
                level = 'v'
                sublevel = 'aa'
                cdr3 = fields.get('junction_aa')
                if cdr3 is not None and len(cdr3) > 0:
                    vcall = fields.get('v_call')
                    if not vcall or len(vcall) == 0: continue
                    vgene = gldb.getDisplayName(germline, vcall, "gene")
                    if not vgene: continue
                    vgene = vgene[0]
                    if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(vgene) is None: cdr3_shared[level][sublevel][vgene] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][vgene]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = {}
                        cdr3_entry = cdr3_shared_level[cdr3]

            elif l == 'v,nucleotide':
                level = 'v'
                sublevel = 'nucleotide'
                cdr3 = fields.get('junction')
                if cdr3 is not None and len(cdr3) > 0:
                    vcall = fields.get('v_call')
                    if not vcall or len(vcall) == 0: continue
                    vgene = gldb.getDisplayName(germline, vcall, "gene")
                    if not vgene: continue
                    vgene = vgene[0]
                    if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(vgene) is None: cdr3_shared[level][sublevel][vgene] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][vgene]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = {}
                        cdr3_entry = cdr3_shared_level[cdr3]

            elif l == 'vj,aa':
                level = 'vj'
                sublevel = 'aa'
                cdr3 = fields.get('junction_aa')
                if cdr3 is not None and len(cdr3) > 0:
                    vcall = fields.get('v_call')
                    if not vcall or len(vcall) == 0: continue
                    vgene = gldb.getDisplayName(germline, vcall, "gene")
                    if not vgene: continue
                    vgene = vgene[0]
                    jcall = fields.get('j_call')
                    if not jcall or len(jcall) == 0: continue
                    jgene = gldb.getDisplayName(germline, jcall, "gene")
                    if not jgene: continue
                    jgene = jgene[0]
                    combo = vgene + '|' + jgene
                    if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(combo) is None: cdr3_shared[level][sublevel][combo] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][combo]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = {}
                        cdr3_entry = cdr3_shared_level[cdr3]

            elif l == 'vj,nucleotide':
                level = 'vj'
                sublevel = 'nucleotide'
                cdr3 = fields.get('junction')
                if cdr3 is not None and len(cdr3) > 0:
                    vcall = fields.get('v_call')
                    if not vcall or len(vcall) == 0: continue
                    vgene = gldb.getDisplayName(germline, vcall, "gene")
                    if not vgene: continue
                    vgene = vgene[0]
                    jcall = fields.get('j_call')
                    if not jcall or len(jcall) == 0: continue
                    jgene = gldb.getDisplayName(germline, jcall, "gene")
                    if not jgene: continue
                    jgene = jgene[0]
                    combo = vgene + '|' + jgene
                    if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(combo) is None: cdr3_shared[level][sublevel][combo] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][combo]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = {}
                        cdr3_entry = cdr3_shared_level[cdr3]

            else:
                # unknown level so just ignore
                continue
            # increment counts for entry
            if cdr3_entry is not None:
                # print("cdr3_entry: ", cdr3_entry.get(rep_id))
                if cdr3_entry.get(rep_id) is None:
                    cdr3_entry[rep_id] = { 'count': 0, 'total_count': 0 }
                cdr3_entry[rep_id]['count'] = cdr3_entry[rep_id]['count'] + 1
                cdr3_entry[rep_id]['total_count'] = cdr3_entry[rep_id]['total_count'] + defaults.get_duplicate_count(fields)
               
                if inputDict.get(defaults.groups_key) is not None:
                    groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
                    if groupList:
                        for group in groupList:
                            if defaults.has_rearrangement_filter(inputDict, group):
                                if defaults.apply_filter(inputDict, group, fields):
                                    if cdr3_entry.get(group) is None:
                                        cdr3_entry[group] = {}
                                    if cdr3_entry[group].get(rep_id) is None:
                                        cdr3_entry[group][rep_id] = { 'count': 0, 'total_count': 0 }
                                    cdr3_entry[group][rep_id]['count'] = cdr3_entry[group][rep_id]['count'] + 1
                                    cdr3_entry[group][rep_id]['total_count'] = cdr3_entry[group][rep_id]['total_count'] + defaults.get_duplicate_count(fields)
               
                # for group in groups:
                #     group_entry = cdr3_entry.get(group)
                #     if group_entry is None:
                #         # cdr3_entry[group] = { 'count': 0, 'total_count': 0 }
                #         cdr3_entry[group] = {}
                #         # group_entry = cdr3_entry.get(group)
                #     # group_entry['count'] = group_entry['count'] + 1
                #     # group_entry['total_count'] = group_entry['total_count'] + defaults.get_duplicate_count(fields)


def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    stage = ''
    if inputDict.get(defaults.processing_stage_key):
        stage = '.' + inputDict.get(defaults.processing_stage_key)

    # length operations
    if lengthKey in calc['operations']:
        # repertoire counts
        # compute relative frequency and output repertoire values
        for rep_id in metadataDict:
            for level in length_levels:
                compute_relative(cdr3_histograms, rep_id, level)
                compute_relative(cdr3_histograms_productive, rep_id, level)
                output_counts(rep_id + stage + ".junction_" + level + "_length.tsv", cdr3_histograms[rep_id][level], cdr3_histograms[rep_id]['rel_' + level])
                output_counts(rep_id + stage + ".productive.junction_" + level + "_length.tsv", cdr3_histograms_productive[rep_id][level], cdr3_histograms_productive[rep_id]['rel_' + level])

        # group counts
        if inputDict.get(defaults.groups_key) is not None:
            # accumulate totals for group across all repertoires in group
            for group in inputDict[defaults.groups_key]:
                rep_list = [rep_id['repertoire_id'] for rep_id in inputDict[defaults.groups_key][group]['repertoires']]
                for rep_id in rep_list:
                    for level in length_levels:
                        # accumulate totals for group across all repertoires in group
                        group_counts(group_cdr3_total_histograms[group][level], group_cdr3_histograms[group][rep_id][level])
                        group_counts(group_cdr3_total_histograms_productive[group][level], group_cdr3_histograms_productive[group][rep_id][level])
                        # compute relative frequency for each repertoire in group
                        compute_relative(group_cdr3_histograms[group], rep_id, level)
                        compute_relative(group_cdr3_histograms_productive[group], rep_id, level)

                for level in length_levels:
                    # compute relative frequency for group values
                    compute_relative(group_cdr3_total_histograms, group, level)
                    compute_relative(group_cdr3_total_histograms_productive, group, level)
                    # output group distribution
                    output_counts(group + stage + ".group.junction_" + level + "_length.tsv", group_cdr3_total_histograms[group][level], group_cdr3_total_histograms[group]['rel_' + level])
                    output_counts(group + stage + ".group.productive.junction_" + level + "_length.tsv", group_cdr3_total_histograms_productive[group][level], group_cdr3_total_histograms_productive[group]['rel_' + level])
                    # output distribution for each repertoire in a group
                    output_counts_table(group + stage + ".repertoire.count.junction_" + level + "_length.csv", group_cdr3_total_histograms[group][level], inputDict[defaults.groups_key][group]['repertoires'], group_cdr3_histograms[group], level)
                    output_counts_table(group + stage + ".repertoire.frequency.junction_" + level + "_length.csv", group_cdr3_total_histograms[group][level], inputDict[defaults.groups_key][group]['repertoires'], group_cdr3_histograms[group], 'rel_' + level)
                    output_counts_table(group + stage + ".repertoire.count.productive.junction_" + level + "_length.csv", group_cdr3_total_histograms_productive[group][level], inputDict[defaults.groups_key][group]['repertoires'], group_cdr3_histograms_productive[group], level)
                    output_counts_table(group + stage + ".repertoire.frequency.productive.junction_" + level + "_length.csv", group_cdr3_total_histograms_productive[group][level], inputDict[defaults.groups_key][group]['repertoires'], group_cdr3_histograms_productive[group], 'rel_' + level)


    # AA/NT distribution operations
    if distributionKey in calc['operations']:
        #  group_dist_counters_total
        #  group_dist_counters_productive_total
        # group counts
        if inputDict.get(defaults.groups_key) is not None:
            for group in inputDict[defaults.groups_key]:
                group_total_distribution(group, group_dist_counters_total, inputDict[defaults.groups_key][group]['repertoires'], group_dist_counters[group], 'aa', aa_list)
                group_total_distribution(group, group_dist_counters_productive_total, inputDict[defaults.groups_key][group]['repertoires'], group_dist_counters_productive[group], 'aa', aa_list)
                group_total_distribution(group, group_dist_counters_total, inputDict[defaults.groups_key][group]['repertoires'], group_dist_counters[group], 'nucleotide', nt_list)
                group_total_distribution(group, group_dist_counters_productive_total, inputDict[defaults.groups_key][group]['repertoires'], group_dist_counters_productive[group], 'nucleotide', nt_list)
                
                output_distribution(group + stage + ".junction_aa_distribution.tsv", group_dist_counters_total[group]['aa'], aa_list)
                output_distribution(group + stage + ".productive.junction_aa_distribution.tsv", group_dist_counters_productive_total[group]['aa'], aa_list)
                output_distribution(group + stage + ".junction_nt_distribution.tsv", group_dist_counters_total[group]['nucleotide'], nt_list)
                output_distribution(group + stage + ".productive.junction_nt_distribution.tsv", group_dist_counters_productive_total[group]['nucleotide'], nt_list)

        # repertoire counts
        for rep_id in metadataDict:
            output_distribution(rep_id + stage + ".junction_aa_distribution.tsv", dist_counters[rep_id]['aa'], aa_list)
            output_distribution(rep_id + stage + ".productive.junction_aa_distribution.tsv", dist_counters_productive[rep_id]['aa'], aa_list)
            output_distribution(rep_id + stage + ".junction_nt_distribution.tsv", dist_counters[rep_id]['nucleotide'], nt_list)
            output_distribution(rep_id + stage + ".productive.junction_nt_distribution.tsv", dist_counters_productive[rep_id]['nucleotide'], nt_list)

        # relative calculations
#        for group in groups:
#            compute_relative_distribution(group, 'aa', aa_list)
#            compute_relative_distribution(group, 'nucleotide', nt_list)

        # write output
#        for group in groups:
            # output specification for process metadata
#            if (not outputSpec['files'].get(group + "_cdr3_distribution")): outputSpec['files'][group + "_cdr3_distribution"] = {}
#            outputSpec['groups'][group]['cdr3_distribution'] = { "files": group + "_cdr3_distribution", "type": "output" }

    # share/unique sequence operations
    if sharedKey in calc['operations']:
        #print(cdr3_shared)
        for level in calc['levels']:
            if level == 'aa':
                generate_share_summary(inputDict, level, None)
                write_share_summary(inputDict, metadataDict, outputSpec, level, None)
                #generate_share_comparison(inputDict, outputSpec, level, None)
            if level == 'nucleotide':
                generate_share_summary(inputDict, level, None)
                write_share_summary(inputDict, metadataDict, outputSpec, level, None)
                #generate_share_comparison(inputDict, outputSpec, level, None)
            if level == 'v,aa':
                generate_share_summary(inputDict, 'v', 'aa')
                write_share_summary(inputDict, metadataDict, outputSpec, 'v', 'aa')
                #generate_share_comparison(inputDict, outputSpec, 'v', 'aa')
            if level == 'v,nucleotide':
                generate_share_summary(inputDict, 'v', 'nucleotide')
                write_share_summary(inputDict, metadataDict, outputSpec, 'v', 'nucleotide')
                #generate_share_comparison(inputDict, outputSpec, 'v', 'nucleotide')
            if level == 'vj,aa':
                generate_share_summary(inputDict, 'vj', 'aa')
                write_share_summary(inputDict, metadataDict, outputSpec, 'vj', 'aa')
                #generate_share_comparison(inputDict, outputSpec, 'vj', 'aa')
            if level == 'vj,nucleotide':
                generate_share_summary(inputDict, 'vj', 'nucleotide')
                write_share_summary(inputDict, metadataDict, outputSpec, 'vj', 'nucleotide')
                #generate_share_comparison(inputDict, outputSpec, 'vj', 'nucleotide')

    # shared/unique pairwise comparison
    if compareKey in calc['operations']:
        groups = inputDict.get(defaults.groups_key)
        filePrefix = None
#        if groups is None:
#            return
#        filePrefix = 'group'
#        if len(groups) == 2: filePrefix = '_'.join(groups.keys())
        for level in calc['levels']:
            if level == 'aa':
                if groups:
                    for group in groups:
                        read_share_detail(inputDict, group, level, None)
                for rep_id in metadataDict:
                    read_share_detail(inputDict, rep_id, level, None)
                #print(cdr3_shared)
                generate_share_summary(inputDict, level, None)
                repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, level, None)
                repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, level, None)
                generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, level, None)
            if level == 'nucleotide':
                if groups:
                    for group in groups:
                        read_share_detail(inputDict, group, level, None)
                for rep_id in metadataDict:
                    read_share_detail(inputDict, rep_id, level, None)
                generate_share_summary(inputDict, level, None)
                repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, level, None)
                repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, level, None)
                generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, level, None)
            if level == 'v,aa':
                if groups:
                    for group in groups:
                        read_share_detail(inputDict, group, 'v', 'aa')
                for rep_id in metadataDict:
                    read_share_detail(inputDict, rep_id, 'v', 'aa')
                generate_share_summary(inputDict, 'v', 'aa')
                repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, 'v', 'aa')
                repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, 'v', 'aa')
                generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, 'v', 'aa')
            if level == 'v,nucleotide':
                if groups:
                    for group in groups:
                        read_share_detail(inputDict, group, 'v', 'nucleotide')
                for rep_id in metadataDict:
                    read_share_detail(inputDict, rep_id, 'v', 'nucleotide')
                generate_share_summary(inputDict, 'v', 'nucleotide')
                repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, 'v', 'nucleotide')
                repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, 'v', 'nucleotide')
                generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, 'v', 'nucleotide')
            if level == 'vj,aa':
                if groups:
                    for group in groups:
                        read_share_detail(inputDict, group, 'vj', 'aa')
                for rep_id in metadataDict:
                    read_share_detail(inputDict, rep_id, 'vj', 'aa')
                generate_share_summary(inputDict, 'vj', 'aa')
                repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, 'vj', 'aa')
                repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, 'vj', 'aa')
                generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, 'vj', 'aa')
            if level == 'vj,nucleotide':
                if groups:
                    for group in groups:
                        read_share_detail(inputDict, group, 'vj', 'nucleotide')
                for rep_id in metadataDict:
                    read_share_detail(inputDict, rep_id, 'vj', 'nucleotide')
                generate_share_summary(inputDict, 'vj', 'nucleotide')
                repertoire_share_matrix(inputDict, metadataDict, outputSpec, filePrefix, 'vj', 'nucleotide')
                repertoire_share_detail(inputDict, metadataDict, outputSpec, filePrefix, 'vj', 'nucleotide')
                generate_share_comparison(inputDict, metadataDict, outputSpec, filePrefix, 'vj', 'nucleotide')
