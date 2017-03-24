"""
CDR3 calculation module
"""

# repsum modules
from .version import __version__
import defaults
import metadata
import gldb
import json
import math
import numpy

# module operations
lengthKey = "length"
sharedKey = "shared"

gene_hierarchy = None
invert_hierarchy = {}

cdr3_histograms = {}

def increment_count(anArray, sequence):
    if len(sequence) == 0: return
    if len(sequence) >= len(anArray):
        for i in range(len(anArray), len(sequence)+1, 1): anArray.append(0.0)
    anArray[len(sequence)] += 1

def group_average(inputDict, sampleGroup, level):
    groups = inputDict[defaults.groupsKey]
    samples = groups[sampleGroup]['samples']

    cdr3_histograms[sampleGroup][level + '_std'] = []

    # average/std
    N = float(len(samples))
    for sample in samples:
        sampleName = metadata.sample_for_uuid(inputDict, sample)
        #print(sampleGroup, sampleName)

        sampleArray = cdr3_histograms[sampleName][level]
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

def compute_relative(group, level):
    lenArray = cdr3_histograms[group][level]
    totalCount = 0.0
    for i in range(0, len(lenArray), 1): totalCount += lenArray[i]
    if totalCount == 0: totalCount = 1.0
    relArray = []
    for i in range(0, len(lenArray), 1): relArray.append(lenArray[i] / totalCount)
    cdr3_histograms[group]['rel_' + level] = relArray

#
# shared/unique sequence operations
#
cdr3_shared = { "aa": {}, "nucleotide": {} }
share_summary = {}
share_summary_cdr3 = {}

def level_share_summary(inputDict, cdr3_shared_level, share_summary_level, share_summary_cdr3_level):
    # generate summary for each group
    # samples and files are just a count
    # sample groups produce an array with sharing counts among samples
    groups = inputDict[defaults.groupsKey]
    for key in cdr3_shared_level.keys():
        group_set = cdr3_shared_level[key]
        for group in group_set:
            entry = share_summary_level.get(group)
            if entry is None:
                if (groups[group]['type'] == 'sampleGroup'):
                    share_summary_level[group] = []
                    share_summary_cdr3_level[group] = {}
                else:
                    share_summary_level[group] = 0
                    share_summary_cdr3_level[group] = set()
            if (groups[group]['type'] == 'sampleGroup'):
                samples = groups[group]['samples']
                count = 0
                for sample in samples:
                    sampleName = metadata.sample_for_uuid(inputDict, sample)
                    if sampleName in group_set: count += 1
                if count >= len(share_summary_level[group]): 
                    for i in range(len(share_summary_level[group]), count+1, 1): share_summary_level[group].append(0)
                share_summary_level[group][count] = share_summary_level[group][count] + 1
                if share_summary_cdr3_level[group].get(count) is None: share_summary_cdr3_level[group][count] = set()
                share_summary_cdr3_level[group][count].add(key)
            else:
                share_summary_level[group] = share_summary_level[group] + 1
                share_summary_cdr3_level[group].add(key)

def generate_share_summary(inputDict, level, sublevel):
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

def write_share_summary(inputDict, level, sublevel):
    cdr3_shared_level = cdr3_shared[level]
    share_summary_level = share_summary[level]
    share_summary_cdr3_level = share_summary_cdr3[level]
    fileTxt = level
    if sublevel is not None:
        cdr3_shared_level = cdr3_shared_level[sublevel]
        share_summary_level = share_summary_level[sublevel]
        share_summary_cdr3_level = share_summary_cdr3_level[sublevel]
        fileTxt = level + '_' + sublevel

    cdr3_text = defaults.headerNames['CDR3_AA']
    if level == 'nucleotide': cdr3_text = defaults.headerNames['CDR3_SEQ']
    if sublevel == 'nucleotide': cdr3_text = defaults.headerNames['CDR3_SEQ']

    groups = inputDict[defaults.groupsKey]
    writer = open("summary_cdr3_" + fileTxt + "_sharing.tsv", 'w')

    # inner function for summary group values
    def write_summary(group, summary_level):
        if (groups[group]['type'] == 'sampleGroup'):
            for i in range(1, len(summary_level[group]), 1): writer.write('\t' + str(summary_level[group][i]))
        else:
            writer.write('\t' + str(summary_level[group]))
        writer.write('\n')

    # write summary
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

    # write detail
    for group in groups:
        filename = group + "_cdr3_" + fileTxt + "_sharing.tsv"
        writer = open(filename, 'w')
        writer.write(cdr3_text)
        if sublevel is not None: writer.write("\tLEVEL")
        if (groups[group]['type'] == 'sampleGroup'):
            writer.write("\tSAMPLES\tCOUNT\n")
        else:
            writer.write("\tCOUNT\n")

        # inner function for detail group values
        def write_detail(sl, summary_level, summary_cdr3_level, cdr3_level):
            if (groups[group]['type'] == 'sampleGroup'):
                for i in range(1, len(summary_cdr3_level[group]), 1):
                    cdr3_set = summary_cdr3_level[group][i]
                    for cdr3 in cdr3_set:
                        writer.write(cdr3)
                        if sl: writer.write('\t' + sl)
                        group_set = cdr3_level[cdr3]
                        samples = groups[group]['samples']
                        sampleSet = []
                        for sample in samples:
                            sampleName = metadata.sample_for_uuid(inputDict, sample)
                            if sampleName in group_set: sampleSet.append(sampleName)
                        writer.write('\t' + ','.join(sampleSet))
                        for sample in sampleSet: writer.write('\t' + str(0))
                        writer.write('\n')
            else:
                cdr3_set = summary_cdr3_level[group]
                for cdr3 in cdr3_set:
                    writer.write(cdr3)
                    if sl: writer.write('\t' + sl)
                    writer.write('\t' + str(0) + '\n')

        if sublevel is None:
            write_detail(None, share_summary_level, share_summary_cdr3_level, cdr3_shared_level)
        else:
            for sl in share_summary_cdr3_level:
                if share_summary_cdr3_level[sl].get(group) is None: continue
                write_detail(sl, share_summary_level[sl], share_summary_cdr3_level[sl], cdr3_shared_level[sl])

        writer.close()

def generate_share_comparison(inputDict, level, sublevel):
    # pairwise group comparison
    groups = inputDict[defaults.groupsKey]
    sampleGroups = []
    singleGroups = []
    for group in groups:
        if (groups[group]['type'] == 'sampleGroup'): sampleGroups.append(group)
        else: singleGroups.append(group)

    fileTxt = level
    if sublevel is not None: fileTxt = level + '_' + sublevel

    cdr3_text = defaults.headerNames['CDR3_AA']
    if level == 'nucleotide': cdr3_text = defaults.headerNames['CDR3_SEQ']
    if sublevel == 'nucleotide': cdr3_text = defaults.headerNames['CDR3_SEQ']

    # single groups are directly compared to each other
    singleShared = None
    singleDiff = None
    writer1 = open("group_comparison_cdr3_" + fileTxt + "_sharing.tsv", 'w')
    writer2 = open("group_diff_cdr3_" + fileTxt + "_sharing.tsv", 'w')

    def single_comparison():
        writer1.write('SHARED')
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
                A = share_summary_cdr3[level][rowGroup]
                B = share_summary_cdr3[level][colGroup]
                singleShared[i,j] = len(A & B)
                singleDiff[i,j] = len(A - B)
                writer1.write('\t' + str(int(singleShared[i,j])))
                writer2.write('\t' + str(int(singleDiff[i,j])))
            writer1.write('\n')
            writer2.write('\n')
    
    def single_sublevel_comparison():
        writer1.write('GROUP_A\tLEVEL_A\tGROUP_B\tLEVEL_B\t' + cdr3_text + '\tCOUNT\n')
        writer2.write('GROUP_A\tLEVEL_A\tGROUP_B\tLEVEL_B\tDIFF\n')
        numSingle = len(singleGroups)
        for i in range(0, numSingle, 1):
            rowGroup = singleGroups[i]
            for rowLevel in share_summary_cdr3[level][sublevel]:
                if not share_summary_cdr3[level][sublevel][rowLevel].get(rowGroup): continue
                for j in range(0, numSingle, 1):
                    colGroup = singleGroups[j]
                    if metadata.sample_contains_file(inputDict, rowGroup, colGroup): continue
                    for colLevel in share_summary_cdr3[level][sublevel]:
                        if not share_summary_cdr3[level][sublevel][colLevel].get(colGroup): continue
                        A = share_summary_cdr3[level][sublevel][rowLevel][rowGroup]
                        B = share_summary_cdr3[level][sublevel][colLevel][colGroup]
                        singleShared = A & B
                        singleDiff = A - B
                        if len(singleShared) > 0:
                            for cdr3 in singleShared:
                                writer1.write(rowGroup + '\t' + rowLevel + '\t' + colGroup + '\t' + colLevel)
                                writer1.write('\t' + cdr3 + '\t' + str(0))
                                writer1.write('\n')
                            writer2.write(rowGroup + '\t' + rowLevel + '\t' + colGroup + '\t' + colLevel)
                            writer2.write('\t' + str(len(singleDiff)))
                            writer2.write('\n')

    if sublevel is None: single_comparison()
    else: single_sublevel_comparison()
    writer1.close()
    writer2.close()

    # sample groups are compared at each level
    writer1 = open("sampleGroup_comparison_cdr3_" + fileTxt + "_sharing.tsv", 'w')
    writer2 = open("sampleGroup_diff_cdr3_" + fileTxt + "_sharing.tsv", 'w')

    def group_comparison():
        numGroups = len(sampleGroups)
        for row_i in range(0, numGroups, 1):
            rowGroup = sampleGroups[row_i]
            numRows = len(share_summary[level][rowGroup]) - 1
            for row_j in range(0, numGroups, 1):
                if row_i == row_j: continue
                colGroup = sampleGroups[row_j]
                numCols = len(share_summary[level][colGroup]) - 1
                sharedMatrix = numpy.zeros([numRows, numCols])
                diffMatrix = numpy.zeros([numRows, numCols])
                writer1.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer1.write('SHARED')
                writer2.write('GROUPS\t' + rowGroup + '\t' + colGroup + '\n')
                writer2.write('DIFF')
                for j in range(0, numCols, 1):
                    writer1.write('\t' + str(j+1))
                    writer2.write('\t' + str(j+1))
                writer1.write('\n')
                writer2.write('\n')
                for i in range(0, numRows, 1):
                    writer1.write(str(i+1))
                    writer2.write(str(i+1))
                    for j in range(0, numCols, 1):
                        A = share_summary_cdr3[level][rowGroup][i+1]
                        B = share_summary_cdr3[level][colGroup][j+1]
                        sharedMatrix[i,j] = len(A & B)
                        diffMatrix[i,j] = len(A - B)
                        writer1.write('\t' + str(int(sharedMatrix[i,j])))
                        writer2.write('\t' + str(int(diffMatrix[i,j])))
                    writer1.write('\n')
                    writer2.write('\n')

    def group_level_comparison():
        writer1.write('GROUP_A\tLEVEL_A\tSHARE_LEVEL_A\tGROUP_B\tLEVEL_B\tSHARE_LEVEL_B\t' + cdr3_text + '\tCOUNT\n')
        writer2.write('GROUP_A\tLEVEL_A\tGROUP_B\tLEVEL_B\tDIFF\n')
        numGroups = len(sampleGroups)
        for row_i in range(0, numGroups, 1):
            rowGroup = sampleGroups[row_i]
            for rowLevel in share_summary_cdr3[level][sublevel]:
                if not share_summary_cdr3[level][sublevel][rowLevel].get(rowGroup): continue
                numRows = len(share_summary[level][sublevel][rowLevel][rowGroup]) - 1
                for row_j in range(0, numGroups, 1):
                    if row_i == row_j: continue
                    colGroup = sampleGroups[row_j]
                    for colLevel in share_summary_cdr3[level][sublevel]:
                        if not share_summary_cdr3[level][sublevel][colLevel].get(colGroup): continue
                        numCols = len(share_summary[level][sublevel][colLevel][colGroup]) - 1
                        for i in range(0, numRows, 1):
                            for j in range(0, numCols, 1):
                                if not share_summary_cdr3[level][sublevel][rowLevel][rowGroup].get(i+1): continue
                                if not share_summary_cdr3[level][sublevel][colLevel][colGroup].get(j+1): continue
                                A = share_summary_cdr3[level][sublevel][rowLevel][rowGroup][i+1]
                                B = share_summary_cdr3[level][sublevel][colLevel][colGroup][j+1]
                                groupShared = A & B
                                gropuDiff = A - B
                                if len(groupShared) > 0:
                                    for cdr3 in groupShared:
                                        writer1.write(rowGroup + '\t' + rowLevel + '\t' + str(i+1))
                                        writer1.write('\t' + colGroup + '\t' + colLevel + '\t' + str(j+1))
                                        writer1.write('\t' + cdr3 + '\t' + str(0))
                                        writer1.write('\n')
                                    
    if sublevel is None: group_comparison()
    else: group_level_comparison()
    writer1.close()
    writer2.close()

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # length operations
    groups = inputDict[defaults.groupsKey]
    for group in groups: cdr3_histograms[group] = { 'aa': [], 'nucleotide': [] }

    gene_hierarchy = gldb.getHierarchyBy(inputDict[defaults.organismKey])
    gldb.invertHierarchy(gene_hierarchy, invert_hierarchy, None, None)
    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    # length operations
    if lengthKey in calc['operations']:
        groups = inputDict[defaults.groupsKey]
        for group in groupSet:
            if (groups[group]['type'] == 'sampleGroup'): continue
            cdr3 = fields[headerMapping[defaults.headerNames['CDR3_AA']]]
            if cdr3 is not None: increment_count(cdr3_histograms[group]['aa'], cdr3)
            cdr3 = fields[headerMapping[defaults.headerNames['CDR3_SEQ']]]
            if cdr3 is not None: increment_count(cdr3_histograms[group]['nucleotide'], cdr3)

    # share/unique sequence operations
    if sharedKey in calc['operations']:
        for l in calc['levels']:
            sublevel = None
            level = None
            if l == 'aa':
                cdr3 = fields[headerMapping[defaults.headerNames['CDR3_AA']]]
                if cdr3 is not None and len(cdr3) > 0:
                    cdr3_entry = cdr3_shared[l].get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared[l][cdr3] = set()
                        cdr3_entry = cdr3_shared[l].get(cdr3)
                    groups = inputDict[defaults.groupsKey]
                    for group in groupSet: cdr3_entry.add(group)
            elif l == 'nucleotide':
                cdr3 = fields[headerMapping[defaults.headerNames['CDR3_SEQ']]]
                if cdr3 is not None and len(cdr3) > 0:
                    cdr3_entry = cdr3_shared[l].get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared[l][cdr3] = set()
                        cdr3_entry = cdr3_shared[l].get(cdr3)
                    groups = inputDict[defaults.groupsKey]
                    for group in groupSet: cdr3_entry.add(group)
            elif l == 'v,aa':
                level = 'v'
                sublevel = 'aa'
                cdr3 = fields[headerMapping[defaults.headerNames['CDR3_AA']]]
                if cdr3 is not None and cdr3!= 'None' and len(cdr3) > 0:
                    vcall = fields[headerMapping[defaults.headerNames['V_CALL']]]
                    if not vcall or vcall == 'None': continue
                    vgene = invert_hierarchy[vcall]
                    if not vgene: continue
                    vgene = vgene.keys()[0]
                    if not vgene: continue
                    if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(vgene) is None: cdr3_shared[level][sublevel][vgene] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][vgene]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = set()
                        cdr3_entry = cdr3_shared_level[cdr3]
                    groups = inputDict[defaults.groupsKey]
                    for group in groupSet: cdr3_entry.add(group)
            elif l == 'v,nucleotide':
                level = 'v'
                sublevel = 'nucleotide'
                cdr3 = fields[headerMapping[defaults.headerNames['CDR3_SEQ']]]
                if cdr3 is not None and cdr3 != 'None' and len(cdr3) > 0:
                    vcall = fields[headerMapping[defaults.headerNames['V_CALL']]]
                    if not vcall or vcall == 'None': continue
                    #print(fields)
                    vgene = invert_hierarchy[vcall]
                    if not vgene: continue
                    vgene = vgene.keys()[0]
                    if not vgene: continue
                    if cdr3_shared.get(level) is None: cdr3_shared[level] = {}
                    if cdr3_shared[level].get(sublevel) is None: cdr3_shared[level][sublevel] = {}
                    if cdr3_shared[level][sublevel].get(vgene) is None: cdr3_shared[level][sublevel][vgene] = {}
                    cdr3_shared_level = cdr3_shared[level][sublevel][vgene]
                    cdr3_entry = cdr3_shared_level.get(cdr3)
                    if cdr3_entry is None:
                        cdr3_shared_level[cdr3] = set()
                        cdr3_entry = cdr3_shared_level[cdr3]
                    groups = inputDict[defaults.groupsKey]
                    for group in groupSet: cdr3_entry.add(group)
#            elif level == 'v,nucleotide':
#            elif level == 'vj,aa':

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    # length operations
    if lengthKey in calc['operations']:
        groups = inputDict[defaults.groupsKey]
        for group in groups:
            compute_relative(group, 'aa')
            compute_relative(group, 'nucleotide')
        for group in groups:
            # output specification for process metadata
            if (not outputSpec['files'].get(group + "_cdr3_length")): outputSpec['files'][group + "_cdr3_length"] = {}
            outputSpec['groups'][group]['cdr3_length'] = { "files": group + "_cdr3_length", "type": "output" }

            # aa histogram
            filename = group + "_cdr3_aa_length.tsv"
            writer = open(filename, 'w')
            cntArray = cdr3_histograms[group]['aa']
            relArray = cdr3_histograms[group]['rel_aa']
            if (groups[group]['type'] == 'sampleGroup'):
                group_average(inputDict, group, 'aa')
                group_average(inputDict, group, 'rel_aa')
                stdArray = cdr3_histograms[group]['aa_std']
                stdRelArray = cdr3_histograms[group]['rel_aa_std']
                writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_COUNT_STD\tCDR3_RELATIVE\tCDR3_RELATIVE_STD\n')
                for i in range(0, len(cntArray), 1):
                    writer.write(str(i) + '\t' + str(cntArray[i]) + '\t' + str(stdArray[i])
                                 + '\t' + str(relArray[i]) + '\t' + str(stdRelArray[i]) + '\n')
            else:
                writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_RELATIVE\n')
                for i in range(0, len(cntArray), 1):
                    writer.write(str(i) + '\t' + str(cntArray[i]) + '\t' + str(relArray[i]) + '\n')
            writer.close()
            outputSpec['files'][group + "_cdr3_length"]['aa'] = { "value": filename, "description":"CDR3 AA Length Histogram", "type":"tsv" }

            # nucleotide histogram
            filename = group + "_cdr3_nucleotide_length.tsv"
            writer = open(filename, 'w')
            cntArray = cdr3_histograms[group]['nucleotide']
            relArray = cdr3_histograms[group]['rel_nucleotide']
            if (groups[group]['type'] == 'sampleGroup'):
                group_average(inputDict, group, 'nucleotide')
                group_average(inputDict, group, 'rel_nucleotide')
                stdArray = cdr3_histograms[group]['nucleotide_std']
                stdRelArray = cdr3_histograms[group]['rel_nucleotide_std']
                writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_COUNT_STD\tCDR3_RELATIVE\tCDR3_RELATIVE_STD\n')
                for i in range(0, len(cntArray), 1):
                    writer.write(str(i) + '\t' + str(cntArray[i]) + '\t' + str(stdArray[i])
                                 + '\t' + str(relArray[i]) + '\t' + str(stdRelArray[i]) + '\n')
            else:
                writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_RELATIVE\n')
                for i in range(0, len(cntArray), 1):
                    writer.write(str(i) + '\t' + str(cntArray[i]) + '\t' + str(relArray[i]) + '\n')
            writer.close()
            outputSpec['files'][group + "_cdr3_length"]['nucleotide'] = { "value": filename, "description":"CDR3 Nucleotide Length Histogram", "type":"tsv" }

    # share/unique sequence operations
    if sharedKey in calc['operations']:
        #print(cdr3_shared)
        for level in calc['levels']:
            if level == 'aa':
                generate_share_summary(inputDict, level, None)
                write_share_summary(inputDict, level, None)
                generate_share_comparison(inputDict, level, None)
            if level == 'nucleotide':
                generate_share_summary(inputDict, level, None)
                write_share_summary(inputDict, level, None)
                generate_share_comparison(inputDict, level, None)
            if level == 'v,aa':
                generate_share_summary(inputDict, 'v', 'aa')
                write_share_summary(inputDict, 'v', 'aa')
                generate_share_comparison(inputDict, 'v', 'aa')
            if level == 'v,nucleotide':
                generate_share_summary(inputDict, 'v', 'nucleotide')
                write_share_summary(inputDict, 'v', 'nucleotide')
                generate_share_comparison(inputDict, 'v', 'nucleotide')
