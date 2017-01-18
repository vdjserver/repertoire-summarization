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

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    groups = inputDict[defaults.groupsKey]
    for group in groups: cdr3_histograms[group] = { 'aa': [], 'nucleotide': [] }
    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    groups = inputDict[defaults.groupsKey]
    for group in groupSet:
        if (groups[group]['type'] == 'sampleGroup'): continue        
        cdr3 = fields[headerMapping[defaults.headerNames['CDR3_AA']]]
        if cdr3 is not None: increment_count(cdr3_histograms[group]['aa'], cdr3)
        cdr3 = fields[headerMapping[defaults.headerNames['CDR3_SEQ']]]
        if cdr3 is not None: increment_count(cdr3_histograms[group]['nucleotide'], cdr3)

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
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
