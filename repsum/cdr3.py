"""
CDR3 calculation module
"""

# repsum modules
from .version import __version__
import defaults
import metadata
import gldb
import json

cdr3_histograms = {}

def incrementCount(anArray, sequence):
    if len(sequence) == 0: return
    if len(sequence) >= len(anArray):
        for i in range(len(anArray), len(sequence)+1, 1): anArray.append(0)
    anArray[len(sequence)] += 1

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    groups = inputDict[defaults.groupsKey]
    for group in groups: cdr3_histograms[group] = { 'aa': [], 'nucleotide': [] }
    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    cdr3 = fields[headerMapping[defaults.headerNames['CDR3_AA']]]
    if cdr3 is not None:
        for group in groupSet: incrementCount(cdr3_histograms[group]['aa'], cdr3)
    cdr3 = fields[headerMapping[defaults.headerNames['CDR3_SEQ']]]
    if cdr3 is not None:
        for group in groupSet: incrementCount(cdr3_histograms[group]['nucleotide'], cdr3)

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    groups = inputDict[defaults.groupsKey]
    for group in groups:
        # output specification for process metadata
        if (not outputSpec['files'].get(group + "_cdr3_length")): outputSpec['files'][group + "_cdr3_length"] = {}
        outputSpec['groups'][group]['cdr3_length'] = { "files": group + "_cdr3_length", "type": "output" }

        # aa histogram
        filename = group + "_cdr3_aa_length.tsv"
        writer = open(filename, 'w')
        lenArray = cdr3_histograms[group]['aa']
        totalCount = 0.0
        for i in range(0, len(lenArray), 1): totalCount += lenArray[i]
        if totalCount == 0: totalCount = 1
        writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_RELATIVE\n')
        for i in range(0, len(lenArray), 1): writer.write(str(i) + '\t' + str(lenArray[i]) + '\t' + str(lenArray[i] / totalCount) + '\n')
        writer.close()
        outputSpec['files'][group + "_cdr3_length"]['aa'] = { "value": filename, "description":"CDR3 AA Length Histogram", "type":"tsv" }

        # nucleotide histogram
        filename = group + "_cdr3_nucleotide_length.tsv"
        writer = open(filename, 'w')
        lenArray = cdr3_histograms[group]['nucleotide']
        totalCount = 0.0
        for i in range(0, len(lenArray), 1): totalCount += lenArray[i]
        if totalCount == 0: totalCount = 1
        writer.write('CDR3_LENGTH\tCDR3_COUNT\tCDR3_RELATIVE\n')
        for i in range(0, len(lenArray), 1): writer.write(str(i) + '\t' + str(lenArray[i]) + '\t' + str(lenArray[i] / totalCount) + '\n')
        writer.close()
        outputSpec['files'][group + "_cdr3_length"]['nucleotide'] = { "value": filename, "description":"CDR3 Nucleotide Length Histogram", "type":"tsv" }
