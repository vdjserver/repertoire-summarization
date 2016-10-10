"""
Gene segment calculation module
"""

# repsum modules
from .version import __version__
import utils
import defaults
import metadata
import gldb

segment_counters = {}

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # TODO: should check for required columns

    # setup data structures for groups
    groups = inputDict[defaults.groupsKey]
    for group in groups: segment_counters[group] = utils.IncrementMapWrapper()

    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    for group in groupSet:
        gene = fields[headerMapping[defaults.headerNames['V_CALL']]]
        if gene is not None: segment_counters[group].increment(gene)
        gene = fields[headerMapping[defaults.headerNames['D_CALL']]]
        if gene is not None: segment_counters[group].increment(gene)
        gene = fields[headerMapping[defaults.headerNames['J_CALL']]]
        if gene is not None: segment_counters[group].increment(gene)

def finalize_calculation_module(inputDict, metadataDict, outputSpec):
    """Finalize and save the calculations"""
    outputSpec['groups'] = {}
    groups = inputDict[defaults.groupsKey]
    for group in groups:
        print("group: " + group)
        if (not outputSpec['groups'].get(group)): outputSpec['groups'][group] = {}
        if (not outputSpec['files'].get(group + "_gene_segment_usage")): outputSpec['files'][group + "_gene_segment_usage"] = {}
        segment_counters[group].printMap()
        organism = 'human'
        hierarchy = gldb.getHierarchyBy(organism)
        JSON = utils.jsonify_hierarchy(hierarchy,organism,segment_counters[group].get_map(),"value")
        filename = group + "_segment_counts.json"
        outputSpec['groups'][group]['gene_segment_usage'] = { "files": group + "_gene_segment_usage", "type": "output" }
        outputSpec['files'][group + "_gene_segment_usage"]['absolute_counts'] = filename
        writer = open(filename, 'w')
        writer.write(JSON)
        writer.close()
