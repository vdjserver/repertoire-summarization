"""
Gene segment calculation module
"""

# repsum modules
from .version import __version__
import utils
import defaults
import metadata
import gldb
import json

segment_counters = {}

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    groups = inputDict[defaults.groupsKey]
    for group in groups: segment_counters[group] = {}

    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
# # UNCOMMENT FOR DUPCOUNT
#
#    if defaults.headerNames['DUPCOUNT'] in headerMapping:	# Check if template counts are available.
#        for group in groupSet:
#            for headerName in ['V_CALL', 'D_CALL', 'J_CALL']:
#                gene = fields[headerMapping[defaults.headerNames[headerName]]]
#                if gene is not None:
#                    if gene in segment_counters[group]: segment_counters[group][gene] += fields[headerMapping[defaults.headerNames['DUPCOUNT']]]
#                    else: segment_counters[group][gene] = fields[headerMapping[defaults.headerNames['DUPCOUNT']]]
#    else:
#        for group in groupSet:
#            for headerName in ['V_CALL', 'D_CALL', 'J_CALL']:
#                gene = fields[headerMapping[defaults.headerNames[headerName]]]
#                if gene is not None:
#                    if gene in segment_counters[group]: segment_counters[group][gene] += 1
#                    else: segment_counters[group][gene] = 1
    """Perform calculation from given fields"""
    for group in groupSet:
        for headerName in ['V_CALL', 'D_CALL', 'J_CALL']:
            gene = fields[headerMapping[defaults.headerNames[headerName]]]
            if gene is not None:
                if gene in segment_counters[group]: segment_counters[group][gene] += 1
                else: segment_counters[group][gene] = 1

def makeTree(hierarchy, segment_counters):
    siblings = []
    total_count = 0
    for label, sub_hierarchy in hierarchy.items():
        if len(sub_hierarchy) > 0:
            children = makeTree(sub_hierarchy, segment_counters)
            count = 0
            for child in children:
                count += child['count']
            total_count += count
            siblings.append({
                'label': label,
                'children': children,
                'count': count,
            })
        else:
            if label in segment_counters:
                count = segment_counters[label]
                total_count += count
                siblings.append({
                    'label': label,
                    'count': count,
                })
            else:
                siblings.append({
                    'label': label,
                    'count': 0,
                })
    for sibling in siblings:
        if total_count > 0:
            sibling['fraction'] = float(sibling['count'])/float(total_count)
        else:
            sibling['fraction'] = 0.0
    return siblings

def finalize_calculation_module(inputDict, metadataDict, outputSpec):
    """Finalize and save the calculations"""
    groups = inputDict[defaults.groupsKey]
    for group in groups:
        print("group: " + group)
        organism = 'human'
        hierarchy = gldb.getHierarchyBy(organism)
        tree = makeTree({'human': hierarchy}, segment_counters[group])
        JSON = json.dumps(tree[0], indent=2)
        filename = group + "_segment_counts.json"
        writer = open(filename, 'w')
        writer.write(JSON)
        writer.close()

