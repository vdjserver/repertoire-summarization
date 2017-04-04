"""
Gene segment calculation module
"""

# repsum modules
from .version import __version__
import defaults
import metadata
import gldb
import json
import math

segment_counters = {}

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    groups = inputDict[defaults.groupsKey]
    for group in groups: segment_counters[group] = {}
    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    groups = inputDict[defaults.groupsKey]
    for group in groupSet: 
        if (groups[group]['type'] == 'sampleGroup'): continue
        for headerName in ['V_CALL', 'D_CALL', 'J_CALL']:
            gene = fields[headerMapping[defaults.headerNames[headerName]]]
            if gene is not None:
                if gene in segment_counters[group]: segment_counters[group][gene] += defaults.get_dupcount(headerMapping, fields)
                else: segment_counters[group][gene] = defaults.get_dupcount(headerMapping, fields)

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

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    groups = inputDict[defaults.groupsKey]
    organism = 'human'
    hierarchy = gldb.getHierarchyBy(organism)
    for group in groups:
        if (groups[group]['type'] == 'sampleGroup'): continue
        tree = make_tree({organism: hierarchy}, segment_counters[group], 0)[0]
        segment_counters[group]['TREE'] = tree

    for group in groups:
        #print("group: " + group)
        if (groups[group]['type'] == 'sampleGroup'):
            samples = groups[group]['samples']
            sampleTrees = {}
            for sample in samples:
                sampleName = metadata.sample_for_uuid(inputDict, sample)
                #print(group, sampleName)
                sampleTrees[sampleName] = [ segment_counters[sampleName]['TREE'] ]

            tree = group_tree(inputDict, group, sampleTrees.keys(), {organism: hierarchy}, sampleTrees)[0]
            segment_counters[group]['TREE'] = tree

#        if 'absolute' not in calc[defaults.calcOpsKey]:
#            tree = prune_tree(tree, 'absolute')
#        if 'relative' not in calc[defaults.calcOpsKey]:
#            tree = prune_tree(tree, 'relative')

        tree = segment_counters[group]['TREE']
        JSON = json.dumps(tree, indent=2)
        filename = group + "_segment_counts.json"
        writer = open(filename, 'w')
        writer.write(JSON)
        writer.close()
        # output specification for process metadata
        if (not outputSpec['files'].get(group + "_gene_segment_usage")): outputSpec['files'][group + "_gene_segment_usage"] = {}
        outputSpec['groups'][group]['gene_segment_usage'] = { "files": group + "_gene_segment_usage", "type": "output" }
        outputSpec['files'][group + "_gene_segment_usage"]['counts'] = { "value": filename, "description":"Gene Segment Usage", "type":"json" }


