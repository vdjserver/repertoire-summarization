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
import numpy

# module operations
usageKey = "usage"
comboKey = "combo"

segment_counters = {}
combo_counters = {}
summary_combo_counters = {}

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    groups = inputDict[defaults.groupsKey]
    for group in groups:
        segment_counters[group] = {}
        combo_counters[group] = {}
    return

def process_record(inputDict, metadataDict, currentFile, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    groups = inputDict[defaults.groupsKey]
    for group in groupSet: 
        if (groups[group]['type'] == 'sampleGroup'): continue
        for headerName in ['V_CALL', 'D_CALL', 'J_CALL']:
            gene = fields[headerMapping[defaults.headerNames[headerName]]]
            if gene is not None:
                if gene in segment_counters[group]: segment_counters[group][gene] += defaults.get_dupcount(headerMapping, fields)
                else: segment_counters[group][gene] = defaults.get_dupcount(headerMapping, fields)

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
                if combo_counters[group][level] is None: continue
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
            if combo_counters[group] is None: continue
            if combo_counters[group][level] is None: continue
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
            row_entry = combo_counters[rowGroup][level]
            A = set(row_entry.keys())
            string_array = []
            count_array = []
            diff_string_array = []
            diff_count_array = []
            for colGroup in singleGroups:
                if metadata.sample_contains_file(inputDict, rowGroup, colGroup): continue
                col_entry = combo_counters[colGroup][level]
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
    groups = inputDict[defaults.groupsKey]
    hierarchy = gldb.getHierarchyBy(inputDict[defaults.organismKey])
    for group in groups:
        if (groups[group]['type'] == 'sampleGroup'): continue
        tree = make_tree({inputDict[defaults.organismKey]: hierarchy}, segment_counters[group], 0)[0]
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

            tree = group_tree(inputDict, group, sampleTrees.keys(), {inputDict[defaults.organismKey]: hierarchy}, sampleTrees)[0]
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

    if comboKey in calc['operations']:
        generate_combo_summary(inputDict, outputSpec, calc)
        generate_combo_detail(inputDict, outputSpec, calc)
        generate_combo_comparison(inputDict, outputSpec, calc)
                    
