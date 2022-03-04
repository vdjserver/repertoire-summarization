"""
Mutational analysis
"""

#
# mutations.py
# Mutational analysis
#
# VDJServer Analysis Portal
# Repertoire calculations and comparison
# https://vdjserver.org
#
# Copyright (C) 2022 The University of Texas Southwestern Medical Center
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
#import repcalc.gldb as gldb
import json
import math
import numpy
import csv
import airr

# module operations
countKey = "count"
module_levels = [ "repertoire", "clone" ]
#module_levels = [ "repertoire", "clone", "allele", "gene", "subgroup" ]

mutation_counts = {}
group_mutation_counts = {}

total_names = ['seq_count','total_count','mu_unique_r','mu_unique_s']
count_names = total_names + ['mu_count_r_aa','mu_count_r_nt','mu_count_s_aa','mu_count_s_nt']
fwr1_names = ['fwr1_r_nt','fwr1_r_aa','fwr1_s_nt','fwr1_s_aa']
cdr1_names = ['cdr1_r_nt','cdr1_r_aa','cdr1_s_nt','cdr1_s_aa']
fwr2_names = ['fwr2_r_nt','fwr2_r_aa','fwr2_s_nt','fwr2_s_aa']
cdr2_names = ['cdr2_r_nt','cdr2_r_aa','cdr2_s_nt','cdr2_s_aa']
fwr3_names = ['fwr3_r_nt','fwr3_r_aa','fwr3_s_nt','fwr3_s_aa']
region_names = fwr1_names + cdr1_names + fwr2_names + cdr2_names + fwr3_names
pos_names = []
for i in range(1,105):
    f = 'mu_count_' + str(i) + '_r'
    pos_names.append(f + '_aa')
    pos_names.append(f + '_nt')
    f = 'mu_count_' + str(i) + '_s'
    pos_names.append(f + '_aa')
    pos_names.append(f + '_nt')
transfer_names = count_names + region_names + pos_names
group_names = []
for field in transfer_names:
    group_names.append(field + '_avg')
    group_names.append(field + '_std')

def add_mutation_count(counts, row):
    """Add mutation counts to counters"""
    duplicate_count = defaults.get_duplicate_count(row)
    # total counts
    counts['seq_count'] += 1
    counts['total_count'] += duplicate_count

    # position counts
    for i in range(1,105):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts[f + '_aa'] += duplicate_count
            counts['mu_count_r_aa'] += duplicate_count
            counts[f + '_nt'] += int(row[f]) * duplicate_count
            counts['mu_count_r_nt'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts[f + '_aa'] += duplicate_count
            counts['mu_count_s_aa'] += duplicate_count
            counts[f + '_nt'] += int(row[f]) * duplicate_count
            counts['mu_count_s_nt'] += int(row[f]) * duplicate_count

    # region counts
    for i in range(1,27):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['fwr1_r_aa'] += duplicate_count
            counts['fwr1_r_nt'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['fwr1_s_aa'] += duplicate_count
            counts['fwr1_s_nt'] += int(row[f]) * duplicate_count

    for i in range(27,39):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['cdr1_r_aa'] += duplicate_count
            counts['cdr1_r_nt'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['cdr1_s_aa'] += duplicate_count
            counts['cdr1_s_nt'] += int(row[f]) * duplicate_count

    for i in range(39,56):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['fwr2_r_aa'] += duplicate_count
            counts['fwr2_r_nt'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['fwr2_s_aa'] += duplicate_count
            counts['fwr2_s_nt'] += int(row[f]) * duplicate_count

    for i in range(56,66):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['cdr2_r_aa'] += duplicate_count
            counts['cdr2_r_nt'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['cdr2_s_aa'] += duplicate_count
            counts['cdr2_s_nt'] += int(row[f]) * duplicate_count

    for i in range(66,105):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['fwr3_r_aa'] += duplicate_count
            counts['fwr3_r_nt'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['fwr3_s_aa'] += duplicate_count
            counts['fwr3_s_nt'] += int(row[f]) * duplicate_count

def add_repertoire_count(rep_id, row):
    # initialize counts per repertoire
    rep_counts = mutation_counts['repertoire'].get(rep_id)
    if rep_counts is None:
        mutation_counts['repertoire'][rep_id] = { 'repertoire_id':rep_id }
        rep_counts = mutation_counts['repertoire'][rep_id]
        for r in pos_names:
            rep_counts[r] = 0
        for r in region_names:
            rep_counts[r] = 0
        for r in count_names:
            rep_counts[r] = 0

    # add to repertoire
    add_mutation_count(rep_counts, row)
    return

def add_clone_count(rep_id, row):
    # need clone_d
    clone_id = row['clone_id']
    if clone_id is None:
        return

    # initialize counts per clone
    clone_counts = mutation_counts['clone'].get(rep_id)
    if clone_counts is None:
        mutation_counts['clone'][rep_id] = {}
    clone_counts = mutation_counts['clone'][rep_id].get(clone_id)
    if clone_counts is None:
        mutation_counts['clone'][rep_id][clone_id] = { 'repertoire_id':rep_id, 'clone_id':clone_id }
        clone_counts = mutation_counts['clone'][rep_id][clone_id]
        for r in pos_names:
            clone_counts[r] = 0
        for r in region_names:
            clone_counts[r] = 0
        for r in count_names:
            clone_counts[r] = 0

    # add to clone
    add_mutation_count(clone_counts, row)
    return

def unique_counts(counts):
    # count of unique mutated positions
    for row_id in counts:
        for i in range(1,105):
            f = 'mu_count_' + str(i) + '_r'
            if counts[row_id][f + '_aa'] > 0:
                counts[row_id]['mu_unique_r'] += 1
            f = 'mu_count_' + str(i) + '_s'
            if counts[row_id][f + '_aa'] > 0:
                counts[row_id]['mu_unique_s'] += 1

def compute_std(entry, repertoire_counters, field, N, groupDict, normalize):
    """Compute average and standard deviation"""
    group_total = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            if normalize:
                rep_value = float(repertoire_counters[rep_id][field]) / float(repertoire_counters[rep_id]['total_count'])
            else:
                rep_value = float(repertoire_counters[rep_id][field])
            group_total += rep_value
    if N > 0:
        group_avg = group_total / N
    else:
        group_avg = 0.0
    group_std = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            if normalize:
                rep_value = float(repertoire_counters[rep_id][field]) / float(repertoire_counters[rep_id]['total_count'])
            else:
                rep_value = float(repertoire_counters[rep_id][field])
            group_std += (group_avg - rep_value) * (group_avg - rep_value)
    if N > 1:
        group_std = math.sqrt(group_std / (N - 1.0))
    else:
        group_std = 0.0
    entry[field + '_avg'] = group_avg
    entry[field + '_std'] = group_std

def compute_group_usage(groupDict, counters, repertoire_counters):
    """Calculate group average and deviation for mutation counters"""
    N = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            N += 1.0
    counters['repertoire_group_id'] = groupDict['repertoire_group_id']
    counters['N'] = N
    for field in transfer_names:
        if field in total_names:
            compute_std(counters, repertoire_counters, field, N, groupDict, False)
        else:
            compute_std(counters, repertoire_counters, field, N, groupDict, True)

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # Mutation counts organized by repertoire
    for level in module_levels:
        mutation_counts[level] = {}

    if inputDict.get(defaults.groups_key) is not None:
        for group in inputDict[defaults.groups_key]:
            group_mutation_counts[group] = {}

def process_record(inputDict, metadataDict, currentFile, calc, row):
    """Perform calculation from given fields"""

    # get repertoire
    rep_id = row.get('repertoire_id')
    if rep_id is None:
        return
    if metadataDict.get(rep_id) is None:
        return

    # Mutation counts
    if countKey in calc['operations']:
        add_repertoire_count(rep_id, row)
        add_clone_count(rep_id, row)

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    # Mutation counts
    if countKey in calc['operations']:
        # generate unique counts
        unique_counts(mutation_counts['repertoire'])
        for rep_id in mutation_counts['clone']:
            unique_counts(mutation_counts['clone'][rep_id])

        # repertoire counts
        names = ['repertoire_id'] + transfer_names
        writer = csv.DictWriter(open('mutational_report.repertoire.csv', 'w'), fieldnames=names)
        writer.writeheader()
        for rep_id in mutation_counts['repertoire']: 
            writer.writerow(mutation_counts['repertoire'][rep_id])

        # clone counts
        names = ['repertoire_id', 'clone_id'] + transfer_names
        writer = csv.DictWriter(open('mutational_report.clone.csv', 'w'), fieldnames=names)
        writer.writeheader()
        for rep_id in mutation_counts['clone']: 
            for clone_id in mutation_counts['clone'][rep_id]: 
                writer.writerow(mutation_counts['clone'][rep_id][clone_id])

        # group counts
        if inputDict.get(defaults.groups_key) is not None:
            names = ['repertoire_group_id', 'N'] + group_names
            writer = csv.DictWriter(open('mutational_report.group.csv', 'w'), fieldnames=names)
            writer.writeheader()
            for group in inputDict[defaults.groups_key]:
                compute_group_usage(inputDict[defaults.groups_key][group], group_mutation_counts[group], mutation_counts['repertoire'])
                writer.writerow(group_mutation_counts[group])
