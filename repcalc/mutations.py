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
from Bio.Seq import Seq

# module operations
countKey = "count"
frequencyKey = "frequency"
module_levels = [ "rearrangement", "repertoire", "clone", "allele", "gene", "subgroup", "repertoire_group" ]

rearrangement_writers = {}

mutation_counts = {}
mutation_frequency = {}
group_mutation_counts = {}
group_mutation_frequency = {}

# aa alignments with gaps
alignment_names = ['cdr1_aa', 'cdr2_aa', 'fwr1_aa', 'fwr2_aa', 'fwr3_aa']
# counts of mutable positions in regions
region_total_names_nt = ['mu_total_count_fwr1', 'mu_total_count_cdr1', 'mu_total_count_fwr2', 'mu_total_count_cdr2', 'mu_total_count_fwr3']
region_total_names_aa = ['mu_total_count_fwr1_aa', 'mu_total_count_cdr1_aa', 'mu_total_count_fwr2_aa', 'mu_total_count_cdr2_aa', 'mu_total_count_fwr3_aa']
alignment_total_names = region_total_names_nt + region_total_names_aa

# global counts
total_names = ['mu_seq_count','mu_duplicate_count','mu_total_count','mu_total_count_aa','mu_unique_r','mu_unique_s']
count_names = total_names + ['mu_count_r_aa','mu_count_r','mu_count_s_aa','mu_count_s'] + alignment_total_names
freq_names_nt = ['mu_freq_r', 'mu_freq_s']
freq_names_aa = ['mu_freq_r_aa', 'mu_freq_s_aa']
freq_names = freq_names_nt + freq_names_aa

# mu_count_X_r/s - region/position counts for R/S mutations
# mu_total_count_X - normalization factor per position X
# mu_freq_* - mu_count_* / mu_total_count_*
region_names_nt = ['mu_count_fwr1_r','mu_count_fwr1_s','mu_count_cdr1_r','mu_count_cdr1_s','mu_count_fwr2_r','mu_count_fwr2_s','mu_count_cdr2_r','mu_count_cdr2_s','mu_count_fwr3_r','mu_count_fwr3_s']
region_names_aa = ['mu_count_fwr1_r_aa','mu_count_fwr1_s_aa','mu_count_cdr1_r_aa','mu_count_cdr1_s_aa','mu_count_fwr2_r_aa','mu_count_fwr2_s_aa','mu_count_cdr2_r_aa','mu_count_cdr2_s_aa','mu_count_fwr3_r_aa','mu_count_fwr3_s_aa']
region_names = region_names_nt + region_names_aa
aa_pos_names = []
nt_pos_names = []
pos_total_names = []
for i in range(1,105):
    f = 'mu_count_' + str(i) + '_r'
    aa_pos_names.append(f + '_aa')
    nt_pos_names.append(f)
    f = 'mu_count_' + str(i) + '_s'
    aa_pos_names.append(f + '_aa')
    nt_pos_names.append(f)
    pos_total_names.append('mu_total_count_' + str(i))
    pos_total_names.append('mu_total_count_' + str(i) + '_aa')
pos_names = pos_total_names + nt_pos_names + aa_pos_names

# region frequencies
region_freq_names_nt = ['mu_freq_fwr1_r','mu_freq_fwr1_s','mu_freq_cdr1_r','mu_freq_cdr1_s','mu_freq_fwr2_r','mu_freq_fwr2_s','mu_freq_cdr2_r','mu_freq_cdr2_s','mu_freq_fwr3_r','mu_freq_fwr3_s']
region_freq_names_aa = ['mu_freq_fwr1_r_aa','mu_freq_fwr1_s_aa','mu_freq_cdr1_r_aa','mu_freq_cdr1_s_aa','mu_freq_fwr2_r_aa','mu_freq_fwr2_s_aa','mu_freq_cdr2_r_aa','mu_freq_cdr2_s_aa','mu_freq_fwr3_r_aa','mu_freq_fwr3_s_aa']
region_freq_names = region_freq_names_nt + region_freq_names_aa

# position frequencies
pos_freq_names_nt = []
pos_freq_names_aa = []
for i in range(1,105):
    f = 'mu_freq_' + str(i) + '_r'
    pos_freq_names_aa.append(f + '_aa')
    pos_freq_names_nt.append(f)
    f = 'mu_freq_' + str(i) + '_s'
    pos_freq_names_aa.append(f + '_aa')
    pos_freq_names_nt.append(f)
pos_freq_names = pos_freq_names_nt + pos_freq_names_aa

# fields written for rearrangement level
row_transfer_names = alignment_names + count_names + region_names + aa_pos_names
# fields used in summary counters
transfer_names = count_names + region_names + pos_names
freq_transfer_names = freq_names + region_freq_names + pos_freq_names

# average and std for groups
group_names = []
for field in transfer_names:
    group_names.append(field + '_avg')
    group_names.append(field + '_std')
    group_names.append(field + '_N')
group_freq_names = []
for field in freq_transfer_names:
    group_freq_names.append(field + '_avg')
    group_freq_names.append(field + '_std')
    group_freq_names.append(field + '_N')

# count mutable nt and aa positions
# add aa translation with gaps
def add_alignment_fields_to_row(row):
    """Add amino acid alignment sequences and counts to row"""
    for n in alignment_names:
        nt_n = n.replace('_aa','')
        seq_nt = row[nt_n].replace('-','.')
        seq = seq_nt.replace('.','N')
        row['mu_seq_nt_' + n] = seq
        for i in range(0,len(seq)):
            if seq[i] != 'N':
                #row['mu_total_count'] += 1
                pass
        try:
            seq_aa = str(Seq(seq).translate())
        except Exception as e:
            print('Cannot translate field', nt_n, ':', seq)
            raise e
        row[n] = ''
        for i in range(0,len(seq_aa)):
            if seq_aa[i] == 'X':
                if seq_nt[i*3] == '.' and seq_nt[i*3+1] == '.' and seq_nt[i*3+2] == '.':
                    row[n] += '.'
                else:
                    row[n] += 'X'
            else:
                row[n] += seq_aa[i]
        len_name = 'mu_total_count_' + nt_n
        row[len_name] = 0
        for c in seq_nt:
            if c != '.':
                row[len_name] += 1
        len_name = 'mu_total_count_' + n
        row[len_name] = 0
        for c in row[n]:
            if c != '.':
                row[len_name] += 1

# nt mutation counts provided as input
# calculate aa mutation counts, total counts
def add_mutation_count_to_row(row):
    """Add mutation counts to row"""

    # multiplicity counts
    # TODO: do we need these, they seem redundant?
    row['mu_seq_count'] = 1
    row['mu_duplicate_count'] = defaults.get_duplicate_count(row)

    # position and total counts
    row['mu_count_r_aa'] = 0
    row['mu_count_r'] = 0
    row['mu_count_s_aa'] = 0
    row['mu_count_s'] = 0
    for i in range(1,105):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            row[f + '_aa'] = 1
            row['mu_count_r'] += int(row[f])
            row['mu_count_r_aa'] += 1
        else:
            row[f + '_aa'] = 0
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            row[f + '_aa'] = 1
            row['mu_count_s'] += int(row[f])
            row['mu_count_s_aa'] += 1
        else:
            row[f + '_aa'] = 0

    # region counts
    for f in region_names:
        row[f] = 0
    for i in range(1,27):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            row['mu_count_fwr1_r_aa'] += 1
            row['mu_count_fwr1_r'] += int(row[f])
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            row['mu_count_fwr1_s_aa'] += 1
            row['mu_count_fwr1_s'] += int(row[f])

    for i in range(27,39):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            row['mu_count_cdr1_r_aa'] += 1
            row['mu_count_cdr1_r'] += int(row[f])
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            row['mu_count_cdr1_s_aa'] += 1
            row['mu_count_cdr1_s'] += int(row[f])

    for i in range(39,56):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            row['mu_count_fwr2_r_aa'] += 1
            row['mu_count_fwr2_r'] += int(row[f])
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            row['mu_count_fwr2_s_aa'] += 1
            row['mu_count_fwr2_s'] += int(row[f])

    for i in range(56,66):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            row['mu_count_cdr2_r_aa'] += 1
            row['mu_count_cdr2_r'] += int(row[f])
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            row['mu_count_cdr2_s_aa'] += 1
            row['mu_count_cdr2_s'] += int(row[f])

    for i in range(66,105):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            row['mu_count_fwr3_r_aa'] += 1
            row['mu_count_fwr3_r'] += int(row[f])
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            row['mu_count_fwr3_s_aa'] += 1
            row['mu_count_fwr3_s'] += int(row[f])

# accumulate counts into summary counters with multiplicity
def add_mutation_count(counts, row):
    """Add mutation counts to counters"""
    duplicate_count = defaults.get_duplicate_count(row)
    # total counts
    counts['mu_seq_count'] += 1
    counts['mu_duplicate_count'] += duplicate_count

    # position counts
    for i in range(1,105):
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts[f + '_aa'] += duplicate_count
            counts['mu_count_r_aa'] += duplicate_count
            counts[f] += int(row[f]) * duplicate_count
            counts['mu_count_r'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts[f + '_aa'] += duplicate_count
            counts['mu_count_s_aa'] += duplicate_count
            counts[f] += int(row[f]) * duplicate_count
            counts['mu_count_s'] += int(row[f]) * duplicate_count

    # region and global counts
    for i in range(1,27):
        if row['fwr1_aa'][i - 1] != '.':
            # mutable position
            counts['mu_total_count'] += 3 * duplicate_count
            counts['mu_total_count_aa'] += duplicate_count
            counts['mu_total_count_fwr1'] += 3 * duplicate_count
            counts['mu_total_count_fwr1_aa'] += duplicate_count
            counts['mu_total_count_' + str(i)] += 3 * duplicate_count
            counts['mu_total_count_' + str(i) + '_aa'] += duplicate_count
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['mu_count_fwr1_r_aa'] += duplicate_count
            counts['mu_count_fwr1_r'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['mu_count_fwr1_s_aa'] += duplicate_count
            counts['mu_count_fwr1_s'] += int(row[f]) * duplicate_count

    for i in range(27,39):
        if row['cdr1_aa'][i - 27] != '.':
            # mutable position
            counts['mu_total_count'] += 3 * duplicate_count
            counts['mu_total_count_aa'] += duplicate_count
            counts['mu_total_count_cdr1'] += 3 * duplicate_count
            counts['mu_total_count_cdr1_aa'] += duplicate_count
            counts['mu_total_count_' + str(i)] += 3 * duplicate_count
            counts['mu_total_count_' + str(i) + '_aa'] += duplicate_count
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['mu_count_cdr1_r_aa'] += duplicate_count
            counts['mu_count_cdr1_r'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['mu_count_cdr1_s_aa'] += duplicate_count
            counts['mu_count_cdr1_s'] += int(row[f]) * duplicate_count

    for i in range(39,56):
        if row['fwr2_aa'][i - 39] != '.':
            # mutable position
            counts['mu_total_count'] += 3 * duplicate_count
            counts['mu_total_count_aa'] += duplicate_count
            counts['mu_total_count_fwr2'] += 3 * duplicate_count
            counts['mu_total_count_fwr2_aa'] += duplicate_count
            counts['mu_total_count_' + str(i)] += 3 * duplicate_count
            counts['mu_total_count_' + str(i) + '_aa'] += duplicate_count
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['mu_count_fwr2_r_aa'] += duplicate_count
            counts['mu_count_fwr2_r'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['mu_count_fwr2_s_aa'] += duplicate_count
            counts['mu_count_fwr2_s'] += int(row[f]) * duplicate_count

    for i in range(56,66):
        if row['cdr2_aa'][i - 56] != '.':
            # mutable position
            counts['mu_total_count'] += 3 * duplicate_count
            counts['mu_total_count_aa'] += duplicate_count
            counts['mu_total_count_cdr2'] += 3 * duplicate_count
            counts['mu_total_count_cdr2_aa'] += duplicate_count
            counts['mu_total_count_' + str(i)] += 3 * duplicate_count
            counts['mu_total_count_' + str(i) + '_aa'] += duplicate_count
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['mu_count_cdr2_r_aa'] += duplicate_count
            counts['mu_count_cdr2_r'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['mu_count_cdr2_s_aa'] += duplicate_count
            counts['mu_count_cdr2_s'] += int(row[f]) * duplicate_count

    for i in range(66,105):
        if row['fwr3_aa'][i - 66] != '.':
            # mutable position
            counts['mu_total_count'] += 3 * duplicate_count
            counts['mu_total_count_aa'] += duplicate_count
            counts['mu_total_count_fwr3'] += 3 * duplicate_count
            counts['mu_total_count_fwr3_aa'] += duplicate_count
            counts['mu_total_count_' + str(i)] += 3 * duplicate_count
            counts['mu_total_count_' + str(i) + '_aa'] += duplicate_count
        f = 'mu_count_' + str(i) + '_r'
        if int(row[f]) > 0:
            counts['mu_count_fwr3_r_aa'] += duplicate_count
            counts['mu_count_fwr3_r'] += int(row[f]) * duplicate_count
        f = 'mu_count_' + str(i) + '_s'
        if int(row[f]) > 0:
            counts['mu_count_fwr3_s_aa'] += duplicate_count
            counts['mu_count_fwr3_s'] += int(row[f]) * duplicate_count


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

def generate_frequency(row, count_row):
    # global frequencies
    tot_nt = 0
    tot_aa = 0
    for n in alignment_names:
        nt_n = n.replace('_aa','')
        tot_aa += float(count_row['mu_total_count_' + n])
        tot_nt += float(count_row['mu_total_count_' + nt_n])
    if tot_nt == 0:
        row['mu_freq_r'] = 0
        row['mu_freq_s'] = 0
    else:
        row['mu_freq_r'] = float(count_row['mu_count_r']) / tot_nt
        row['mu_freq_s'] = float(count_row['mu_count_s']) / tot_nt
    if tot_aa == 0:
        row['mu_freq_r_aa'] = 0
        row['mu_freq_s_aa'] = 0
    else:
        row['mu_freq_r_aa'] = float(count_row['mu_count_r_aa']) / tot_aa
        row['mu_freq_s_aa'] = float(count_row['mu_count_s_aa']) / tot_aa

    # region frequencies
    for n in alignment_names:
        nt_n = n.replace('_aa','')
        tot_aa = float(count_row['mu_total_count_' + n])
        if tot_aa == 0:
            row['mu_freq_' + nt_n + '_r_aa'] = 0
            row['mu_freq_' + nt_n + '_s_aa'] = 0
        else:
            row['mu_freq_' + nt_n + '_r_aa'] = float(count_row['mu_count_' + nt_n + '_r_aa']) / tot_aa
            row['mu_freq_' + nt_n + '_s_aa'] = float(count_row['mu_count_' + nt_n + '_s_aa']) / tot_aa
        tot_nt = float(count_row['mu_total_count_' + nt_n])
        if tot_nt == 0:
            row['mu_freq_' + nt_n + '_r'] = 0
            row['mu_freq_' + nt_n + '_s'] = 0
        else:
            row['mu_freq_' + nt_n + '_r'] = float(count_row['mu_count_' + nt_n + '_r']) / tot_nt
            row['mu_freq_' + nt_n + '_s'] = float(count_row['mu_count_' + nt_n + '_s']) / tot_nt

    # position frequencies
    for i in range(1,105):
        tot_cnt = float(count_row['mu_total_count_' + str(i)])
        if tot_cnt == 0:
            row['mu_freq_' + str(i) + '_r'] = 0
            row['mu_freq_' + str(i) + '_r_aa'] = 0
            row['mu_freq_' + str(i) + '_s'] = 0
            row['mu_freq_' + str(i) + '_s_aa'] = 0
        else:
            row['mu_freq_' + str(i) + '_r'] = float(count_row['mu_count_' + str(i) + '_r']) / tot_cnt
            row['mu_freq_' + str(i) + '_r_aa'] = float(count_row['mu_count_' + str(i) + '_r_aa']) / tot_cnt
            row['mu_freq_' + str(i) + '_s'] = float(count_row['mu_count_' + str(i) + '_s']) / tot_cnt
            row['mu_freq_' + str(i) + '_s_aa'] = float(count_row['mu_count_' + str(i) + '_s_aa']) / tot_cnt

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

def old_compute_std(entry, repertoire_counters, field, total_field, N, groupDict, normalize):
    """Compute average and standard deviation"""
    group_total = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_counters.get(rep_id) is not None:
            if normalize:
                rep_value = float(repertoire_counters[rep_id][field]) / float(repertoire_counters[rep_id]['mu_total_count'])
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
                rep_value = float(repertoire_counters[rep_id][field]) / float(repertoire_counters[rep_id]['mu_total_count'])
            else:
                rep_value = float(repertoire_counters[rep_id][field])
            group_std += (group_avg - rep_value) * (group_avg - rep_value)
    if N > 1:
        group_std = math.sqrt(group_std / (N - 1.0))
    else:
        group_std = 0.0
    entry[field + '_avg'] = group_avg
    entry[field + '_std'] = group_std

def compute_std(entry, repertoire_values, field, repertoire_counters, total_field, tot_N, groupDict, normalize):
    """Compute average and standard deviation"""
    group_total = 0.0
    N = tot_N
    tot_cnt = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_values.get(rep_id) is not None:
            if float(repertoire_counters[rep_id][total_field]) == 0:
                # repertoire does not have data
                N -= 1.0
                rep_value = 0.0
            else:
                rep_value = float(repertoire_values[rep_id][field])
            group_total += rep_value
    if N > 0:
        group_avg = group_total / N
    else:
        group_avg = 0.0
    group_std = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_values.get(rep_id) is not None:
            if float(repertoire_counters[rep_id][total_field]) == 0:
                # repertoire does not have data
                rep_value = 0.0
            else:
                rep_value = float(repertoire_values[rep_id][field])
                group_std += (group_avg - rep_value) * (group_avg - rep_value)
    if N > 1:
        group_std = math.sqrt(group_std / (N - 1.0))
    else:
        group_std = 0.0
    entry[field + '_avg'] = group_avg
    entry[field + '_std'] = group_std
    entry[field + '_N'] = N

#def compute_group_usage(groupDict, counters, repertoire_counters):
#    """Calculate group average and deviation for mutation counters"""
#    N = 0.0
#    for rep in groupDict['repertoires']:
#        rep_id = rep['repertoire_id']
#        if repertoire_counters.get(rep_id) is not None:
#            N += 1.0
#    counters['repertoire_group_id'] = groupDict['repertoire_group_id']
#    counters['N'] = N
#    for field in transfer_names:
#        if field in total_names:
#            compute_std(counters, repertoire_counters, field, N, groupDict, False)
#        else:
#            compute_std(counters, repertoire_counters, field, N, groupDict, True)

def compute_group_frequency(groupDict, freqs, repertoire_freqs, repertoire_counters):
    """Calculate group average and deviation for mutation frequencies"""
    tot_N = 0.0
    for rep in groupDict['repertoires']:
        rep_id = rep['repertoire_id']
        if repertoire_freqs.get(rep_id) is not None:
            tot_N += 1.0
    freqs['repertoire_group_id'] = groupDict['repertoire_group_id']
    freqs['N'] = tot_N

    # global
    for field in freq_names_nt:
        compute_std(freqs, repertoire_freqs, field, repertoire_counters, 'mu_total_count', tot_N, groupDict, False)
    for field in freq_names_aa:
        compute_std(freqs, repertoire_freqs, field, repertoire_counters, 'mu_total_count_aa', tot_N, groupDict, False)

    # regions
    for field in region_freq_names_nt:
        total_name = field.replace('_r','').replace('_s','').replace('mu_freq','mu_total_count')
        compute_std(freqs, repertoire_freqs, field, repertoire_counters, total_name, tot_N, groupDict, False)
    for field in region_freq_names_aa:
        total_name = field.replace('_r','').replace('_s','').replace('mu_freq','mu_total_count')
        compute_std(freqs, repertoire_freqs, field, repertoire_counters, total_name, tot_N, groupDict, False)

    # positions
    for field in pos_freq_names_nt:
        total_name = field.replace('_r','').replace('_s','').replace('mu_freq','mu_total_count')
        compute_std(freqs, repertoire_freqs, field, repertoire_counters, total_name, tot_N, groupDict, False)
    for field in pos_freq_names_aa:
        total_name = field.replace('_r','').replace('_s','').replace('_aa','').replace('mu_freq','mu_total_count')
        compute_std(freqs, repertoire_freqs, field, repertoire_counters, total_name, tot_N, groupDict, False)

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # Mutation counts organized by repertoire
    for level in module_levels:
        if level == 'rearrangement':
            continue
        mutation_counts[level] = {}
        mutation_frequency[level] = {}

    if inputDict.get(defaults.groups_key) is not None:
        for group in inputDict[defaults.groups_key]:
            group_mutation_counts[group] = {}
            group_mutation_frequency[group] = {}

def process_record(inputDict, metadataDict, currentFile, calc, row):
    """Perform calculation from given fields"""

    # get repertoire
    rep_id = row.get('repertoire_id')
    if rep_id is None:
        return
    if metadataDict.get(rep_id) is None:
        return

    # Rearrangement level
    if 'rearrangement' in calc['levels']:
        if rearrangement_writers.get(rep_id) is None:
            if inputDict.get(defaults.processing_stage_key) is not None:
                rearrangement_writers[rep_id] = airr.derive_rearrangement(rep_id + inputDict.get(defaults.processing_stage_key) + '.mutations.airr.tsv', currentFile, row_transfer_names)
            else:
                rearrangement_writers[rep_id] = airr.derive_rearrangement(rep_id + '.mutations.airr.tsv', currentFile, row_transfer_names)

    # compute additional fields
    add_alignment_fields_to_row(row)
    add_mutation_count_to_row(row)

    # Mutation counts
    if 'rearrangement' in calc['levels']:
        rearrangement_writers[rep_id].write(row)
    if 'repertoire' in calc['levels'] or 'repertoire_group' in calc['levels']:
        add_repertoire_count(rep_id, row)
    if 'clone' in calc['levels']:
        add_clone_count(rep_id, row)

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""
    # Mutation counts
    if countKey in calc['operations']:
        # generate unique counts
        #unique_counts(mutation_counts['repertoire'])
        #for rep_id in mutation_counts['clone']:
        #    unique_counts(mutation_counts['clone'][rep_id])

        # repertoire counts
        if 'repertoire' in calc['levels']:
            names = ['repertoire_id'] + transfer_names
            filename = 'repertoire.count.mutational_report.csv'
            if inputDict.get(defaults.processing_stage_key) is not None:
                filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
            writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
            writer.writeheader()
            for rep_id in mutation_counts['repertoire']: 
                writer.writerow(mutation_counts['repertoire'][rep_id])

        # clone counts
        if 'clone' in calc['levels']:
            names = ['repertoire_id', 'clone_id'] + transfer_names
            filename = 'clone.count.mutational_report.csv'
            if inputDict.get(defaults.processing_stage_key) is not None:
                filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
            writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
            writer.writeheader()
            for rep_id in mutation_counts['clone']: 
                for clone_id in mutation_counts['clone'][rep_id]: 
                    writer.writerow(mutation_counts['clone'][rep_id][clone_id])

        # repertoire counts by group
        if 'repertoire_group' in calc['levels']:
            if inputDict.get(defaults.groups_key) is not None:
                for group in inputDict[defaults.groups_key]:
                    names = ['repertoire_id'] + transfer_names
                    filename = 'repertoire.count.mutational_report.csv'
                    if inputDict.get(defaults.processing_stage_key) is not None:
                        filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
                    filename = group + '.' + filename
                    writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
                    writer.writeheader()
                    for rep in inputDict[defaults.groups_key][group]['repertoires']:
                        rep_id = rep['repertoire_id']
                        if mutation_counts['repertoire'].get(rep_id):
                            writer.writerow(mutation_counts['repertoire'][rep_id])

        # TODO: clone counts by group
        # TODO: other levels

    # Mutation frequencies
    if frequencyKey in calc['operations']:
        # repertoire frequencies
        if 'repertoire' in calc['levels']:
            names = ['repertoire_id'] + freq_transfer_names
            filename = 'repertoire.frequency.mutational_report.csv'
            if inputDict.get(defaults.processing_stage_key) is not None:
                filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
            writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
            writer.writeheader()
            for rep_id in mutation_counts['repertoire']: 
                mutation_frequency['repertoire'][rep_id] = { 'repertoire_id':rep_id }
                generate_frequency(mutation_frequency['repertoire'][rep_id], mutation_counts['repertoire'][rep_id])
                writer.writerow(mutation_frequency['repertoire'][rep_id])

        # clone frequencies
        if 'clone' in calc['levels']:
            names = ['repertoire_id', 'clone_id'] + freq_transfer_names
            filename = 'clone.frequency.mutational_report.csv'
            if inputDict.get(defaults.processing_stage_key) is not None:
                filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
            writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
            writer.writeheader()
            for rep_id in mutation_counts['clone']:
                for clone_id in mutation_counts['clone'][rep_id]:
                    clone_freq = { 'repertoire_id':rep_id, 'clone_id':clone_id }
                    generate_frequency(clone_freq, mutation_counts['clone'][rep_id][clone_id])
                    writer.writerow(clone_freq)

        # group frequencies
        if 'repertoire_group' in calc['levels']:
            if inputDict.get(defaults.groups_key) is not None:
                names = ['repertoire_group_id', 'N'] + group_freq_names
                filename = 'repertoire_group.frequency.mutational_report.csv'
                if inputDict.get(defaults.processing_stage_key) is not None:
                    filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
                writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
                writer.writeheader()
                for group in inputDict[defaults.groups_key]:
                    compute_group_frequency(inputDict[defaults.groups_key][group], group_mutation_frequency[group], mutation_frequency['repertoire'], mutation_counts['repertoire'])
                    writer.writerow(group_mutation_frequency[group])

        # repertoire frequencies by group
        if 'repertoire_group' in calc['levels']:
            if inputDict.get(defaults.groups_key) is not None:
                for group in inputDict[defaults.groups_key]:
                    names = ['repertoire_id'] + freq_transfer_names
                    filename = 'repertoire.frequency.mutational_report.csv'
                    if inputDict.get(defaults.processing_stage_key) is not None:
                        filename = inputDict.get(defaults.processing_stage_key) + '.' + filename
                    filename = group + '.' + filename
                    writer = csv.DictWriter(open(filename, 'w'), fieldnames=names)
                    writer.writeheader()
                    for rep in inputDict[defaults.groups_key][group]['repertoires']:
                        rep_id = rep['repertoire_id']
                        if mutation_counts['repertoire'].get(rep_id):
                            mutation_frequency['repertoire'][rep_id] = { 'repertoire_id':rep_id }
                            generate_frequency(mutation_frequency['repertoire'][rep_id], mutation_counts['repertoire'][rep_id])
                            writer.writerow(mutation_frequency['repertoire'][rep_id])

        # TODO: other levels
