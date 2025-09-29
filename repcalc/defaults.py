"""
Default settings
"""

#
# defaults.py
# Default settings
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

# Info
#from repsum import __author__, __version__, __date__

import repcalc.gldb as gldb

# Keys in specification file
metadata_file_key = "metadata"
groups_file_key = "groups"
germline_file_key = "germline"
germline_key = "gldb"
repertoire_key = "repertoires"
full_metadata_key = "full_metadata"
rearrangement_files_key = "rearrangement_files"
processing_stage_key = "processing_stage"
organismKey = "organism"
filesKey = "files"
fileMetadataKey = "fileMetadata"
groups_key = "repertoire_groups"
summaryKey = "summary"
vdjmlKey = "vdjml"
changeoKey = "changeo"
calculationsKey = "calculations"
calcTypeKey = "type"
calcOpsKey = "operations"
calcFilters = "filters"

# Calculation modules
calculationModules = {
    "filters": {"filename":"filters", "require_germline":True},
    "gene_segment": {"filename":"gene_segment", "require_germline":True},
    "CDR3": {"filename":"cdr3", "require_germline":True},
    "diversity": {"filename":"diversity"},
    "mutations": {"filename":"mutations"},
    "clonal_assignment": {"filename":"tcr_clone", "require_germline":True},
    "lineage": {"filename":"lineage"}
}

# Header column names for summary file
headerNames = {
    "V_CALL": "V gene",
    "D_CALL": "D gene",
    "J_CALL": "J gene",
    "DUPCOUNT": "DUPCOUNT",
    "NonFuncCondition1": "Out-of-frame junction",
    "NonFuncCondition2": "Missing CYS",
    "NonFuncCondition3": "Missing TRP/PHE",
    "NonFuncCondition4": "Stop Codon?",
    "NonFuncCondition5": "Indels Found",
    "NonFuncCondition6": "Only Frame-Preserving Indels Found",
    "CDR3_AA": "CDR3 AA (imgt)",
    "CDR3_SEQ": "CDR3 NA (imgt)",
}



def get_dupcount(headerMapping, fields):
    dupcount_field = headerMapping.get(headerNames['DUPCOUNT'])
    if dupcount_field is None: return 1
    else: return int(fields[dupcount_field])

def get_duplicate_count(fields):
    if fields.get('duplicate_count') is None: return 1
    else: return int(fields['duplicate_count'])

def has_rearrangement_filter(inputDict, group):
    g = inputDict[groups_key][group]
    if g.get('filter') and g.get('filter').get('Rearrangement'):
        return True
    return False

def apply_filter(inputDict, group, fields):
    valid = True

    # helper function
    def logical_op(op, field_value, check_value):
        if op == '=':
            return field_value == check_value
        elif op == '!=':
            return field_value != check_value
        elif op == '<':
            return field_value < check_value
        elif op == '<=':
            return field_value <= check_value
        elif op == '>':
            return field_value > check_value
        elif op == '>=':
            return field_value >= check_value
        elif op == 'contains':
            return check_value in field_value
        else:
            return False

    # helper function
    def eval_exp(fields, exp):
        fname = exp['content']['field']

        # subgroup and gene are derived fields
        if fname in ['v_subgroup', 'd_subgroup', 'j_subgroup', 'c_subgroup']:
            fcall = fname.replace('_subgroup', '_call')
            value = fields.get(fcall)
            germline = inputDict[germline_key]
            ad = germline['allele_descriptions'].get(value)
            if ad:
                value = gldb.getSubgroup(ad)
            else:
                value = None
        elif fname in ['v_gene', 'd_gene', 'j_gene', 'c_gene']:
            fcall = fname.replace('_gene', '_call')
            value = fields.get(fcall)
            germline = inputDict[germline_key]
            ad = germline['allele_descriptions'].get(value)
            if ad:
                value = gldb.getGene(ad)
            else:
                value = None
        else:
            value = fields.get(exp['content']['field'])

        if not value:
            return False
        return logical_op(exp['op'], value, exp['content']['value'])

    # does group have a rearrangement filter
    g = inputDict[groups_key][group]
    if g.get('filter') and g.get('filter').get('Rearrangement'):
        f = g['filter']['Rearrangement']

        if f['op'].upper() == 'AND':
            valid = eval_exp(fields, f['content'][0]) and eval_exp(fields, f['content'][1]);
        elif f['op'].upper() == 'OR':
            valid = eval_exp(fields, f['content'][0]) or eval_exp(fields, f['content'][1]);
        else:
            valid = eval_exp(fields, f);

    return valid

