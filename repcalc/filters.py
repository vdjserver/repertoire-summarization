"""
Filter rearrangements
"""

#
# filters.py
# Filter rearrangements
#
# VDJServer Analysis Portal
# Repertoire calculations and comparison
# https://vdjserver.org
#
# Copyright (C) 2025 The University of Texas Southwestern Medical Center
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
import repcalc.gldb as gldb

import json
import math
import numpy
import csv
import airr

# module operations

file_writers = {}

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    if inputDict.get(defaults.groups_file_key) is None:
        print('ERROR: No repertoire groups are defined.', file=sys.stderr)
        sys.exit(1)

    file_writers = {}

def process_record(inputDict, metadataDict, currentFile, calc, fields):
    """Perform calculation from given fields"""

    # get repertoire
    rep_id = fields['repertoire_id']
    if metadataDict.get(rep_id) is None:
        return

    # open output files
    groupList = metadata.groupsWithRepertoire(inputDict, rep_id)
    if file_writers.get(currentFile) is None:
        file_writers[currentFile] = {}
        for group in groupList:
            filename = currentFile.replace('.airr.tsv', '.filters.' + group + '.airr.tsv')
            if filename.endswith(".gz"):
                filename = filename.replace('.gz','')
            print('Output file: ' + filename)
            file_writers[currentFile][group] = airr.derive_rearrangement(filename, currentFile)

    # write data
    for group in groupList:
        if defaults.apply_filter(inputDict, group, fields):
            file_writers[currentFile][group].write(fields)

def finalize_calculation_module(inputDict, metadataDict, outputSpec, calc):
    """Finalize and save the calculations"""

