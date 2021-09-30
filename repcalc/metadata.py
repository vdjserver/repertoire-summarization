"""
Metadata functions
"""

#
# metadata.py
# Metadata utility functions
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

import sys

# repsum modules
#import defaults

def filename(file):
    """Extract filename from metadata record"""
    if (file):
        # trim compression extension
        name = file['value']['name']
        sname = name.split('.')
        if (sname[-1] == 'zip'): del sname[-1]
        if (sname[-1] == 'gz'): del sname[-1]
        if (sname[-1] == 'bz2'): del sname[-1]
        name = '.'.join(sname)
        return name
    else:
        return None
