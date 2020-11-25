"""
Metadata functions
"""

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
