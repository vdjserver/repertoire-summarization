"""
Metadata functions
"""

import sys

# repsum modules
import defaults

def file_with_uuid(metadataDict, uuid):
        """Return file metadata record for given uuid"""
	for file in metadataDict[defaults.filesKey]:
		if (file['uuid'] == uuid): return file
		
	return None

def filename(file):
        """Extract filename from metadata record"""
	if (file): return file['value']['name']
	else: return None

def filenames_from_uuids(metadataDict, uuidSet):
        """Return dictionary of file name for given set of uuids"""
        nameDict = {}
        for uuid in uuidSet:
                mfile = file_with_uuid(metadataDict, uuid)
                if (mfile): nameDict[uuid] = filename(mfile)
        return nameDict

