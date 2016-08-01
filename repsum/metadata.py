"""
Metadata functions
"""

import sys

# repsum modules
import defaults

def file_with_uuid(metadataDict, uuid):
	for file in metadataDict[defaults.files_key]:
		if (file['uuid'] == uuid): return file
		
	return None

def filename(file):
	if (file): return file['value']['name']
	else: return None

