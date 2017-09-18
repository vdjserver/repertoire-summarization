"""
Metadata functions
"""

import sys

# repsum modules
import defaults

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
	else: return None

def filenames_from_uuids(metadataDict, uuidSet):
        """Return dictionary of file name for given set of uuids"""
        nameDict = {}
        for uuid in uuidSet: nameDict[uuid] = filename(metadataDict[defaults.fileMetadataKey][uuid])
        return nameDict

def sample_for_uuid(inputDict, uuid):
        groups = inputDict[defaults.groupsKey]
        for group in groups:
		if (groups[group]['type'] == 'sample'):
			for sample in groups[group]['samples']:
                                if sample == uuid: return group
        return None

# Determine if a sample contains a file
# sample/file group id can be passed in either parameter
def sample_contains_file(inputDict, group1, group2):
        if group1 == group2: return True
        if not inputDict[defaults.groupsKey].get(group1): return False
        if not inputDict[defaults.groupsKey].get(group2): return False
        group1_type = inputDict[defaults.groupsKey][group1]['type']
        if group1_type == 'sampleGroup': return False
        group2_type = inputDict[defaults.groupsKey][group2]['type']
        if group2_type == 'sampleGroup': return False
        if group1_type == 'sample' and group2_type == 'file':
                for sampleKey in inputDict[defaults.groupsKey][group1]['samples']:
                        for sampleFile in inputDict[defaults.groupsKey][group1]['samples'][sampleKey]:
                                if sampleFile == group2: return True
        elif group1_type == 'file' and group2_type == 'sample':
                for sampleKey in inputDict[defaults.groupsKey][group2]['samples']:
                        for sampleFile in inputDict[defaults.groupsKey][group2]['samples'][sampleKey]:
                                if sampleFile == group1: return True
        return False
