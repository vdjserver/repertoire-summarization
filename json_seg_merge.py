#!/usr/bin/env python


import argparse
import argparse
from segment_utils import IncrementMapWrapper
from utils import readFileIntoString,printMap
import json
from imgt_utils import imgt_db


#assuming the same organism get a key_value map
#with terminal leaves and their values
def get_json_terminals_and_counts(j_data,existing_map):
	if(existing_map is None):
		existing_map=dict()
	if('children' in j_data):
		for child in j_data['children']:
			existing_map=get_json_terminals_and_counts(child,existing_map)
	elif(('label' in j_data) and ('value' in j_data)):
		key=j_data['label']
		value=int(j_data['value'])
		existing_map[key]=value
		return existing_map
	else:
		raise Exception("WARNING, either no children or no label/value found in ",j_data,"!")
	return existing_map




#use the first label found
def readOrganismFromJSON(j):
	#print "got a json"
	if('label' in j):
		#print "label fouind"
		return j['label']
	else:
		#print "label not found"
		return None
	

#program to merge JSON hierachy counts
if (__name__=="__main__"):
	parser=argparse.ArgumentParser(description='Merge multiple JSON hierachy counts into a single JSON.  Write the merged result to stdout')
	parser.add_argument('db_base_dir',type=str,nargs=1,help="path to the root of the VDJ server database")
	parser.add_argument('json_hier_in',type=str,nargs='+',help="path(s) to JSON hierachies to merge.  At least one is requried!  NOTE : all must be of one organism and of the same organism!  Allowed organisms are only those in the indicated VDJ database")
	args=parser.parse_args()
	if(args):
		#first verify that they're all the same organism
		#and if they're an allowed organism
		imgtdb_obj=imgt_db(args.db_base_dir[0])
		allowed_organisms=imgtdb_obj.getOrganismList()
		jsons_to_read=args.json_hier_in
		org_dict=dict()
		last_organism=None
		for j in jsons_to_read:
			json_string=readFileIntoString(j)
			json_obj=json.loads(json_string)
			j_organism=readOrganismFromJSON(json_obj)
			org_dict[j]=j_organism
			if(j_organism==None):
				raise Exception("Error, failure in determining organism from file "+j)
			elif(j_organism not in allowed_organisms):
				raise Exception("Error, organism ("+j_organism+") from file "+j)
			last_organism=j_organism
		for j in org_dict:
			if(org_dict[j]!=last_organism):
				raise Exception("Error, organism mismatch at file ",j," does not match ",last_organism,"!")
		#print "Using organism",last_organism
		imw=IncrementMapWrapper()
		for j in jsons_to_read:
			json_string=readFileIntoString(j)
			json_obj=json.loads(json_string)
			kvMap=get_json_terminals_and_counts(json_obj,None)
			for k in kvMap:
				val=int(kvMap[k])
				for i in range(val):
					#print "incrementing for value "
					imw.increment(k)
		imw.JSONIFYToFile(imgtdb_obj.getBaseDir(),last_organism,"/dev/stdout")			
	else:
		parser.print_help()


