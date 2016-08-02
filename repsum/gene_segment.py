"""
Gene segment calculation module
"""

# repsum modules
from .version import __version__
import utils
import defaults
import metadata

#class for count map
class IncrementMapWrapper():

	#have a count map
	count_map=None

	#have a default value
	default_value=0

	#constructor
	def __init__(self):
		self.count_map=dict()
		self.default_value=0

	#treating the keys as ints, find the max key
	#return -1 if no keys
	def getMaxKeyWithIntOrdering(self):
		keys=self.count_map.keys()
		int_keys=list()
		for k in keys:
			int_keys.append(int(k))
		int_keys.sort()
		if(len(int_keys)>0):
			return max(int_keys)
		else:
			return (-1)

	#be able to change the default value
	def setDefault(self,d):
		self.default_value=d

	#have a subroutine to return the map
	def get_map(self):
		return  self.count_map


	#have a getter method for the default value
	def get_val(self,key):
		if(key in self.count_map):
			return self.count_map[key]
		else:
			return self.default_value
	
	#increment the given value
	def increment(self,val):
		if(val in self.count_map):
			self.count_map[val]+=1
		else:
			self.count_map[val]=1

	#print it
	def printMap(self):
		print "I am an increment map."
		for k in self.count_map:
			print str(k)+"\t"+str(self.count_map[k])

	#zero out at a key
	def zeroOut(self,key):
		self.count_map[key]=0

	#given another increment map wrapper, merge its counts into this/self 
	def mergeInto(self,other):
		other_map=other.get_map()
		for k in other_map:
			k_count=other_map[k]
			#print "GOT k=",k," and k_count=",k_count," from other map"
			for i in range(k_count):
				self.increment(k)

	#write JSON to file
	def JSONIFYToFile(self,db_base,organism,filePath,filterbyFastaAlleles=False,pickle_file_full_path=None):
		JSON=self.JSONIFYIntoHierarchy(db_base,organism,filterbyFastaAlleles,pickle_file_full_path)
		writer=open(filePath,'w')
		writer.write(JSON)
		writer.close()


	#JSONIFY into hierarchy
	def JSONIFYIntoHierarchy(self,db_base,organism,filterbyFastaAlleles=False,pickle_file_full_path=None):
		hierarchy=getHierarchyBy(db_base+"/"+organism+"/GeneTables/",organism,filterbyFastaAlleles,pickle_file_full_path)
		JSON=jsonify_hierarchy(hierarchy,organism,self.count_map,"value")
		return JSON



segment_counters = {}

def initialize_calculation_module(inputDict, metadataDict, headerMapping):
    """Perform any module initialization"""
    # TODO: should check for required columns

    # setup data structures for groups
    groups = inputDict[defaults.groupsKey]
    for group in groups: segment_counters[group] = IncrementMapWrapper()

    return

def process_record(inputDict, metadataDict, headerMapping, groupSet, calc, fields):
    """Perform calculation from given fields"""
    for group in groupSet:
        gene = fields[headerMapping[defaults.headerNames['V_CALL']]]
        if gene is not None: segment_counters[group].increment(gene)
        gene = fields[headerMapping[defaults.headerNames['D_CALL']]]
        if gene is not None: segment_counters[group].increment(gene)
        gene = fields[headerMapping[defaults.headerNames['J_CALL']]]
        if gene is not None: segment_counters[group].increment(gene)

def finalize_calculation_module(inputDict, metadataDict):
    """Finalize and save the calculations"""
    groups = inputDict[defaults.groupsKey]
    for group in groups:
        print("group: " + group)
        segment_counters[group].printMap()
