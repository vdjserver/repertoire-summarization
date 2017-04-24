"""
Default settings
"""

# Info
#from repsum import __author__, __version__, __date__

# Keys in specification file
metadata_file_key = "metadata"
organismKey = "organism"
filesKey = "files"
fileMetadataKey = "fileMetadata"
groupsKey = "groups"
summaryKey = "summary"
vdjmlKey = "vdjml"
changeoKey = "changeo"
calculationsKey = "calculations"
calcTypeKey = "type"
calcOpsKey = "operations"
calcFilters = "filters"

# Calculation modules
calculationModules = {
    "gene segment usage": {"filename":"gene_segment"},
    "CDR3": {"filename":"cdr3"},
    "diversity": {"filename":"diversity"},
    "mutations": {"filename":"mutations"},
    "clonality": {"filename":"clones"},
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
