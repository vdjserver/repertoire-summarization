"""
Default settings
"""

# Info
#from repsum import __author__, __version__, __date__

# Keys in specification file
metadata_file_key = "metadata"
filesKey = "files"
fileMetadataKey = "fileMetadata"
groupsKey = "groups"
summaryKey = "summary"
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
}



