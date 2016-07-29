"""
Default parameters
"""

# Info
#from repsum import __author__, __version__, __date__

# Keys in specification file
metadata_file_key = "metadata"
vdjml_file_key = "vdjml"
sequence_file_key = "sequence"
summary_file_key = "summary"
calculations_key = "calculations"
calc_type_key = "type"

# Calculation types
calculation_types = [
	"summary",
	"gene segment usage",
	"CDR3",
	"diversity",
	"mutations",
	"clonality",
	"lineage"]
