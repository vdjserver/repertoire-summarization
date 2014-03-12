#########################################################
#README FOR repertoire-summarization SUITE of VDJ SERVER#
#########################################################

The purpose of this suite of software tools is to generate 
summary statistics of REP-SEQ data.

IGBLAST output files (using a certain output format) are 
taken as input and 4 distinct files are generated as outputs:
1) CDR3 length histogram (in kabat and imgt modes),
2) JSON-formatted hierachies of the IMGT alleles and 
corresponding counts of segments in top combinations
with agregrated counts, 3) read-level summary statistics
, and 4 ) vdjml files of the parsed data.

The rep_char.py script is designed to support the above 4 functions.
If run over several igblast outputs that should be merged the
following script can be used to logicaly merge outputs:
[TODO INSERT LIST OF 4 MERGERS HERE]

To get help on any of the scripts run it as seen here, and 
use the "-h" flag for "help"!

-Eddie Salinas
VDJ Server Team Member








