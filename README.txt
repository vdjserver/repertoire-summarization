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
If run over several igblast outputs that should subsequently be 
merged the following scripts can be used to logicaly merge outputs:
1) cdr3_hist.py can logically merge CDR3 length histogram data,
2) json_seg_merge.py can be used to logically merge JSON count data,
3) vdjml_merge.py can be used to merge VDJML files, and 4) merge_tables.py
can be used to merge TSV tables.

To get help on any of the *.pyscripts run it as seen here, and 
use the "-h" flag for "help"!  For the .sh scripts please open them and read the comments.

NOTE : use of rep_char.py require that the python module "vdjml" be usable.

NOTE : rep_char.py expect to be able to find "vdjml_igblast_parse.py".  

Much testing was done with this version :
esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ md5sum ./vdjml_igblast_parse.py 
6f4571fad519527138be25cc62baab6a  ./vdjml_igblast_parse.py

One way to enable finding of the file is  with a soft-link as seen here with the link pointing to the file
esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ ls -alh vdjml_igblast_parse.py 
lrwxrwxrwx 1 esalina2 esalina2 31 Mar 25 10:05 vdjml_igblast_parse.py -> ../vdjml/vdjml_igblast_parse.py



-Eddie Salinas
VDJ Server Team Member








