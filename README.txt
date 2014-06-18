###########################################################
# README FOR repertoire-summarization SUITE of VDJ SERVER #
###########################################################

The purpose of this suite of software tools is to generate 
summary statistics of REP-SEQ data.

IGBLAST output files (using a certain output format) are 
taken as input and 4 distinct files are generated as outputs:
1) CDR3 length histogram (in kabat and imgt modes),
2) JSON-formatted hierachies of the IMGT alleles and 
corresponding counts of segments in top combinations
with agregrated counts, 3) read-level summary statistics (TSV)
, and 4 ) vdjml files of the parsed data.

The rep_char.py script is designed to support the above 4 functions.

If run over several igblast outputs that should subsequently be 
merged the following scripts can be used to logicaly merge outputs:
1) cdr3_hist.py can logically merge CDR3 length histogram data,
2) json_seg_merge.py can be used to logically merge JSON count data,
3) vdjml_merge.py can be used to merge VDJML files, and 4) merge_tables.py
can be used to merge TSV tables.  NOTE : For 1 and #2, since counts 
are being merged, input order is not important.  For #3 and #4, since 
files are being merged and read order too, order is important.  For
#3 and #4 the input order is used in the output.

To get help on any of the *.pyscripts run it as seen here, and 
use the "-h" flag for "help"!  For the .sh scripts please open them and read the comments.

esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ ./rep_char.py -h
usage: rep_char.py [-h] [-igblast_version IGBLAST_VERSION]
                   [-igblast_params IGBLAST_PARAMS]
                   [-igblast_runid IGBLAST_RUNID] [-igblast_uri IGBLAST_URI]
                   [-db_name_igblast DB_NAME_IGBLAST] [-db_ver DB_VER]
                   [-db_species_vdjml DB_SPECIES_VDJML] [-db_uri DB_URI]
                   [-json_out JSON_OUT] [-cdr3_hist_out CDR3_HIST_OUT]
                   [-skip_char]
                   igblast_in vdjml_out char_out vdj_db_root qry_fasta
                   {human,Mus_musculus}

Parse IGBLAST output, write a VDJML output file. Generate read-level
repertoire-characterization data.

positional arguments:
  igblast_in            path to igblast analysis file of results, hits,
                        segments, etc. *** CRITICAL NOTE *** SHOULD BE RUN
                        WITH OUTFMT AS SEEN HERE -outfmt '7 qseqid qgi qacc
                        qaccver qlen sseqid sallseqid sgi sallgi sacc saccver
                        sallacc slen qstart qend sstart send qseq sseq evalue
                        bitscore score length pident nident mismatch positive
                        gapopen gaps ppos frames qframe sframe btop'
  vdjml_out             the path to the output vdjml XML
  char_out              the path to the output TSV file of read-level
                        repertoire characterization data
  vdj_db_root           path to the VDJ directory root REQUIRED
  qry_fasta             path to the input fasta file of query (Rep-Seq) data
                        input to IgBLAST
  {human,Mus_musculus}  the organism IgBLASTed against;must exist under
                        vdj_db_root

optional arguments:
  -h, --help            show this help message and exit
  -igblast_version IGBLAST_VERSION
                        the version of igblast used placed in the VDJML output
  -igblast_params IGBLAST_PARAMS
                        the parameters passed to igblast placed in the VDJML
                        output
  -igblast_runid IGBLAST_RUNID
                        the ID assigned to the run placed in the VDJML output
  -igblast_uri IGBLAST_URI
                        a URI string for the IgBLAST aligner to be placed in
                        the VDJML output
  -db_name_igblast DB_NAME_IGBLAST
                        the name of the IgBLAST database placed in the VDJML
                        output
  -db_ver DB_VER        a version string associated with the database placed
                        in the VDJML output
  -db_species_vdjml DB_SPECIES_VDJML
                        db organism for VDJML output file(s) placed in the
                        VDJML output
  -db_uri DB_URI        a URI associated with the database placed in the VDJML
                        output
  -json_out JSON_OUT    output file for the JSON segment count IMGT hierarchy
  -cdr3_hist_out CDR3_HIST_OUT
                        output file for the CDR3 histogram of lengths (both
                        kabat and imgt systems)
  -skip_char            If this is set, then characterization besides
                        asignment and CDR3 length is skipped




#################################################################
# REPERTOIRE-CHARACTERIZATION REQUIRES THE VDJML PYTHON MODULE. #
#################################################################

NOTE : use of rep_char.py require that the python module "vdjml" be usable. Moreover
version 0.0.12 is required.

NOTE : rep_char.py expects to be able to find the file "vdjml_igblast_parse.py". 
Successful finding may be achieved by using the PYTHONPATH.  See the example here : 

esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ ./rep_char.py -h
Traceback (most recent call last):
  File "./rep_char.py", line 4, in <module>
    from vdjml_igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs

esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ ls -alh /home/data/vdj_server/vdjml/python/vdjml_igblast_parse.py 
-rwxrwxr-x 1 esalina2 esalina2 29K Jun 13 09:22 /home/data/vdj_server/vdjml/python/vdjml_igblast_parse.py

esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ export PYTHONPATH=/home/data/vdj_server/vdjml/python/


#####################################################################
# REPERTOIRE-CHARACTERIZATION REQUIRES A VDJ_DB DIRECTORY STRUCTURE #
#####################################################################

Required by rep_char.py (whose outputs are 1) VDJML files, 
2) cdr3 length histogram(s), 3) TSV output files with huge 
number of columns and 4) JSON-formatted segment counts is 
the VDJ_DB_ROOT which is a directory with files and directories 
placed and expected to be in an organized manner.

The program vdj_db_dl.py may be used to download required data/reference
files from imgt and set up the VDJ_DB directory and structures.

In brief summary, the files/directories required as well as their
uses are listed below:

*  PER organism a directory is expected under the root (you can see 
Mus_musculus and human in the listing below)
*  UNDER each organism is expected a directory "ReferenceDirectorySet" 
is expected (containing FASTA files and auto-downloaded HTML/fasta 
files and BLAST databases and domain-system region lookup data for 
region start/stop data for IMGT and KABAT)
*  UNDER each organism is expected a directory "GeneTables" 
(containing auto-downloaded and parsed HTML files of gene table 
[hierarchical] information which is parsed for creating .JSON 
output files for the chart viewing for knowing the hierarchy 
and how to zoom)
*  UNDER the root is www.imgt.dat which contains under it the 
GENE-DB directory and LIGM-DB directory with files downloaded
from imgt.org.  Files in those contain the reference IG and
TCR reference sequence data, annotation data, as well as other 
data.  See http://www.imgt.org/download/GENE-DB/README.txt and
also http://www.imgt.org/download/LIGM-DB/README for further 
information
*  UNDER each ReferenceDirectorySet directory are the IgBLAST-formatted
databases as well as fasta files (of 17 loci) of Ig and TCR IMGT data.
Also there are the IMGT and KABAT directories having region-delineation
indices for the reference sequences (e.g. CDR1 start/stop positions).

A "tree -d" listing is seen here showing an example :

esalina2@eddiecomp:/home/data/DATABASE/06_05_2014$ tree -d `pwd`
/home/data/DATABASE/06_05_2014
├── human
│   ├── GeneTables
│   └── ReferenceDirectorySet
│       ├── IMGT
│       └── KABAT
├── Mus_musculus
│   ├── GeneTables
│   └── ReferenceDirectorySet
│       ├── IMGT
│       └── KABAT
└── www.imgt.org
    └── download
        ├── GENE-DB
        └── LIGM-DB

14 directories





##############################################
# REPERTOIRE-CHARACTERIZATION TEST DIRECTORY #
##############################################

Under the Test_Data directory lies a driver script and data to test 
some core functionalities.

The driver script is tests.sh.  Besides testing, it can be viewed
to see some example calls/usages!

Test data include : human.IG.fna, human.TR.fna (human IG and TR data)
and Mus_musculus_IG.fna and  Mus_musculus_TR.fna (mouse IG and TR data).
During the course of script execution IgBLAST is invoked on the data
in both IMGT and KABAT domain modes.  Also, during the course of
execution, the vdj_sim.py script uses the VDJ_DB files to 
simulate VDJ recombination and somatic hypermutatation and
these data files are also used as part of the tests.

Tests generate VDJML files, characterization output files 
(TSV-formatted), CDR3 length histograms, and JSON segment count
data.  Output/error log files are also generated.



-Eddie Salinas
VDJ Server Team Member
