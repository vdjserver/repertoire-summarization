###########################################################
# README FOR repertoire-summarization SUITE of VDJ SERVER #
###########################################################

The purpose of this suite of software tools is to generate 
analysis and summary statistics of Ig-Seq data.

IGBLAST output files (using a certain output format) are 
taken as input and 4 distinct files are generated as outputs:
1) CDR3 length histogram (in kabat and imgt modes),
2) JSON-formatted hierachies of the IMGT alleles and 
   corresponding counts of segments in top combinations
with agregrated counts, 
3) read-level summary statistics (TSV) , and 
4) vdjml files of the parsed data.

The rep_char.py script is designed to support the above 4 functions.

If run over several igblast outputs that should subsequently be 
merged the following scripts can be used to logicaly merge outputs:
1) cdr3_hist.py can logically merge CDR3 length histogram data,
2) json_seg_merge.py can be used to logically merge JSON count data,
3) vdjml_merge.py can be used to merge VDJML files, and 
4) merge_tables.py
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



##################################################################
# REPERTOIRE-CHARACTERIZATION REQUIRES ADDITIONAL PYTHON MODULES #
##################################################################


TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE
INSERT NOTES ABOUT REQUIRED MODULES HERE



from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from bs4 import BeautifulSoup


from collections import defaultdict
from datetime import datetime


from os.path import basename
from pprint import pprint
from sets import Set
from subprocess import call
from subprocess import Popen
import argparse
import glob
import json
import math
import multiprocessing
import ntpath
import numpy
import os
import pickle
import pprint
import pprint, pickle
import Queue
import random
import re
import subprocess
import sys
import sys, traceback
import threading
import time
import traceback,sys
import urllib2
import yaml



#####################################################################
# REPERTOIRE-CHARACTERIZATION REQUIRES A VDJ_DB DIRECTORY STRUCTURE #
#####################################################################

Required by rep_char.py is the VDJ_DB_ROOT which is a directory with
files and directories downloaded/created and expected to be in an organized manner.

The program vdj_db_dl.py may be used to download required data/reference
files from www.imgt.org and set up the VDJ_DB directory and structures.
This program uses the downloaded data to generate other required files
in the VDJ_DB_ROOT.

The expected structure of the VDJ_DB_ROOT is described below.

*  UNDER the root is www.imgt.dat which contains under it the 
GENE-DB directory and LIGM-DB directory with files downloaded
from imgt.org.  Files in those contain the reference IG and
TCR reference sequence data, annotation data, as well as other 
data.  See http://www.imgt.org/download/GENE-DB/README.txt and
also http://www.imgt.org/download/LIGM-DB/README for further 
information.  

*  UNDER the www.imgt.org/download/GENE-DB directory is a file of
importance named IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.
It is from this file that the reference sequences are extracted.
Similarly, under the www.imgt.org/download/LIGM-DB directory is the
file imgt.dat.  It is from this file that annotation information
(such as J-segment IMGT-CDR3 end information) is extracted.

*  PER each organism a corresponding directory is expected under 
the root (you can see Mus_musculus and human in the listing below)

*  UNDER each organism directory is expected a directory "ReferenceDirectorySet" 

*  UNDER each "ReferenceDirectorySet" directory there are expected to be
17 fasta files of the form "LOCUS.fna" where "LOCUS" is ranges
over the 17 TR/IG loci : IGHV,IGKV,IGLV,IGHD,IGHJ,IGKJ,IGLJ,
TRAV,TRAJ,TRBV,TRBD,TRBJ,TRDV,TRDD,TRDJ,TRGV,TRGJ.  These files
can be created by the vdj_db_dl script and their contents represent
the fasta reference data from the loci for the organism.  The contents
come from IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP.
Extractions are made by matching on the organism, the locus, and also
the region type "V-REGION", "D-REGION", and "J-REGION".  NOTE
that the "Mus_musculus" directory contains "Mus musculus" and several
"Mus spretus" sequences ; this is because organism matching here is
based on "Mus".

*  UNDER each ReferenceDirectorySet directory are FNA files that have
any IMGT-gaps (".") removed and grouped by V, D, and J and are named
as ORGANISM_SEQTYPE_SEGMENT.fna. This way there are two V files
(for IG and TR), two D files, and two J files.  For these files are 
IgBLAST-formatted.  These are the primary databases used
by IgBLAST/VDJServer for germline inference.

*  UNDER each "ReferenceDirectorySet" directory are IMGT and KABAT
directories that contain lookup data for region start/stop data
for the corresponding schemes.  The lookup files are Vlooup.tsv
and Jlookup.tsv for the V and J segments respectively.  They are
tab-separated values files.  In the Vlookup.tsv file the allele name 
is in the first column and then the start/stop (1-based; not 0-based)
indices for the region FR1, CDR1, FR2, CDR2, FR3, and CDR3 regions.  
A value of -1 should appear when the start/stop value is not defined 
in the sequence (for example for CDR3 end for regions that are
not in the sequence ; for example shortened sequences not covering 
FR1 the first region).  For the Jlookup.tsv file only the allele
name is in the first column and the CDR3 end in the other. BOTH
files have 1-based indices that are in the space of the reference
sequences in the FNA files in the "ReferenceDirectorySet" directory.

*  UNDER each organism is expected a directory "GeneTables" 
(containing auto-downloaded and parsed HTML files of gene table 
[hierarchical] information.  The files are of the form LOCUS.html
and LOCUS.html.orphons.html.  The two files are downloaded from
imgt.org using the URLs 
http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=SPECIES&group=LOCUS
http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=SPECIES&group=LOCUS&orphon
respectively.  This way 34 html files are downloaded altogether.
These files are parsed and used for creating JSON-formatted 
hierarchies (using LOCI, subgroups, genes, and alleles) 
of segment counts.

*  OPTIONALLY ACCOMPANYING the LOCUS.html file sare LOCUS.html.patch files.  
The use and format of these files are described in the vdj_db_dl readme section 
in further detail.

A "tree -d" listing is seen here showing an example VDJ_DB and sub-directories:

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
