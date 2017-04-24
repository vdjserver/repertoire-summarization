#
# Simple script to run RepSum within the docker image
#
# It has a long command line with numerous parameters,
# so modify script to suite your needs.
#

# Change these:

# your input file and IgBlast results file
INPUT=/data/Projects/research/immune/data/igblast-app-test/18.fasta
IGBLAST_FILE=igblast.out
VDJML_FILE=out.vdjml

organism="human"
#organism="Mus_musculus"

# build up the arguments, generally do not need to change this
RSARGS=""
RSARGS="$RSARGS $VDJML_FILE"
RSARGS="$RSARGS out.rc_out.tsv"
RSARGS="$RSARGS $VDJ_DB_ROOT" 
RSARGS="$RSARGS $INPUT" 
RSARGS="$RSARGS $organism" 

# run rep_char
repsum $RSARGS
