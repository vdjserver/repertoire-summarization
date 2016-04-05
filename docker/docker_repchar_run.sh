#
# Simple script to run rep_char within the docker image
#
# It has a long command line with numerous parameters,
# so modify script to suite your needs.
#

# Change these:

# your input file and IgBlast results file
INPUT=/data/s166813/Downloads/Test4.txt
IGBLAST_FILE=igblast.out

organism="human"
#organism="Mus_musculus"

# build up the arguments, generally do not need to change this
RCARGS=""
RCARGS="$RCARGS -combo_out combo.json" 
RCARGS="$RCARGS -json_out segment_counts.json"
RCARGS="$RCARGS -cdr3_hist_out cdr3_hist.tsv" 
RCARGS="$RCARGS -sample_json_out sample.json" 
RCARGS="$RCARGS $IGBLAST_FILE" 
RCARGS="$RCARGS out.vdjml"
RCARGS="$RCARGS out.rc_out.tsv"
RCARGS="$RCARGS $VDJ_DB_ROOT" 
RCARGS="$RCARGS $INPUT" 
RCARGS="$RCARGS $organism" 

# run rep_char
/repsum-root/rep_char.py $RCARGS
