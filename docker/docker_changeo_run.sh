#
# Simple script to run Changeo within the docker image
#

# Change these:

# your input file and IgBlast results file
INPUT=/data/s166813/Downloads/18.fasta
IGBLAST_FILE=18.igblast.out

organism="human"
#organism="Mus_musculus"

seqType="IG"
#seqType="TR"

domainSystem="imgt"

# build up the arguments, generally do not need to change this
#
ARGS="-s $INPUT"
ARGS="$ARGS -i $IGBLAST_FILE"
ARGS="$ARGS -r $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${seqType}_VDJ.fna"
ARGS="$ARGS --regions --scores"

# run changeo
MakeDb.py igblast $ARGS
