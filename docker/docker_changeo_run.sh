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

# parse IgBlast output
MakeDb.py igblast $ARGS

# filter out non-functional sequences
fileBasename="${IGBLAST_FILE%.*}"
ParseDb.py split -d ${fileBasename}_db-pass.tab -f FUNCTIONAL

# assigning clones
DefineClones.py bygroup -d ${fileBasename}_db-pass_FUNCTIONAL-T.tab --act set --model hs5f --sym min --norm len --dist 0.165

#DefineClones.py bygroup -d ${fileBasename}_db-pass_FUNCTIONAL-T.tab --act set --model m1n --sym min --norm len --dist 0.165
