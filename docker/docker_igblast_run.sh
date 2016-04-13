#
# Simple script to run IgBlast within the docker image
#
# It has a long command line with numerous parameters,
# so modify script to suite your needs.
#
# IgBlast sends its output to the screen, so redirect
# to a file to save.
#

# Change these:

# your input file
INPUT=/data/s166813/Downloads/Test4.txt

organism="human"
#organism="Mus_musculus"

seqType="Ig"
#seqType="TCR"

domainSystem="imgt"
#domainSystem="kabat"

# output format
# use this for rep_char
OUTFMT="7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"
# use this for changeo
#OUTFMT="7 std qseq sseq btop"

# build up the arguments, generally do not need change this
ARGS="-query $INPUT"
ARGS="$ARGS -ig_seqtype $seqType"
if [ "$seqType" == "TCR" ]; then seqType="TR"; fi  
if [ "$seqType" == "Ig" ]; then seqType="IG"; fi  
ARGS="$ARGS -organism $organism"
ARGS="$ARGS -auxiliary_data $IGDATA/optional_file/${organism}_gl.aux"
ARGS="$ARGS -germline_db_V $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${organism}_${seqType}_V.fna"
ARGS="$ARGS -germline_db_D $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${organism}_${seqType}_D.fna"
ARGS="$ARGS -germline_db_J $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${organism}_${seqType}_J.fna"
ARGS="$ARGS -domain_system $domainSystem"
ARGS="$ARGS -outfmt "

# run igblast
igblastn $ARGS "$OUTFMT"
