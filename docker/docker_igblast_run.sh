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
INPUT=/data/s166813/Downloads/18.fasta
IGBLAST_OUT=igblast.out
VDJML_OUT=out.vdjml

organism="human"
#organism="mouse"

seqType="Ig"
#seqType="TCR"

domainSystem="imgt"
#domainSystem="kabat"

# output format
# use this for rep_char
OUTFMT="7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"
# use this for changeo
#OUTFMT="7 std qseq sseq btop"

# Make sure to update these when a new IgBlast or DB is used
IGBLAST_VERSION=1.4.0
IGBLAST_URI=http://www.ncbi.nlm.nih.gov/projects/igblast/
VDJ_DB_VERSION=07_11_2014
VDJ_DB_URI=http://wiki.vdjserver.org/vdjserver/index.php/VDJServer_IgBlast_Database

######

# build up the arguments, generally do not need change anything below this line
ARGS="-query $INPUT"
ARGS="$ARGS -ig_seqtype $seqType"
if [ "$seqType" == "TCR" ]; then seqType="TR"; fi  
if [ "$seqType" == "Ig" ]; then seqType="IG"; fi  
ARGS="$ARGS -organism $organism"
ARGS="$ARGS -auxiliary_data $IGDATA/optional_file/${organism}_gl.aux"
if [ "$organism" == "mouse" ]; then organism="Mus_musculus"; fi  
ARGS="$ARGS -germline_db_V $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${organism}_${seqType}_V.fna"
ARGS="$ARGS -germline_db_D $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${organism}_${seqType}_D.fna"
ARGS="$ARGS -germline_db_J $VDJ_DB_ROOT/$organism/ReferenceDirectorySet/${organism}_${seqType}_J.fna"
ARGS="$ARGS -domain_system $domainSystem"
ARGS="$ARGS -outfmt "

# run igblast
igblastn $ARGS "$OUTFMT" > $IGBLAST_OUT

# generate vdjml
RCARGS="$IGBLAST_OUT $VDJML_OUT"
RCARGS="$RCARGS -db_name_igblast ${organism}_${seqType}"
RCARGS="$RCARGS -db_ver $VDJ_DB_VERSION"
RCARGS="$RCARGS -db_species_vdjml $organism"
RCARGS="$RCARGS -db_uri $VDJ_DB_URI"
RCARGS="$RCARGS -igblast_version $IGBLAST_VERSION"
RCARGS="$RCARGS -igblast_uri $IGBLAST_URI"
IGBLAST_PARAMS="$ARGS"

igblast_parse.py $RCARGS -igblast_params "$IGBLAST_PARAMS"
