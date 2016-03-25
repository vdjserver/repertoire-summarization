#
# Trp/Phe check with a variety of sequence inputs
#
# This test is designed to be run inside igblast-app docker image.
# bash test_Trp_Phe_check.sh
#

#
# Run igblast
#

# input file
INPUT=TRP.fasta
organism="human"
seqType="Ig"
domainSystem="imgt"

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
OUTFMT="7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop"

# run igblast
igblastn $ARGS "$OUTFMT" >test_TRP.igblast.out
if [ $? -ne 0 ]; then
    echo "Error running igblastn"
    echo "*** TEST FAILED ***"
    exit
fi

#
# Run rep_char
#

# your input file and IgBlast results file
INPUT=TRP.fasta
IGBLAST_FILE=test_TRP.igblast.out

organism="human"

# build up the arguments, generally do not need to change this
RCARGS=""
RCARGS="$RCARGS -combo_out test_TRP.combo.json" 
RCARGS="$RCARGS -json_out test_TRP.segment_counts.json"
RCARGS="$RCARGS -cdr3_hist_out test_TRP.cdr3_hist.tsv" 
RCARGS="$RCARGS -sample_json_out test_TRP.sample.json" 
RCARGS="$RCARGS $IGBLAST_FILE" 
RCARGS="$RCARGS test_TRP.out.vdjml"
RCARGS="$RCARGS test_TRP.out.rc_out.tsv"
RCARGS="$RCARGS $VDJ_DB_ROOT" 
RCARGS="$RCARGS $INPUT" 
RCARGS="$RCARGS $organism" 

# run rep_char
/repsum-root/rep_char.py $RCARGS

diff test_TRP.out.rc_out.tsv gold.test_TRP.rc_out.tsv >/dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "test_TRP.out.rc_out.tsv output does not match gold.test_TRP.rc_out.tsv"
    echo "*** TEST FAILED ***"
else
    echo "*** TEST PASSED ***"
fi
