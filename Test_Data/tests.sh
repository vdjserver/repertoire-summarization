#!/bin/bash -x 
IGBLAST_EXEC=/usr/local/bin/igblastn
VDJ_DB_ROOT=/home/data/DATABASE/01_22_2014
BLAST_DB_ROOT=$VDJ_DB_ROOT/human/ReferenceDirectorySet/human_IG
DB_V=${BLAST_DB_ROOT}_V.fna
DB_D=${BLAST_DB_ROOT}_D.fna
DB_J=${BLAST_DB_ROOT}_J.fna
QRY=AY671579.1.fna
IGDATA=/usr/local/igblast_from_lonestar/
AUX_PATH=$IGDATA/optional_file/human_gl.aux
OUTPUT=igblast.out
OUTPUT_HUMAN=$OUTPUT.human_readable.txt
echo "Running IgBLASTn ..."
date
$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'    > $OUTPUT
echo "Running IgBLASTn (human readble)..."
date
$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J -auxiliary_data $AUX_PATH -show_translation >$OUTPUT_HUMAN
echo "IgBLASTn finished"
date
echo "Running repertoire characterization..."
RC_IN=$OUTPUT
RC_VDJML=$RC_IN.rc.vdjml
RC_CDR3_HIST=$RC_IN.cdr3_hist.tsv
RC_CDR3_SEG_JSON=$RC_IN.segments.json
RC_OUT=$RC_IN.rc_out.tsv
RC_DB=$VDJ_DB_ROOT
RC_LOG_OUT=${RC_OUT}.out
RC_LOG_ERR=${RC_OUT}.err
../rep_char.py -json_out $RC_CDR3_SEG_JSON -cdr3_hist_out $RC_CDR3_HIST $RC_IN $RC_VDJML $RC_OUT $RC_DB $QRY 1>$RC_LOG_OUT 2>$RC_LOG_ERR
echo "Repertoire characterization complete."
date


