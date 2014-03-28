#!/bin/bash -x 
VDJ_DB_ROOT=/home/data/DATABASE/01_22_2014
IGDATA=/usr/local/igblast-1.2.0
for ORGANISM in Mus_musculus human ;
do
	echo "Now analyzing for organism=$ORGANISM" ;
	if [ "$ORGANISM" == "Mus_musculus" ];
	then
		AUX_PATH=$IGDATA/optional_file/mouse_gl.aux ;
	else
		AUX_PATH=$IGDATA/optional_file/human_gl.aux ;
	fi ;
	BLAST_DB_ORG_ROOT=$VDJ_DB_ROOT/$ORGANISM/ReferenceDirectorySet/${ORGANISM}_
	for SEQ_TYPE in IG TR ; 
	do
		echo -e "\n\n\n\n\n******************************************"
		QRY=$ORGANISM.$SEQ_TYPE.fna
		RC_IN=$OUTPUT
		OUTPUT=$QRY.igblast.out
		OUTPUT_HUMAN=$OUTPUT.human_readable.txt
		BLAST_DB_ROOT=${BLAST_DB_ORG_ROOT}$SEQ_TYPE
		DB_V=${BLAST_DB_ROOT}_V.fna
		DB_D=${BLAST_DB_ROOT}_D.fna
		DB_J=${BLAST_DB_ROOT}_J.fna
		RC_VDJML=$RC_IN.rc.vdjml
		RC_CDR3_HIST=$RC_IN.cdr3_hist.tsv
		RC_CDR3_SEG_JSON=$RC_IN.segments.json
		RC_OUT=$RC_IN.rc_out.tsv
		RC_DB=$VDJ_DB_ROOT
		RC_LOG_OUT=${RC_OUT}.out
		RC_LOG_ERR=${RC_OUT}.err
		echo "BLAST_DB_ROOT is $BLAST_DB_ROOT"
		echo "Now analyzing seq type $SEQ_TYPE for $ORGANISM" ; 
		IGBLAST_EXEC=$IGDATA/bin/igblastn
		if [ "$SEQ_TYPE" == "TR" ]; then
			echo "SEQ_TYPE has the value 'TR'"
			IGB_SEQ_FLAG="TCR"			
		else
			echo "SEQ_TYPE has the value 'IG'"
			IGB_SEQ_FLAG="Ig"
		fi
		if [ -f "$QRY" ]  ;
		then
			echo "FOUND $QRY" ;
			$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'  >$OUTPUT
			$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -show_translation >$OUTPUT_HUMAN 
		else
			echo "QUERY $QRY NOT FOUND" ;
		fi ;
	done ;
done ;



#IGBLAST_EXEC=/usr/local/bin/igblastn
#VDJ_DB_ROOT=/home/data/DATABASE/01_22_2014
#BLAST_DB_ROOT=$VDJ_DB_ROOT/human/ReferenceDirectorySet/human_IG
#DB_V=${BLAST_DB_ROOT}_V.fna
#DB_D=${BLAST_DB_ROOT}_D.fna
#DB_J=${BLAST_DB_ROOT}_J.fna
#QRY=AY671579.1.fna
#IGDATA=/usr/local/igblast_from_lonestar/
#AUX_PATH=$IGDATA/optional_file/human_gl.aux
#OUTPUT=igblast.out
#OUTPUT_HUMAN=$OUTPUT.human_readable.txt
#echo "Running IgBLASTn for $QRY ..."
#date
#$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'    > $OUTPUT
#echo "Running IgBLASTn for $QRY (human readble)..."
#date
#$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J -auxiliary_data $AUX_PATH -show_translation >$OUTPUT_HUMAN
#echo "IgBLASTn finished for $QRY"
#date
#echo "Running repertoire characterization..."
#RC_IN=$OUTPUT
#RC_VDJML=$RC_IN.rc.vdjml
#RC_CDR3_HIST=$RC_IN.cdr3_hist.tsv
#RC_CDR3_SEG_JSON=$RC_IN.segments.json
#RC_OUT=$RC_IN.rc_out.tsv
#RC_DB=$VDJ_DB_ROOT
#RC_LOG_OUT=${RC_OUT}.out
#RC_LOG_ERR=${RC_OUT}.err
#../rep_char.py -json_out $RC_CDR3_SEG_JSON -cdr3_hist_out $RC_CDR3_HIST $RC_IN $RC_VDJML $RC_OUT $RC_DB $QRY 1>$RC_LOG_OUT 2>$RC_LOG_ERR
#echo "Repertoire characterization complete for $QRY"
#QRY=SRR443375.12.GQMC0HM04JC84C.fna
#IGDATA=/usr/local/igblast-1.2.0
#AUX_PATH=$IGDATA/optional_file/human_gl.aux
#OUTPUT=igblast.$QRY.out
#OUTPUT_HUMAN=$OUTPUT.human_readable.txt
#BLAST_DB_ROOT=$VDJ_DB_ROOT/human/ReferenceDirectorySet/human_TR
#DB_V=${BLAST_DB_ROOT}_V.fna
#DB_D=${BLAST_DB_ROOT}_D.fna
#DB_J=${BLAST_DB_ROOT}_J.fna
#IGBLAST_EXEC=$IGDATA/bin/igblastn
#formatted for rep_char
#$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype TCR -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'  >$OUTPUT
##human readable
#$IGBLAST_EXEC -domain_system imgt  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype TCR -auxiliary_data $AUX_PATH -show_translation >$OUTPUT_HUMAN
#RC_IN=$OUTPUT
#RC_VDJML=$RC_IN.rc.vdjml
#RC_CDR3_HIST=$RC_IN.cdr3_hist.tsv
#RC_CDR3_SEG_JSON=$RC_IN.segments.json
#RC_OUT=$RC_IN.rc_out.tsv
#RC_DB=$VDJ_DB_ROOT
#RC_LOG_OUT=${RC_OUT}.out
#RC_LOG_ERR=${RC_OUT}.err
#../rep_char.py -json_out $RC_CDR3_SEG_JSON -cdr3_hist_out $RC_CDR3_HIST $RC_IN $RC_VDJML $RC_OUT $RC_DB $QRY 1>$RC_LOG_OUT 2>$RC_LOG_ERR


date


