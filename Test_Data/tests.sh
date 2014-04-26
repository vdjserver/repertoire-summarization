#!/bin/bash -x 
VDJ_DB_ROOT=/home/data/DATABASE/01_22_2014
#VDJ_DB_ROOT=/home/data/DATABASE/04_22_2014
#IGDATA=/usr/local/igblast-1.2.0
IGDATA=/usr/local/igblast-1.3.0
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
		echo "USING query file = $QRY "
		for DCMODE in kabat imgt ; 
		do
			if [ "$SEQ_TYPE" == "TR" ]; then
				echo "SEQ_TYPE has the value 'TR'"
				IGB_SEQ_FLAG="TCR"			
			else
				echo "SEQ_TYPE has the value 'IG'"
				IGB_SEQ_FLAG="Ig"
			fi
			if [ -f "$QRY" ]  ;
			then
				#echo "Want work for mode=$DCMODE" ;
				OUTPUT=$QRY.igblast.$DCMODE.out
				RC_IN=$OUTPUT
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
				RC_LOG_OUT=${RC_OUT}.log
				RC_LOG_ERR=${RC_OUT}.err
				echo "BLAST_DB_ROOT is $BLAST_DB_ROOT"
				echo "Now analyzing seq type $SEQ_TYPE for $ORGANISM" ; 
				IGBLAST_EXEC=$IGDATA/bin/igblastn
				echo "FOUND $QRY" ;
				echo "Starting IgBLAST with outfmt 7 at "
				date
				time $IGBLAST_EXEC  -num_threads 6   -domain_system $DCMODE  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'  >$OUTPUT
				echo "IgBLAST with outfmt 7 finished at "
				date
				time echo "Now running IgBLAST for human-readable output ..."
				$IGBLAST_EXEC -domain_system $DCMODE  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -show_translation >$OUTPUT_HUMAN 
				echo "IgBLAST for human-readable output finished."
				date
				echo "Now running repertoire characterization ..."
				time ../rep_char.py -json_out $RC_CDR3_SEG_JSON -cdr3_hist_out $RC_CDR3_HIST $RC_IN $RC_VDJML $RC_OUT $RC_DB $QRY $ORGANISM  1>$RC_LOG_OUT 2>$RC_LOG_ERR
			else
				echo "QUERY $QRY NOT FOUND" ;
				echo 'So no analysis for it being done !'
			fi ;

		done ;
	done ;
done ;




