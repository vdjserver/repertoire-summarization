#!/bin/bash -x 
NUM_SIM=${1:-0}
echo "Using NUM_SIM ${NUM_SIM}"
#export VDJSERVER_ROOT=`dirname /home/data/vdj_server/repertoire-summarization/rep_char.py`
#get this from the existing environment . put in bashrc
#VDJ_DB_ROOT=/home/data/DATABASE/07_11_2014
VDJ_DB_ROOT=/home/data/DATABASE/12_09_2014
IGDATA=/usr/local/igblast-1.4.0
NEW_PYTHONPATH=/home/data/vdj_server/vdjml/python/
echo "Adding new PYTHON PATH $NEW_PYTHONPATH"
export PYTHONPATH=$PYTHONPATH:$NEW_PYTHONPATH
echo "PYTHONPATH IS NOW : $PYTHONPATH"
export IGDATA=$IGDATA
echo "USING IGDATA $IGDATA"
echo "USING VDJ_DB_ROOT $VDJ_DB_ROOT"
echo "START TIME"
#IGBLAST_GLOBAL_PARAMS=" -penalty -1 "
date
#http://unix.stackexchange.com/questions/56837/how-to-test-if-a-variable-is-defined-at-all-in-bash-prior-to-version-4-2-with-th
if [ -n "${ORGANISM_SET+1}" ]; then
  #echo "foobar is defined"
  echo "Using (environmental) organism set ${ORGANISM_SET}" ;
else
  #echo "foobar is not defined"
  ORGANISM_SET="human Mus_musculus" ;
fi ;
for ORGANISM in $ORGANISM_SET ;
do
	echo "Now analyzing for organism=$ORGANISM" ;
	if [ "$ORGANISM" == "Mus_musculus" ];
	then
		AUX_PATH=$IGDATA/optional_file/mouse_gl.aux ;
	else
		AUX_PATH=$IGDATA/optional_file/human_gl.aux ;
	fi ;
	BLAST_DB_ORG_ROOT=$VDJ_DB_ROOT/$ORGANISM/ReferenceDirectorySet/${ORGANISM}_
	if [ -n "${SEQ_TYPE_SET+1}" ]; then
		#it's defined
		echo "Using (environmental) SEQ_TYPE_SET '${SEQ_TYPE_SET}'"
	else
		#undefined
		SEQ_TYPE_SET="IG TR"
	fi ;
	for SEQ_TYPE in $SEQ_TYPE_SET ; 
	do
		echo -e "\n\n\n\n\n******************************************"
		QRY=$ORGANISM.$SEQ_TYPE.fna
		echo "USING query file = $QRY "
		for DCMODE in imgt kabat ; 
		do
			if [ "$SEQ_TYPE" == "TR" ]; then
				echo "SEQ_TYPE has the value 'TR'"
				IGB_SEQ_FLAG="TCR"			
			else
				echo "SEQ_TYPE has the value 'IG'"
				IGB_SEQ_FLAG="Ig"
			fi

			#echo "Want work for mode=$DCMODE" ;
			OUTPUT=$QRY.igblast.$DCMODE.out
			RC_IN=$OUTPUT
			OUTPUT_HUMAN=$OUTPUT.human_readable.txt
			BLAST_DB_ROOT=${BLAST_DB_ORG_ROOT}$SEQ_TYPE
			DB_V=${BLAST_DB_ROOT}_V.fna
			DB_D=${BLAST_DB_ROOT}_D.fna
			DB_J=${BLAST_DB_ROOT}_J.fna
			RC_VDJML=$RC_IN.rc.vdjml
			RC_COMBO=$RC_IN.combo.json
			RC_CDR3_HIST=$RC_IN.cdr3_hist.tsv
			RC_CDR3_SEG_JSON=$RC_IN.segments.json
			RC_SAMP_JSON=$RC_IN.sample.json
			RC_OUT=$RC_IN.rc_out.tsv
			RC_DB=$VDJ_DB_ROOT
			RC_LOG_OUT=${RC_OUT}.log
			RC_LOG_ERR=${RC_OUT}.err
			echo "BLAST_DB_ROOT is $BLAST_DB_ROOT"
			echo "Now analyzing seq type $SEQ_TYPE for $ORGANISM" ; 
			IGBLAST_EXEC=$IGDATA/bin/igblastn

			echo "Running VDJ sim for data..."
			SIM_DATA=$ORGANISM.$SEQ_TYPE.$DCMODE.sim.fna
			SIM_D=${VDJ_DB_ROOT}/${ORGANISM}/ReferenceDirectorySet/${ORGANISM}_${SEQ_TYPE}_D.fna
			SIM_V=${VDJ_DB_ROOT}/${ORGANISM}/ReferenceDirectorySet/${ORGANISM}_${SEQ_TYPE}_V.fna
			SIM_J=${VDJ_DB_ROOT}/${ORGANISM}/ReferenceDirectorySet/${ORGANISM}_${SEQ_TYPE}_J.fna
			../vdj_sim.py -num_seqs ${NUM_SIM} -dfasta ${SIM_D} ${SIM_V} ${SIM_J} > ${SIM_DATA}
			SIM_B_OUT=${SIM_DATA}.igblast.out
			SIM_B_OUT_HR=${SIM_B_OUT}.human_readable.txt
			SIM_JSON=$SIM_B_OUT.json
			SIM_CDR3=$SIM_B_OUT.cdr3_hist.tsv
			SIM_VDJML=$SIM_B_OUT.vdjml
			SIM_COMBO=$SIM_B_OUT.combo.json
			SIM_SAMP_JSON=$SIM_B_OUT.sample.json
			SIM_RC_OUT=$SIM_B_OUT.rc_out.tsv
			SIM_RC_OUT_LOG=${SIM_RC_OUT}.log
			SIM_RC_OUT_ERR=${SIM_RC_OUT}.err
			REP_CHAR_OPTS=" -db_species_vdjml $ORGANISM   -igblast_version igblast-1.4.0 -db_uri file://${VDJ_DB_ROOT} "
			if [ "$NUM_SIM" -eq "0" ] ;
			then 
				echo "Skipping rep_char on $SIM_B_OUT because the number of simulated sequences is $NUM_SIM" ;
			else
				echo "Running IgBLAST for sim data $SIM_DATA ..."
				echo "First in human-readable form..."
				$IGBLAST_EXEC  $IGBLAST_GLOBAL_PARAMS  -num_threads 6   -domain_system $DCMODE  -query ${SIM_DATA} -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH  -show_translation   > $SIM_B_OUT_HR &
				$IGBLAST_EXEC $IGBLAST_GLOBAL_PARAMS   -num_threads 6   -domain_system $DCMODE  -query ${SIM_DATA} -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop'  |tee $SIM_B_OUT| ../rep_char.py  ${REP_CHAR_OPTS}  -combo_out $SIM_COMBO    -sample_json_out  $SIM_SAMP_JSON     -json_out $SIM_JSON -cdr3_hist_out $SIM_CDR3 /dev/stdin $SIM_VDJML $SIM_RC_OUT $RC_DB ${SIM_DATA} $ORGANISM  1>$SIM_RC_OUT_LOG 2>$SIM_RC_OUT_ERR
				wait
			fi ;
			if [ -f "$QRY" ]  ;
			then

				echo "FOUND $QRY" ;
				

				#FORM IGBLAST HR				
				IGBLAST_HR="$IGBLAST_EXEC  $IGBLAST_GLOBAL_PARAMS  -domain_system $DCMODE  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -show_translation" ;
				#EXECUTE IT
				$IGBLAST_HR > $OUTPUT_HUMAN & 
				
				#FORM REP_CHAR for reading from stdin
				REP_CHAR_CMD="../rep_char.py   ${REP_CHAR_OPTS}  -combo_out $RC_COMBO   -sample_json_out  $RC_SAMP_JSON     -json_out $RC_CDR3_SEG_JSON -cdr3_hist_out $RC_CDR3_HIST /dev/stdin $RC_VDJML $RC_OUT $RC_DB $QRY $ORGANISM"

				#run the IgBLAST/tee/rep_char pipeline
				$IGBLAST_EXEC  $IGBLAST_GLOBAL_PARAMS   -num_threads 6   -domain_system $DCMODE  -query $QRY -germline_db_V $DB_V -germline_db_D $DB_D -germline_db_J $DB_J  -ig_seqtype $IGB_SEQ_FLAG -auxiliary_data $AUX_PATH -outfmt '7 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop' |tee $OUTPUT|$REP_CHAR_CMD   1>$RC_LOG_OUT 2>$RC_LOG_ERR ;
				echo "IgBLAST and rep_char completed at " ;
				date ;
				wait ;
			else
				echo "QUERY $QRY NOT FOUND" ;
				echo 'So no analysis for it being done !'
			fi ;

		done ;
	done ;
done ;

echo "END TIME "
date


