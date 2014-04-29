#!/bin/bash
LOOP_NUM=1
LINES_ACCUM=0
while [ -n "$1" ]
	do
	#echo "Current Parameter: $1 , Remaining $#"
	#echo "LOOP_NUM $LOOP_NUM"
	INFILE=$1
	if [ $LOOP_NUM == "1" ] ; then
		#just cat the file!
		cat $INFILE
		LINES_ACCUM=`wc -l $INFILE`
		LINES_ACCUM=`expr $LINES_ACCUM - 1`
	else
		#cat the file BUT use awk to ignore the first line (NR==1)
		cat $INFILE |awk '{if(NR>1) {print $0}}'
	fi ;
	#Pass $1 to some bash function or do whatever
	shift
	LOOP_NUM=`expr $LOOP_NUM + 1`
done



