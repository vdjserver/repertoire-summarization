#
# Simple script to run RepCalc within the docker image
#
# It has to create a JSON configuration with numerous
# setting so modify script to suite your needs.
#

# Change these:

StudyMetadata=study_metadata.json

organism="human"
#organism="Mus_musculus"

seqtype="Ig"
#seqtype="TCR"

# input
VDJMLFilesMetadata="6624457842181926425-242ac11c-0001-012 6544013104727846425-242ac11c-0001-012"
SummaryFilesMetadata="6659676574009126425-242ac11c-0001-012 6564027652327206425-242ac11c-0001-012"
ChangeOFilesMetadata="6600964371072806425-242ac11c-0001-012 6515580421228326425-242ac11c-0001-012"
JobSelected="1810667911394946585-242ac11b-0001-007"

# operations
GeneSegmentFlag=1
GeneSegmentOperations="absolute relative combo"
GeneSegmentLevels="vj vdj"
GeneSegmentFilters="productive"

CDR3Flag=1
CDR3Operations="absolute length shared"
CDR3Levels="aa v,aa vj,aa"
CDR3Filters="productive"

DiversityFlag=1
DiversityOperations="shannon"
DiversityLevels="gene aa"
DiversityFilters="productive"

ClonalFlag=0
ClonalOperations="abundance"
ClonalFilters="productive"

MutationalFlag=0
MutationalFilters="productive"

# Now run commands:

# Create initial JSON
repcalc_create_config --init $StudyMetadata $organism $seqtype repcalc_config.json

# add input files to config
vdjmlMeta=($VDJMLFilesMetadata)
summaryMeta=($SummaryFilesMetadata)
changeoMeta=($ChangeOFilesMetadata)
count=0
while [ "x${vdjmlMeta[count]}" != "x" ]
do
	filekey="file${count}"
	repcalc_create_config repcalc_config.json --file $filekey ${vdjmlMeta[count]} ${summaryMeta[count]} ${changeoMeta[count]}
	count=$(( $count + 1 ))
done

# add samples to input config
if [ -n "$JobSelected" ]; then
	repcalc_create_config repcalc_config.json --samples ${JobSelected}
fi

# Gene segment usage
if [[ $GeneSegmentFlag -eq 1 ]]; then
	ARGS=""
	if [ -n "$GeneSegmentLevels" ]; then
	    ARGS="${ARGS} --geneSegmentLevels $GeneSegmentLevels"
	fi
	if [ -n "$GeneSegmentOperations" ]; then
	    ARGS="${ARGS} --geneSegmentOperations $GeneSegmentOperations"
	fi
	if [ -n "$GeneSegmentFilters" ]; then
	    ARGS="${ARGS} --geneSegmentFilters $GeneSegmentFilters"
	fi
	cp repcalc_config.json repcalc_segment_config.json
	repcalc_create_config repcalc_segment_config.json $ARGS

    # Run RepCalc
	echo "Calculate gene segment usage"
	repcalc repcalc_segment_config.json --output output_repcalc_segment_config.json
fi

# CDR3
if [[ $CDR3Flag -eq 1 ]]; then
	rm -f joblist
	touch joblist

	ARGS=""
	if [ -n "$CDR3Levels" ]; then
	    ARGS="${ARGS} --cdr3Levels $CDR3Levels"
	fi
	if [ -n "$CDR3Operations" ]; then
	    ARGS="${ARGS} --cdr3Operations $CDR3Operations"
	fi
	if [ -n "$CDR3Filters" ]; then
	    ARGS="${ARGS} --cdr3Filters $CDR3Filters"
	fi

	# generate single group configs
	scriptList=$(repcalc_create_config repcalc_config.json --cdr3Single)
	#echo $scriptList

	for script in ${scriptList[@]}; do
	    repcalc_create_config $script $ARGS
	done
fi
