#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: 3276777473314001386-242ac116-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data and repertoire metadata
# metadata-list -v -Q '{"name":"adc_cache_study","value.study_id":"3276777473314001386-242ac116-0001-012"}'

curl -o study.tar https://vdj-agave-api.tacc.utexas.edu/postits/v2/937871ca-0d6e-4903-9d2f-c63d4c4638c4-010
tar xf study.tar
gunzip *.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# not needed as it is in the study archive
# curl -k --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"3276777473314001386-242ac116-0001-012"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
