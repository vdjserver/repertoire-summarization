#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: 3276777473314001386-242ac116-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data and repertoire metadata
# metadata-list -v -Q '{"name":"adc_cache_study","value.study_id":"3276777473314001386-242ac116-0001-012"}'

curl -o study.tar https://vdjserver.tapis.io/v3/files/postits/redeem/b7e828ca-f22f-4391-ba24-c183f7432ea9-010
tar xf study.tar 
gunzip *.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# not needed as it is in the study archive
# curl -k --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"3276777473314001386-242ac116-0001-012"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
