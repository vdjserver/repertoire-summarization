#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: 1371444213709729305-242ac11c-0001-012
# repertoire_id: 8064755271837946346-242ac113-0001-012, 8118141715327226346-242ac113-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"8064755271837946346-242ac113-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"8118141715327226346-242ac113-0001-012"}'

curl -o 8064755271837946346-242ac113-0001-012.airr.tsv.gz https://vdjserver.org/postits/v2/732e64e5-1328-46c5-9fb6-1747a15b9f46-010
gunzip 8064755271837946346-242ac113-0001-012.airr.tsv.gz
curl -o 8118141715327226346-242ac113-0001-012.airr.tsv.gz https://vdjserver.org/postits/v2/51c311e1-aa0b-4be0-9a9a-aa7e9d6e5da0-010
gunzip 8118141715327226346-242ac113-0001-012.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# this actually pulls all metadata for the study but technically only the one repertoire
# is needed. This allows us to test that repcalc only processes the one repertoire.
curl -k -H 'content-type:application/json' --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"1371444213709729305-242ac11c-0001-012"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
