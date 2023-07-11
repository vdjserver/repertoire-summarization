#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: PRJNA362309
# repertoire_id: 8418294783347593706-242ac116-0001-012, 8434100262996873706-242ac116-0001-012, 8617753064573833706-242ac116-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"8418294783347593706-242ac116-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"8434100262996873706-242ac116-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"8617753064573833706-242ac116-0001-012"}'

curl -o 8418294783347593706-242ac116-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/bc5ab518-8100-4320-87a1-7679d21a3015-010
gunzip 8418294783347593706-242ac116-0001-012.airr.tsv.gz
curl -o 8434100262996873706-242ac116-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/0a36a167-bb53-4cda-8bf4-5effa656a83f-010
gunzip 8434100262996873706-242ac116-0001-012.airr.tsv.gz
curl -o 8617753064573833706-242ac116-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/f7c96efe-f356-4f79-a4ef-de16397cdb42-010
gunzip 8617753064573833706-242ac116-0001-012.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# this actually pulls all metadata for the study but technically only the one repertoire
# is needed. This allows us to test that repcalc only processes the one repertoire.
curl -k -H 'content-type: application/json' --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"PRJNA362309"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_mouse_germline.airr.json
