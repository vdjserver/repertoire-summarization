#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: PRJEB18631
# repertoire_id: 4869796857698062826-242ac113-0001-012, 4684511968548622826-242ac113-0001-012, 4746531296302862826-242ac113-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"4869796857698062826-242ac113-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"4684511968548622826-242ac113-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"4746531296302862826-242ac113-0001-012"}'

curl -o 4869796857698062826-242ac113-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/6b4ffbc9-7c1b-4f17-b5be-e44818d1ed47-010
gunzip 4869796857698062826-242ac113-0001-012.airr.tsv.gz
curl -o 4684511968548622826-242ac113-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/ef1263a3-09ed-4078-9662-1ce34161a0ac-010
gunzip 4684511968548622826-242ac113-0001-012.airr.tsv.gz
curl -o 4746531296302862826-242ac113-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/f3dbb507-f8e9-4bd9-979b-fc8977044ed0-010
gunzip 4746531296302862826-242ac113-0001-012.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# this actually pulls all metadata for the study but technically only the one repertoire
# is needed. This allows us to test that repcalc only processes the one repertoire.
curl -k -H 'content-type: application/json' --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"PRJEB18631"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_mouse_germline.airr.json
