#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: PRJNA549712
# repertoire_id: 6096696254733292012-242ac114-0001-012, 6079859978638004716-242ac114-0001-012, 6085958832198324716-242ac114-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"6096696254733292012-242ac114-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"6079859978638004716-242ac114-0001-012"}'
# metadata-list -v -Q '{"name":"adc_cache_repertoire","value.repertoire_id":"6085958832198324716-242ac114-0001-012"}'

curl -o 6096696254733292012-242ac114-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/8169e37d-d967-4883-a3c4-c83e4bd5eff3-010
gunzip 6096696254733292012-242ac114-0001-012.airr.tsv.gz
curl -o 6079859978638004716-242ac114-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/fed9847b-1ea2-47b0-9840-d7afc6e40216-010
gunzip 6079859978638004716-242ac114-0001-012.airr.tsv.gz
curl -o 6085958832198324716-242ac114-0001-012.airr.tsv.gz https://vdj-agave-api.tacc.utexas.edu/postits/v2/930adb3b-0106-49f9-b23a-4efed22abe47-010
gunzip 6085958832198324716-242ac114-0001-012.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# this actually pulls all metadata for the study but technically only the one repertoire
# is needed. This allows us to test that repcalc only processes the one repertoire.
curl -k -H 'content-type: application/json' --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"PRJNA549712"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
