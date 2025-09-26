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

# BCR (High-Throughput Mapping of B Cell Receptor Sequences to Antigen Specificity)
curl -o study.tar https://vdjserver.tapis.io/v3/files/postits/redeem/8b3a855e-7433-47f0-a1af-ce71c2dd07bf-010
tar xvf study.tar
gunzip 1159043099869244949-242ac114-0001-012.airr.tsv.gz
gunzip 1159043104164212245-242ac114-0001-012.airr.tsv.gz
gunzip 1159043095574277653-242ac114-0001-012.airr.tsv.gz

# query VDJServer ADC to get repertoire metadata
# this actually pulls all metadata for the study but technically only the one repertoire
# is needed. This allows us to test that repcalc only processes the one repertoire.
curl -k -H 'content-type:application/json' --data '{"filters":{"op":"=","content":{"field":"study.study_id","value":"PRJNA578389"}}}' https://vdjserver.org/airr/v1/repertoire | jq '.' > repertoire.airr.json

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
