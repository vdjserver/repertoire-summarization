#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: PRJNA362309
# repertoire_id: 8418294783347593706-242ac116-0001-012, 8434100262996873706-242ac116-0001-012, 8617753064573833706-242ac116-0001-012

# If the cache is regenerated then download URLs might have changed.

curl -o study.tar https://vdjserver.tapis.io/v3/files/postits/redeem/78bf628a-d16e-4f39-9e5c-1d92a46a6495-010
tar xvf study.tar
gunzip 8418294783347593706-242ac116-0001-012.airr.tsv.gz
gunzip 8434100262996873706-242ac116-0001-012.airr.tsv.gz
gunzip 8617753064573833706-242ac116-0001-012.airr.tsv.gz

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_mouse_germline.airr.json
