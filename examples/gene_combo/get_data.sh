#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: 3276777473314001386-242ac116-0001-012
# repertoire_id: 2648490830777881066-242ac113-0001-012, 2669106673798681066-242ac113-0001-012

# If the cache is regenerated then download URLs might have changed.

# TCR
curl -o study.tar https://vdjserver.tapis.io/v3/files/postits/redeem/b7e828ca-f22f-4391-ba24-c183f7432ea9-010
tar xvf study.tar
gunzip 2648490830777881066-242ac113-0001-012.airr.tsv.gz
gunzip 2669106673798681066-242ac113-0001-012.airr.tsv.gz

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json