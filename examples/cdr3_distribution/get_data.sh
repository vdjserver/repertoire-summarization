#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: 1371444213709729305-242ac11c-0001-012
# repertoire_id: 8064755271837946346-242ac113-0001-012, 8118141715327226346-242ac113-0001-012

# If the cache is regenerated then download URLs might have changed.

curl -o study.tar https://vdjserver.tapis.io/v3/files/postits/redeem/b7e828ca-f22f-4391-ba24-c183f7432ea9-010
tar xvf study.tar 2648490830777881066-242ac113-0001-012.airr.tsv.gz
tar xvf study.tar 2669106673798681066-242ac113-0001-012.airr.tsv.gz
tar xvf study.tar repertoires.airr.json

gunzip 2648490830777881066-242ac113-0001-012.airr.tsv.gz
gunzip 2669106673798681066-242ac113-0001-012.airr.tsv.gz

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
