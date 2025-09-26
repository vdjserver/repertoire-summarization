#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: PRJEB18631
# repertoire_id: 4869796857698062826-242ac113-0001-012, 4684511968548622826-242ac113-0001-012, 4746531296302862826-242ac113-0001-012

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

curl -o study.tar https://vdjserver.tapis.io/v3/files/postits/redeem/7fde810f-83d2-4751-a5ab-25e7ec6f734d-010
tar xvf study.tar 4869796857698062826-242ac113-0001-012.airr.tsv.gz
tar xvf study.tar 4684511968548622826-242ac113-0001-012.airr.tsv.gz
tar xvf study.tar 4746531296302862826-242ac113-0001-012.airr.tsv.gz
tar xvf study.tar repertoires.airr.json

gunzip 4869796857698062826-242ac113-0001-012.airr.tsv.gz
gunzip 4684511968548622826-242ac113-0001-012.airr.tsv.gz
gunzip 4746531296302862826-242ac113-0001-012.airr.tsv.gz

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_mouse_germline.airr.json
