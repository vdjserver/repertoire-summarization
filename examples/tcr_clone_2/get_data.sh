#
# This pulls data from VDJServer ADC and the download cache
#
# study_id: PRJNA300878
# repertoire_id: 

# If the cache is regenerated then download URLs might have changed. Need to
# query the metadata entries from the cache to update

# download url for rearrangement data
# metadata-list -v -Q '{"name":"adc_cache_study","value.study_id":"PRJNA300878"}'

# get full study data
curl -o study.tar https://vdj-agave-api.tacc.utexas.edu/postits/v2/51bfa1ed-2748-4703-8666-fac3f7735189-010
tar xf study.tar
gunzip *.airr.tsv.gz

