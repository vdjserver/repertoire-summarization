#
# The ADC doesn't have mutation data currently, so we will get data
# from our Regeneron COVVAX project. This is currently private data.

# VDJServer project: c0aef9ef-362d-46fd-abba-39ac0945bd66
# RepCalc job id: d9c76906-f138-49d7-814f-490f3a7f7124-007

# You will need to manually run these commands

# vdjserver-tools files download /projects/c0aef9ef-362d-46fd-abba-39ac0945bd66/analyses/d9c76906-f138-49d7-814f-490f3a7f7124-007/d9c76906-f138-49d7-814f-490f3a7f7124-007.zip
# vdjserver-tools files download /projects/c0aef9ef-362d-46fd-abba-39ac0945bd66/files/repertoires.airr.json

unzip d9c76906-f138-49d7-814f-490f3a7f7124-007.zip

cp d9c76906-f138-49d7-814f-490f3a7f7124-007/ad1fbdf8-35e1-4848-baea-fcabff0e9b5a.igblast.makedb.gene.clone.mutations.airr.tsv.gz .
cp d9c76906-f138-49d7-814f-490f3a7f7124-007/d3d1d40d-1192-49f2-98b2-cfc7c3064160.igblast.makedb.gene.clone.mutations.airr.tsv.gz .
cp d9c76906-f138-49d7-814f-490f3a7f7124-007/1216ddee-ac62-4069-b2f4-ee2b144ea4f2.igblast.makedb.gene.clone.mutations.airr.tsv.gz .
gunzip ad1fbdf8-35e1-4848-baea-fcabff0e9b5a.igblast.makedb.gene.clone.mutations.airr.tsv.gz
gunzip d3d1d40d-1192-49f2-98b2-cfc7c3064160.igblast.makedb.gene.clone.mutations.airr.tsv.gz
gunzip 1216ddee-ac62-4069-b2f4-ee2b144ea4f2.igblast.makedb.gene.clone.mutations.airr.tsv.gz

# TODO: we should be able to download this
# get germline database
# cp /from/some/place vdjserver_human_germline.airr.json
