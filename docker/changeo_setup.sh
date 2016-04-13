#
# Some additional steps to setup docker image for Changeo
#

# fasta all in one file for MakeDb
cd $VDJ_DB_ROOT/human/ReferenceDirectorySet
rm -f IG_VDJ.fna TR_VDJ.fna
cat IG*.fna >IG_VDJ.fna
cat TR*.fna >TR_VDJ.fna

# fasta all in one file for MakeDb
cd $VDJ_DB_ROOT/Mus_musculus/ReferenceDirectorySet
rm -f IG_VDJ.fna TR_VDJ.fna
cat IG*.fna >IG_VDJ.fna
cat TR*.fna >TR_VDJ.fna

echo 'install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "alakazam", "shazam"), repos=c("https://cran.revolutionanalytics.com/"))' | R --no-save
