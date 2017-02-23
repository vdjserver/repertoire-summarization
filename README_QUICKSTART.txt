#######################################################################
# README QUICK_START FOR repertoire-summarization SUITE of VDJ SERVER #
#######################################################################

Verify/perform the following for quick start and testing!

1) make sure that pyVDJML is installed.
	Download and install VDJML 1.1.0 python bindings

2) make sure that IgBLAST v1.3.0 is installed.  Make sure that the offical
	NCBI V,D,J databases for mouse and human are both ready for use.  Make sure
	that its package contains not only igblastn, but also makeblastd and blastx.
	The coupled makeblastd and blastx should be used with the vdj_db_dl database download/update script.

3) make sure that the following python modules are installed
	Bioptyon v1.63
	BeautifulSoup4 v4.3.2

4) edit the dl_map file (if necessary) making sure that the IgBLAST variables are
correct (names/keys on the first/left columns ; values on the right/second columns).

5) run the vdj_db_dl.py script to download IMGT data and prepare a VDJ_server database 
(including region annotation information and GeneTable information).

6) re-run the vdj_db_dl.py script but this time using the "-analyze_only" flag to analyze
the gene tables.  Apply patches as necessary and re-run the "analyze_only" (as often as needed) to verify
the patches are used as desired.  A handful (e.g. 4 or 5) patches may be necessary.

6) enter into the Test_Data directory and read the tests_README about editing
the tests.sh script environment variables (as needed ; IGDATA and VDJ_DB_ROOT are
the variables most likely to need editing).

7) run in that directory the tests.sh script to exercise the code.  Note any
error output as well as the contents of any *.err files or *.out files.
If no errors are seen in either the output of the .err or the .out files (from igblast) then 
execution has succeded without errors.  See the *.tsv file for read-characterization
data, CDR3 length/histogram data, the *.json files for segment/hit count data, and
the *.vdjml files for VDJML output.

