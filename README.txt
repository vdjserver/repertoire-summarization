#########################################################
#README FOR repertoire-summarization SUITE of VDJ SERVER#
#########################################################

The purpose of this suite of software tools is to generate 
summary statistics of REP-SEQ data.

IGBLAST output files (using a certain output format) are 
taken as input and 4 distinct files are generated as outputs:
1) CDR3 length histogram (in kabat and imgt modes),
2) JSON-formatted hierachies of the IMGT alleles and 
corresponding counts of segments in top combinations
with agregrated counts, 3) read-level summary statistics
, and 4 ) vdjml files of the parsed data.

The rep_char.py script is designed to support the above 4 functions.
If run over several igblast outputs that should subsequently be 
merged the following scripts can be used to logicaly merge outputs:
1) cdr3_hist.py can logically merge CDR3 length histogram data,
2) json_seg_merge.py can be used to logically merge JSON count data,
3) vdjml_merge.py can be used to merge VDJML files, and 4) merge_tables.py
can be used to merge TSV tables.

To get help on any of the *.pyscripts run it as seen here, and 
use the "-h" flag for "help"!  For the .sh scripts please open them and read the comments.

NOTE : use of rep_char.py require that the python module "vdjml" be usable.

NOTE : rep_char.py expects to be able to find the file"vdjml_igblast_parse.py". 
Successful finding may be achieved by using the PYTHONPATH

esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ ./rep_char.py -h
Traceback (most recent call last):
  File "./rep_char.py", line 4, in <module>
    from vdjml_igblast_parse import scanOutputToVDJML,makeParserArgs,makeVDJMLDefaultMetaAndFactoryFromArgs
esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ export PYTHONPATH=/home/data/vdj_server/vdjml/python/
esalina2@eddiecomp:/home/data/vdj_server/repertoire-summarization$ ./rep_char.py -h
usage: rep_char.py [-h] [-igblast_version IGBLAST_VERSION]
                   [-igblast_params IGBLAST_PARAMS]
                   [-igblast_runid IGBLAST_RUNID] [-db_name DB_NAME]
                   [-db_ver DB_VER] [-db_species {human,Mus_musculus}]
                   [-db_uri DB_URI] [-igblast_dc {imgt,kabat}]
                   [-json_out JSON_OUT] [-cdr3_hist_out CDR3_HIST_OUT]
                   [-skip_char]
                   igblast_in vdjml_out char_out vdj_db qry_fasta
                   {human,Mus_musculus}

Parse IGBLAST output, write a VDJML output file. Generate read-level
repertoire-characterization data.

positional arguments:
  igblast_in            path to igblast analysis file of results, hits,
                        segments, etc. *** CRITICAL NOTE *** SHOULD BE RUN
                        WITH OUTFMT AS SEEN HERE -outfmt '7 qseqid qgi qacc
                        qaccver qlen sseqid sallseqid sgi sallgi sacc saccver
                        sallacc slen qstart qend sstart send qseq sseq evalue
                        bitscore score length pident nident mismatch positive
                        gapopen gaps ppos frames qframe sframe btop'
  vdjml_out             the path to the output vdjml
  char_out              the path to the output TSV file of read-level
                        repertoire characterization data
  vdj_db                path to the VDJ root REQUIRED
  qry_fasta             path to the input fasta file of query (Rep-Seq) data
                        input to IgBLAST
  {human,Mus_musculus}  db organism

optional arguments:
  -h, --help            show this help message and exit
  -igblast_version IGBLAST_VERSION
                        the version of igblast used
  -igblast_params IGBLAST_PARAMS
                        the parameters passed to igblast
  -igblast_runid IGBLAST_RUNID
                        the ID assigned to the run
  -db_name DB_NAME      the name of the IgBLAST database
  -db_ver DB_VER        a version string associated with the database
  -db_species {human,Mus_musculus}
                        db organism for VDJML output file(s)
  -db_uri DB_URI        a URI associated with the database
  -igblast_dc {imgt,kabat}
                        the domain classification system used
  -json_out JSON_OUT    output file for the JSON segment count IMGT hierarchy
  -cdr3_hist_out CDR3_HIST_OUT
                        output file for the CDR3 histogram of lengths (both
                        kabat and imgt systems)
  -skip_char            If this is set, then characterization besides
                        asignment and CDR3 length is skipped
 

-Eddie Salinas
VDJ Server Team Member
