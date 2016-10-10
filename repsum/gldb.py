"""
Germline Database
"""

import glob
import pickle
import os
import sys

# repsum modules
import imgt_utils as imgt
import utils

GERMLINE_DB_ROOT = None
IMGT_DB = None

def init_germline_db_root(db_path=None):
    """Determine the path to germline database files"""
    global GERMLINE_DB_ROOT
    global IMGT_DB

    if db_path is None:
        if "VDJ_DB_ROOT" in os.environ: GERMLINE_DB_ROOT = os.environ['VDJ_DB_ROOT']
    else: GERMLINE_DB_ROOT = db_path
    if GERMLINE_DB_ROOT is None:
        print("ERROR: Germline database path is not defined.")
        sys.exit()
    else:
        if not os.path.exists(GERMLINE_DB_ROOT):
            print("ERROR: Germline database path is not a valid path: " + GERMLINE_DB_ROOT)
            sys.exit()
    # TODO: additional checks

    # create IMGT database object
    IMGT_DB = imgt.imgt_db(GERMLINE_DB_ROOT)


def germline_db_root():
    """Return germline database root path"""
    if GERMLINE_DB_ROOT is None:
        print("ERROR: Germline database path is not defined.")
        sys.exit()
    else: return GERMLINE_DB_ROOT

#rooted at an organism get the hierarchy from the gene tables
def getHierarchyBy(org_name):
    geneTablesDirectoryOfHTMLFiles = germline_db_root()+"/"+org_name+"/GeneTables/"
    fullPklPath = IMGT_DB.getPickleFullPath()
    print "Trying to use " ,fullPklPath
    if(not(fullPklPath==None)):
        #print "not none"
        if(os.path.exists(fullPklPath)):
            #print "it exists"
            unpickled_data=utils.pickleRead(fullPklPath)
            #print unpickled_data
            #print "len of data is ",len(unpickled_data)
            subset_for_extraction=unpickled_data[0]
            for item in subset_for_extraction:
                if(item==org_name):
                    return subset_for_extraction[item]
                #print item
            #return unpickled_data[org_name]
    loci=imgt.get_loci_list()
    hierarchy=utils.tree()
    for locus in loci:
        htmlGeneTablesGlob=geneTablesDirectoryOfHTMLFiles+"/*"+locus+"*.html"
        geneTableHTMLFiles=glob.glob(htmlGeneTablesGlob)
        fastaAlleleList=None
        if(filterbyFastaAlleles):
            fastaFilesGlob=geneTablesDirectoryOfHTMLFiles+"/../ReferenceDirectorySet/*"+locus+"*.fna"
            print "the fasta glob is ",fastaFilesGlob
            fastaFiles=glob.glob(fastaFilesGlob)
            fastaMainMap=dict()
            for fastaFile in fastaFiles:
                fastaString=readFileIntoString(fastaFile)
                fastaList=read_fasta_string(fastaString)
                fastaMap=read_fasta_into_map(fastaList)	
                for fastaKey in fastaMap:
                    fastaMainMap[fastaKey]=fastaMap[fastaKey]
                print "done loading into main map..."
                fastaListOfNames=getIMGTNameListFromFastaMap(fastaMainMap)
                print "FASTA LIST OF NAMES:"
                fastaListOfNames.sort()
                printList(fastaListOfNames)
                print "got name list.... its length is",len(fastaListOfNames)
                fastaAlleleList=allelifyList(fastaListOfNames)
                print "The Fasta Allele list is ",fastaAlleleList
            html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[0],locus,fastaAlleleList)
            locusHierarchyData=html_data[0]
            html_data=hierarchyTreeFromGenetableURL("file://"+geneTableHTMLFiles[1],locus,fastaAlleleList,locusHierarchyData)
            locusHierarchyData=html_data[0]
            hierarchy[locus]=locusHierarchyData
    return hierarchy

