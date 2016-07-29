#!/usr/bin/env python

import os
import re
import glob
from os.path import basename
from segment_utils import getFastaListOfDescs,getQueryIndexGivenSubjectIndexAndAlignment,getAdjustedCDR3StartFromRefDirSetAllele,getADJCDR3EndFromJAllele,getEmptyRegCharMap,alleleIsTR,getTheFrameForThisReferenceAtThisPosition,getTheFrameForThisJReferenceAtThisPosition
from imgt_utils import imgt_db,get_cdr3_end
from utils import *
import vdjml
from vdjml_utils import getTopVDJItems,getHitInfo
from igblast_parse import rev_comp_dna
import argparse
from alignment import alignment,CodonAnalysis,codonAnalyzer








#utility to reconstruct alignments from a BTOP and data map and add it to the 
#data map and return the data map
def addAlignmentsPreCDR3(dataMap,alleleName,imgtdb_obj,organism,query_record):
    btop=dataMap['btop']
    query_seq=str(query_record.seq)
    if(dataMap['is_inverted']):
        #print "inversion is necessary...."
        query_seq=rev_comp_dna(query_seq)
    query_seq=query_seq[int(dataMap['q. start']-1):int(dataMap['q. end'])]
    subject_seq=imgtdb_obj.getRefDirSetFNASeqGivenOrgAndAllele(alleleName,organism)
    subject_seq=subject_seq[int(dataMap['s. start']-1):int(dataMap['s. end'])]
    qsAlgn=buildAlignmentWholeSeqs(btop,query_seq,subject_seq)
    dataMap['query seq']=qsAlgn[0]
    dataMap['subject seq']=qsAlgn[2]
    return dataMap



#wrapper for getting CDR3 length using pyVDJML objects as input
def CDR3LengthAnalysisVDMLOBJ(read_result_obj,meta,organism,imgtdb_obj,query_record):
    #for V and J require:
    # 1) allele names
    # 2) alignment from BTOP reconstruction both Q and S
    # 3) q. start and q. end and s. start and s. end
    # 4) read inversion flags
    #print "got into cdr3 hist wrapper"
    segment_combinations=read_result_obj.segment_combinations()
    #if(len(segment_combinations)
    #print "the length is ",len(segment_combinations)
    VDJMap=getTopVDJItems(read_result_obj,meta)
    vAllele=VDJMap['V']
    jAllele=VDJMap['J']
    empty_map=dict()
    modes=get_domain_modes()
    for mode in modes:
        empty_map[mode]=(-1)
        empty_map[mode+'_from']=(-1)
        empty_map[mode+'_to']=(-1)
        empty_map[mode+'_to']=(-1)
        empty_map['qry_rev']=False
    if(vAllele!=None and jAllele!=None):
        vData=getHitInfo(read_result_obj,meta,vAllele)
        vData=addAlignmentsPreCDR3(vData,vAllele,imgtdb_obj,organism,query_record)
        jData=getHitInfo(read_result_obj,meta,jAllele)
        jData=addAlignmentsPreCDR3(jData,jAllele,imgtdb_obj,organism,query_record)
        #printMap(vData)
        #printMap(jData)
        cdr3_analysis_map=CDR3LengthAnalysis(vData,jData,organism,imgtdb_obj,query_record)
        #printMap(cdr3_analysis_map)
        return cdr3_analysis_map
    else:
        #print "insufficient data!"
        return empty_map
    pass


#class for CDR3 histogram
class histoMapClass:


    #counter
    count_map=dict()
    modes=None

    #constructor
    def __init__(self,modes):
        self.count_map=dict()
        self.modes=modes
        for mode in modes:
            self.count_map[mode]=dict()

    #verify integer string
    def appearsAsInt(self,s):
        ire=re.compile(r'^\-?\d+$')
        if(re.match(ire,s)):
            return True
        else:
            return False


    #read a file into the object!
    def read_from_file(self,infile):
        line_num=1
        try:
            reader=open(infile,'r')
            self.modes=list()
            for line in reader:
                line=line.strip()
                line_pieces=line.split('\t')
                #print "got line=",line," and line_pieces=",line_pieces
                if(line_num==1 and len(line_pieces)>=2):
                    for m in range(1,len(line_pieces)):
                        self.modes.append(line_pieces[m].strip())
                elif(line_num==1 and len(line_pieces)<2):
                    raise Exception('Error, invalid header in histogram file '+str(infile))
                elif(len(line_pieces)==len(self.modes)+1):
                    #got data!
                    for lp in line_pieces:
                        if(not(self.appearsAsInt(lp))):
                            raise Exception("Error, data ",lp," on line # ",line_num," appears as non-integral!  Integer values are expected!")
                    val=line_pieces[0]
                    for m in range(0,len(self.modes)):
                        #print "for val=",val," and mode=",self.modes[m]," need to inc count ",line_pieces[m+1]
                        for i in range(int(line_pieces[m+1])):
                            mode=self.modes[m]
                            #print "calling inc with mode=",mode," val=",val
                            self.inc(mode,val)
                else:
                    print "warning, bad data in file ",infile," line ",line_num
                line_num+=1             
            reader.close()
        except Exception as my_exception:
            print my_exception," line #",line_num


    #increment
    def inc(self,mode,val):
        val=int(val)
        #print "INC CALLED"
        if(not val in self.count_map[mode]):
            self.count_map[mode][val]=1
        else:
            self.count_map[mode][val]+=1

    #get min val over all vals
    def gminVal(self):
        min_val=1000
        for mode in self.modes:
            for val in self.count_map[mode]:
                if(val<min_val):
                    min_val=val
        return min_val

    #get max val over all vals
    def gmaxVal(self):
        max_val=(-1)
        for mode in self.modes:
            for val in self.count_map[mode]:
                if(val>max_val):
                    max_val=val
        return max_val
                            


    #merge from other hist
    def merge(self,other):
        #print "MERGE IS CALLED"
        other_map=other.count_map
        #print "other_map is ",other_map
        for mode in other_map:
            if(not mode in self.count_map):
                self.count_map[mode]=dict()
            for val in other_map[mode]:
                #print "for mode=",mode," at val=",val
                other_count=other_map[mode][val]
                #print "the count is ",other_count
                for i in range(other_count):
                    self.inc(mode,val)

    #print maps basic
    def printMaps(self):
        for mode in self.modes:
            #print "SHOWING MAP FOR ",mode
            for val in self.count_map[mode]:
                print val,"\t",self.count_map[mode][val]


    #write as STREAM JSON
    def JSONIFYStreamAsString(self):
        self.modes.sort()
        mode_stream_list=list()
        for mode in self.modes:
            #print "Now JSONIFying for MODE=",mode
            min_val=self.gminVal()
            max_val=self.gmaxVal()
            #print "min and max are ",min_val,max_val
            if((type(min_val)!=int) or (type(max_val)!=int)):
                min_val=(-1)
                max_val=(-1)
                print "Warning, could not determine min/max values for CDR3 lengths!  Output may not be defined!"           
            x_y_vals_array=list()
            for v in range(min(min_val,max_val),max(min_val,max_val)+1):
                x_y_str=""
                x_y_str+="\"x\": "+str(v)
                x_y_str+=","
                actual_val=0
                if(v in self.count_map[mode]):
                    actual_val=self.count_map[mode][v]
                x_y_str+="\"y\": "+str(actual_val)
                x_y_vals_array.append(x_y_str)
            obj_sep="}\n,\n{"
            x_y_vals_array_str=obj_sep.join(x_y_vals_array)
            stream_str="{\n\"key\":\""+str(mode)+"\",\n"
            stream_str+="\"values\": [ \n{"+x_y_vals_array_str+"}\n]\n}"
            mode_stream_list.append(stream_str)
        stream_sep=","
        JSON="["+stream_sep.join(mode_stream_list)+"]"
        return JSON
#The sample data is 
#[
# {
#  "key": "KABAT",
#  "values": [
#   {
#    "x": 0,
#    "y": 0.14202716750712038
#   },
#   {
#    "x": 1,
#    "y": 0.17686738381528422
#   },
#   {
#    "x": 2,
#    "y": 0.2520439938176909
#   },


    #write a basic histogram to a file!
    def writeToFile(self,filename=None,defOneAsMin=True):
        if filename and filename != '-':
            output = open(os.path.abspath(filename), 'w')
        else:
            output = sys.stdout

        min_val=self.gminVal()
        max_val=self.gmaxVal()
        if (type(min_val)!=int) or (type(max_val)!=int):
            min_val=(-1)
            max_val=(-1)
            print "Warning, could not determine min/max values for CDR3 lengths!  Output may not be defined!"
        if defOneAsMin:
            min_val=1

        self.modes.sort()
        # write header
        output.write("CDR3_LENGTH")
        for mode in self.modes:
            output.write("\t"+mode)
        output.write("\n")

        for v in range(min(min_val,max_val),max(min_val,max_val)+1):
            if(v!=(0)):
                output.write(str(v))
                for mode in self.modes:
                    actual_val=0
                    if(v in self.count_map[mode]):
                        actual_val=self.count_map[mode][v]
                    output.write("\t"+str(actual_val))
                output.write("\n")
        if output is not sys.stdout:
            output.close()



#cache for CDR3 start/end
rsmap=dict()
remap=dict()
domain_list=get_domain_modes()
for d in domain_list:
       rsmap[d]=dict()
       remap[d]=dict()




def getEmptyCDR3Map():
    domain_modes=get_domain_modes()
    cdr3_hist=dict()
    for dm in domain_modes:
        cdr3_hist[dm+'_from']=(-1)
        cdr3_hist[dm+'_to']=(-1)
        cdr3_hist[dm]=(-1)
    return cdr3_hist
        
        







def VJRearrangementInFrameTest(vInfo,jInfo,imgtdb_obj,organism):
    if(vInfo==None or jInfo==None):
        #need valid data to test. return false in this case
        return False
    #use V and J frame
    if(jInfo==None):
        #no J means no prod. rearrangment???
        return False
    s_end=int(vInfo['s. end'])
    q_end=int(vInfo['q. end'])
    refName=vInfo['subject ids']
    s_end_frame=getTheFrameForThisReferenceAtThisPosition(refName,organism,imgtdb_obj,s_end)
    q_end_frame=s_end_frame
    #print "========================================================="
    #print "\n\n\nLooking at ",vInfo['query id']," for FRAME test...."
    #print "q_start_frame (spos=",s_end,") ",q_end_frame
    q_end_j=int(jInfo['q. end'])
    q_bgn_j=int(jInfo['q. start'])
    #print "First test...."
    if(q_end_j<=q_end):
        #print "WAY TOO SHORT or BAD ALIGNMENT"
        #print "IF Q END IN J IS LESS THAN OR EQUAL TO Q END IN V\n\n\n"
        return False
    else:
        #print "Second test"
        j_bgn_j=int(jInfo['s. start'])
        j_bgn_frame=getTheFrameForThisJReferenceAtThisPosition(jInfo['subject ids'],organism,imgtdb_obj,j_bgn_j)
        #print "v_end_frame is ",s_end_frame
        #print "j_bgn_frame is ",j_bgn_frame
        if(j_bgn_frame!=None):
            #print "passed into final...."
            #print "V sub and qry end are : ",s_end," and ",q_end
            #print "J sub and query bgn are : ",j_bgn_j," and ",q_bgn_j
            expected_frame_based_on_V=(q_end_frame+(q_bgn_j-q_end))%3
            #print "expected_frame_based_on_V=",expected_frame_based_on_V,"at read position=",q_bgn_j
            expected_frame_based_on_J=j_bgn_frame
            #print "expected_frame_based_on_J=",expected_frame_based_on_J,"at read position=",q_bgn_j
            if(expected_frame_based_on_J==expected_frame_based_on_V):
                #if both V and J impose the same frame, then the read should NOT be filtered
                #print "V extrapolated frame (",expected_frame_based_on_V,") is the same as J frame so return TRUE for frame test."
                return True
            else:
                #print "V extrapolated frame (",expected_frame_based_on_V,") is NOT the same as J frame so return TRUE for frame test."
                return False
        return None





#given info maps for V and J and the return a 
#dictionary with kabat and imgt lengths
#return other related information as well
def CDR3LengthAnalysis(vMap,jMap,organism,imgtdb_obj,read_rec):
    #currentQueryName=str(currentQueryName.strip())
    currentV=vMap['subject ids']
    currentJ=jMap['subject ids']
    cdr3_hist=getEmptyCDR3Map()
    cdr3_hist['Missing CYS']=True
    cdr3_hist['Missing TRP/PHE']=True
    cdr3_hist['Out-of-frame junction']=None
    cdr3_hist['Out-of-frame CDR3']=None
    if(vMap['query id'].find("reversed|")==0):
        cdr3_hist['qry_rev']=True
    else:
        cdr3_hist['qry_rev']=False

    if(looksLikeAlleleStr(currentV) and looksLikeAlleleStr(currentJ)):
        #print "WE'RE IN BUSINESS!"

        domain_modes=get_domain_modes()
        for dm in domain_modes:
            cdr3_hist[dm]=(-1)
            cdr3_hist[dm+'_from']=(-1)
            cdr3_hist[dm+'_to']=(-1)
            if(alleleIsTR(currentV) and dm=="kabat"):
                #no kabat analysis for CDR3 TR!
                return cdr3_hist
            #print "processing in ",dm
            if(not currentV in rsmap[dm]):
                #print currentV,"not in lookup for dm=",dm
                ref_cdr3_start=getAdjustedCDR3StartFromRefDirSetAllele(currentV,imgtdb_obj,organism,dm)
                rsmap[dm][currentV]=ref_cdr3_start
            else:
                ref_cdr3_start=rsmap[dm][currentV]
            #print "After initial retrieval, the ref CDR3 start for ",currentV," is ",ref_cdr3_start
            if(not currentJ in remap[dm]):
                #print currentJ,"not in lookup for dm=",dm
                ref_cdr3_end=getADJCDR3EndFromJAllele(currentJ,imgtdb_obj,organism,dm)
                remap[dm][currentJ]=ref_cdr3_end
            else:
                ref_cdr3_end=remap[dm][currentJ]
            if(ref_cdr3_start!=(-1) and ref_cdr3_end!=(-1)):
                #print "ref_cdr3_end is ",ref_cdr3_end
                if(dm=="imgt"):
                    cdr3_hist['Out-of-frame junction']=not(VJRearrangementInFrameTest(vMap,jMap,imgtdb_obj,organism))
                vq_aln=vMap['query seq']
                vs_aln=vMap['subject seq']
                vq_f=int(vMap['q. start'])
                vq_t=int(vMap['q. end'])
                vs_f=int(vMap['s. start'])
                vs_t=int(vMap['s. end'])
                jq_aln=jMap['query seq']
                js_aln=jMap['subject seq']
                jq_f=int(jMap['q. start'])
                jq_t=int(jMap['q. end'])
                js_f=int(jMap['s. start'])
                js_t=int(jMap['s. end'])
                
                # Get the start of CDR3
                ref_cdr3_start-=1
                qry_cdr3_start=getQueryIndexGivenSubjectIndexAndAlignment(vq_aln,vs_aln,vq_f,vq_t,vs_f,vs_t,ref_cdr3_start)
                #print "For ",currentV,dm," the ref CDR3 start is ",ref_cdr3_start
                #print "For ",currentV,dm," the qry CDR3 start is ",qry_cdr3_start

                # Check for Cys
                if (qry_cdr3_start != (-1) and dm=='imgt'):
                    # There are scenarios (issue #24) whereby the IgBlast alignment ends early
                    # and does not contain the sequence with the Cys.
                    # Thus we need to use the read sequence versus the alignment sequence.
                    qw=str(read_rec.seq)
                    if(cdr3_hist['qry_rev']):
                        qw=rev_comp_dna(qw)
                    #print "read_rec=",qw

                    qry_cys_start = qry_cdr3_start-3
                    qry_cys_end = qry_cdr3_start
                    qry_test_cys=qw[qry_cys_start:qry_cys_end]
                    #print "qry_cys_start",qry_cys_start
                    #print "qry_cys_end",qry_cys_end
                    #print "qry_test_cys=",qry_test_cys
                    if (len(qry_test_cys) == 3):
                        if (codonAnalyzer.is_unambiguous_codon(qry_test_cys) == True):
                            query_cys_trx=codonAnalyzer.fastTrans(qry_test_cys)
                            #print "The fast trans (with query=",vMap['query id']," and ref=",currentV,")=",query_cys_trx
                            if(query_cys_trx=='C'):
                                cdr3_hist['Missing CYS']=False
                            else:
                                cdr3_hist['Missing CYS']=True

                # CDR3 starts 1bp after alignment
                if(qry_cdr3_start!=(-1)):
                    qry_cdr3_start+=1

                # determine CDR3 end codon
                locus=currentV[0:4]
                cdr3_end_codon=get_cdr3_end(locus)
                #print "locus=",locus,"codon=",cdr3_end_codon

                # Check for TRP
                ref_cdr3_trp_start=ref_cdr3_end+1 #cause CDR3 end is 1bp before the TRP start
                #print "ref_cdr3_trp_start=",ref_cdr3_trp_start
                qry_trp_start=getQueryIndexGivenSubjectIndexAndAlignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t,ref_cdr3_trp_start,"left")
                #print "qry_trp_start=",qry_trp_start
                if (qry_trp_start!=(-1) and dm=='imgt'):
                    # use query sequence because alignment may not contain Trp/Phe
                    qw=str(read_rec.seq)
                    if(cdr3_hist['qry_rev']):
                        qw=rev_comp_dna(qw)
                    #print "read_rec=",qw

                    qry_trp_end = qry_trp_start+2
                    qry_test_trp=qw[qry_trp_start-1:qry_trp_end]
                    #print "qry_trp_start",qry_trp_start
                    #print "qry_trp_end",qry_trp_end
                    #print "qry_test_trp=",qry_test_trp
                    if (len(qry_test_trp) == 3):
                        if (codonAnalyzer.is_unambiguous_codon(qry_test_trp) == True):
                            qry_trp_trx=codonAnalyzer.fastTrans(qry_test_trp)
                            #print "The fast trans (with query=",vMap['query id']," and ref=",currentV,")=",qry_trp_trx
                            if(qry_trp_trx==cdr3_end_codon):
                                cdr3_hist['Missing TRP/PHE']=False
                            else:
                                cdr3_hist['Missing TRP/PHE']=True

                # CDR3 end
                ref_cdr3_end+=1
                qry_cdr3_end=getQueryIndexGivenSubjectIndexAndAlignment(jq_aln,js_aln,jq_f,jq_t,js_f,js_t,ref_cdr3_end,"left")
                #print "For ",currentV,dm," the ref CDR3 end is ",ref_cdr3_end
                #print "For ",currentV,dm," the qry CDR3 end is ",qry_cdr3_end
                if(qry_cdr3_end!=(-1)):
                    qry_cdr3_end-=1

                # Verify CDR3 start and end before reporting
                if(qry_cdr3_start!=(-1) and qry_cdr3_end!=(-1)):
                    if(dm=="imgt"):
                        if(((qry_cdr3_end-qry_cdr3_start+1)%3)==0):
                            cdr3_hist['Out-of-frame CDR3']=False
                        else:
                            cdr3_hist['Out-of-frame CDR3']=True
                    #query_coding_seq=query_seq_map[currentQueryName]
                    #coding_seq=query_coding_seq[(qry_cdr3_start-1):(qry_cdr3_end-1)]
                    if(qry_cdr3_start<=qry_cdr3_end):
                        #good
                        cdr3_len=qry_cdr3_end-qry_cdr3_start+1
                        #print "mode=",dm,"read=",vMap['query id']
                        #print "For ",currentV," the ref_cdr3_bgn is ",ref_cdr3_start,"!"
                        #print "For ",currentJ," the ref_cdr3_end is ",ref_cdr3_end,"!"
                        #print "read start and end are ",qry_cdr3_start," and ",qry_cdr3_end,"\n\n\n\n"
                        cdr3_hist[dm]=cdr3_len
                        cdr3_hist[dm+'_from']=qry_cdr3_start
                        cdr3_hist[dm+'_to']=qry_cdr3_end
                    else:
                        #messed up alignment presumably due to overlap! or in the RARE case where J aligns before V in the read!
                        #print "messed up alignment presumably due to overlap! or in the RARE case where J aligns before V in the read!"
                        pass
                else:
                    #print "BADQRYMAP Failure to map to query for mode=",dm," V=",currentV," J=",currentJ," read=",vMap['query id'],"  REFSTART=",ref_cdr3_start,"QRYSTART=",qry_cdr3_start,"REFEND=",ref_cdr3_end,"QRYEND=",qry_cdr3_end
                    pass
            else:
                #print "BADREFMAP mode=",dm," VALLELE=",currentV," and JALLELE=",currentJ," refVCDR3=(-1) for ",currentV," = ",ref_cdr3_start," or refJCDR3 ",currentJ," = ",ref_cdr3_end
                pass
    else:
        print "Ref names ",currentV," and ",currentJ," don't appear alleleic!"
        pass
    #print "***************************\n"
    #print "RETURNING THIS CDR3_HIST for ",vMap['query id'],":"
    #printMap(cdr3_hist)
    #print "***************************\n"  
    return cdr3_hist


#program to merge CDR3 histograms for kabat and imgt modes
if (__name__=="__main__"):
    parser=argparse.ArgumentParser(description='Merge multiple CDR3 length histograms into a single histogram.  Write the merged result to stdout')
    parser.add_argument('-i', dest='infiles', type=str, required=True, nargs='+', help="path(s) to CDR3 histograms to merge.  At least one is requried!")
    parser.add_argument('-o', dest='outfile', type=str, required=False, help="Write the output as a 'stream-layer' formatted JSON string")
    parser.add_argument('-j', dest='j', action='store_true', help="Write the output as a 'stream-layer' formatted JSON string")
    args=parser.parse_args()

    main_hist=histoMapClass(get_domain_modes())
    for infile in args.infiles:
        temp_hist=histoMapClass(get_domain_modes())
        temp_hist.read_from_file(infile)
        main_hist.merge(temp_hist)
    if args.j: # print JSON
        print main_hist.JSONIFYStreamAsString()
    else:
        if args.outfile:
            main_hist.writeToFile(args.outfile,defOneAsMin=False)
        else:
            main_hist.writeToFile(defOneAsMin=False) #sys.stdout

