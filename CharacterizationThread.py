#!/usr/bin/env python

import multiprocessing
import threading
from vdjml_utils import getVDJServerRegionAlignmentFromLargerVAlignmentPyObj,getHitInfo
from alignment import alignment
from segment_utils import getTheFrameForThisReferenceAtThisPosition,getVRegionStartAndStopGivenRefData
import re
from utils import printMap

#thread class for parallel charactgerization of regions
class CharacterizationThread(threading.Thread):

	def clear_result(self):
		self.result=None


	def __init__(self,queue,result_queue,alignment_output_queue):
		#print "IN INIT OF CharacterizationThread"
		threading.Thread.__init__(self)
		self.queue=queue
		self.result_queue=result_queue
		self.alignment_output_queue=alignment_output_queue

	def run(self):
		while(True):
			try:
				#print "to call CHAR_QUEUE get"
				char_job_dict=self.queue.get()
				#print "called CHAR_QUEUE get"
				read_result_obj=char_job_dict['read_result_obj']
				meta=char_job_dict['meta']
				organism=char_job_dict['organism']
				mode=char_job_dict['mode']
				region=char_job_dict['region']
				refName=char_job_dict['refName']
				imgtdb_obj=char_job_dict['imgtdb_obj']
				read_rec=char_job_dict['read_rec']
				key_base=char_job_dict['key_base']
				region=char_job_dict['region']
				noneSeg_flag=char_job_dict['noneSeg_flag']
				readName=char_job_dict['read_name']
				regionAlignment=None
				if(not(noneSeg_flag)):
					regionAlignment=getVDJServerRegionAlignmentFromLargerVAlignmentPyObj(read_result_obj,meta,organism,mode,region,imgtdb_obj,False,read_rec)
					if(regionAlignment!=None):
						regionAlignment.setName(region+"_"+mode+"_"+readName)
						s_start=regionAlignment.s_start
						s_end=regionAlignment.s_end
						ref_region_interval=getVRegionStartAndStopGivenRefData(refName,organism,imgtdb_obj,region,mode)
						largerVInfo=getHitInfo(read_result_obj,meta,refName,None,imgtdb_obj,organism)
						partialFlag=None
						if(s_start<=ref_region_interval[0] and s_end>=ref_region_interval[1]):
							#full region!
							partialFlag=False
						else:
							#partial region
							partialFlag=True
						if(not(partialFlag)):
							regionAlignment.setCompleteRegion(True)
						firstFrame=getTheFrameForThisReferenceAtThisPosition(refName,organism,imgtdb_obj,s_start)
						regionAlignment.setSFM(firstFrame)
						reg_ann_show_msg=False
						reg_ann_msg="characterization for region="+region+" mode="+mode+" for read="+read_rec.id
						#reg_ann_msg=None
						reg_ann=regionAlignment.characterize(reg_ann_msg,reg_ann_show_msg)
						#if(mode=="imgt" and re.search("IS9DOBB01B0XCX",str(read_rec.id)) and region=="CDR2"):
						#	print "got ann for ",read_rec," for region=CDR2"
						#	printMap(reg_ann)
						#	print "The alignment nice : "
						#	print regionAlignment.getNiceString()
						#	#sys.exit(0)
						#	print "\n\n\n\n\n\n\n\n\n"
						annMap=dict()
						for key in reg_ann:
							#annMap[key_base+region+"_"+key]=reg_ann[key]
							annMap[region.upper()+" "+key+" ("+mode+")"]=reg_ann[key]
						annMap[key_base+region+'_qry_aln']=regionAlignment.q_aln
						annMap[key_base+region+'_qry_srt']=regionAlignment.q_start
						annMap[key_base+region+'_qry_end']=regionAlignment.q_end
						annMap[key_base+region+'_ref_aln']=regionAlignment.s_aln
						annMap[key_base+region+'_frm_msk']=regionAlignment.s_frame_mask
						annMap[key_base+region+'_ref_srt']=int(s_start)
						annMap[key_base+region+'_ref_end']=int(s_end)
						#annMap[key_base+region+'_partial']=partialFlag
						self.result=annMap
				#print "to call CHAR_QUEUE task_done..."
				self.result_queue.put(self.result)
				self.alignment_output_queue.put(regionAlignment)
				self.queue.task_done()
				#print "called CHAR_QUEUE task_done..."
			except:
				import traceback,sys
				print "Exception in user code:"
				print '-'*60
				traceback.print_exc(file=sys.stdout)
				print '-'*60
				import os
				os._exit(1)
			

	

	



