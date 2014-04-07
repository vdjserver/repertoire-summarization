#!/usr/bin/env python

import multiprocessing
import threading
from vdjml_utils import getVDJServerRegionAlignmentFromLargerVAlignmentPyObj
from alignment import alignment

#thread class for parallel charactgerization of regions
class CharacterizationThread(threading.Thread):
	def clear_result(self):
		self.result=None
	def get_result(self):
		return self.result
	def __init__(self,queue):
		#print "IN INIT OF CharacterizationThread"
		threading.Thread.__init__(self)
		self.queue=queue
	def run(self):
		while(True):
			try:
				char_job_dict=self.queue.get()
				read_result_obj=char_job_dict['read_result_obj']
				meta=char_job_dict['meta']
				organism=char_job_dict['organism']
				mode=char_job_dict['mode']
				region=char_job_dict['region']
				imgtdb_obj=char_job_dict['imgtdb_obj']
				read_rec=char_job_dict['read_rec']
				key_base=char_job_dict['key_base']
				region=char_job_dict['region']
				noneSeg_flag=char_job_dict['noneSeg_flag']
				if(not(noneSeg_flag)):
					regionAlignment=getVDJServerRegionAlignmentFromLargerVAlignmentPyObj(read_result_obj,meta,organism,mode,region,imgtdb_obj,False,read_rec)
					if(regionAlignment!=None):
						reg_ann_show_msg=False
						#reg_ann_msg="characterization for region="+region+" mode="+mode+" for read="+read_rec.id
						reg_ann_msg=None
						reg_ann=regionAlignment.characterize(reg_ann_msg,reg_ann_show_msg)
						#if(mode=="imgt" and read_rec.id=="FR3_STOP" and region=="FR3"):
						#	print "got ann "
						#	printMap(reg_ann)
						#	print "The alignment nice : "
						#	print regionAlignment.getNiceString()
						#	sys.exit(0)
						annMap=dict()
						for key in reg_ann:
							annMap[key_base+region+"_"+key]=reg_ann[key]
						annMap[key_base+region+'_qry_aln']=regionAlignment.q_aln
						annMap[key_base+region+'_qry_srt']=regionAlignment.q_start
						annMap[key_base+region+'_qry_end']=regionAlignment.q_end
						annMap[key_base+region+'_ref_aln']=regionAlignment.s_aln
						annMap[key_base+region+'_frm_msk']=regionAlignment.s_frame_mask
						self.result=annMap
				self.queue.task_done()
			except:
				import traceback,sys
				print "Exception in user code:"
				print '-'*60
				traceback.print_exc(file=sys.stdout)
				print '-'*60
				import os
				os._exit(1)
			

	

	



