#!/usr/bin/env python

import fileinput
from codon_analysis import codonCounter

class ags_manager:

	codonCounterObj=None
	ags6_nums=None
	ags5_nums=None
	name=None

	def __init__(self,init_name):
		self.name=init_name
		self.codonCounterObj=codonCounter(init_name)
		self.ags6_nums=self.codonCounterObj.get_ags6_nums()
		self.ags5_nums=self.codonCounterObj.get_ags5_nums()
		
		

	def receive_numbered_mut(self,numbered_mut):
		mut_pieces=self.appearsToBeNumberedMut(numbered_mut)
		if(mut_pieces!=None):
			numbered_pos=mut_pieces[1]
			if(numbered_pos in self.ags6_nums):
				self.codonCounterObj.increment_ags6_rep_muts(numbered_pos)
			if(numbered_pos in self.ags5_nums):
				self.codonCounterObj.increment_ags5_rep_muts(numbered_pos)			
			self.codonCounterObj.sampleRepMuts.increment(numbered_pos)

	
	def appearsToBeNumberedMut(self,nm):
		return self.codonCounterObj.appearsToBeNumberedMut(nm)


	def compute_ags(self,mode="AGS6"):
		#print "mgr named ",self.name,"giving ags from counter obj named ",self.codonCounterObj.name
		if(mode=="AGS6"):
			ags6_score=self.codonCounterObj.computeAGS()
			sampAGSTot=float(self.codonCounterObj.computeAGS6TotRM())
			sampTot=float(self.codonCounterObj.computeSampTotRM())
			return [ags6_score,sampAGSTot,sampTot]
		elif(mode=="AGS5"):
			ags5_score=self.codonCounterObj.computeAGS5()
			sampAGSTot=float(self.codonCounterObj.computeAGS5TotRM())
			sampTot=float(self.codonCounterObj.computeSampTotRM())
			return [ags5_score,sampAGSTot,sampTot]
		else:
			raise Exception("ERROR, UNHANDLED AGS MODE = "+mode)






if (__name__=="__main__"):
	#myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")
	myAGSMgr=ags_manager()
	print "now to enter...."
	for line in fileinput.input():
    		#print "got '"+line.strip()+"'"
		pm=line.strip()
		if(myAGSMgr.appearsToBeNumberedMut(pm)):
			#print "got a pm ",pm," ",myAGSMgr.appearsToBeNumberedMut(pm)
			myAGSMgr.receive_numbered_mut(pm)
		else:
			#print "not a pm"
			pass
	print myAGSMgr.compute_ags()






