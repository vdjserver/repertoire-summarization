#!/usr/bin/env python

import fileinput
from codon_analysis import codonCounter

class ags_manager:
	"""
	A class that is 'fed' numbered mutations (e.g. A31BC)
	keeps track of counts, then, when requested computes
	AGS and NMO scores
	"""
	ags6_count=0
	ags5_count=0
	rm_count=0
	ags6_nums=None
	ags5_nums=None
	name=None
	def __init__(self,init_name):
		self.name=init_name
		#self.codonCounterObj=codonCounter(init_name)
		self.ags6_nums=["31B","40","56","57","81","89"]
		self.ags5_nums=["31B","40","56","57","81"]
	
	def receive_numbered_mut(self,numbered_mut):
		"""add a mutation to count"""
		mut_pieces=self.appearsToBeNumberedMut(numbered_mut)
		if(mut_pieces!=None):
			numbered_pos=mut_pieces[1]
			if(numbered_pos in self.ags6_nums):
				self.ags6_count+=1
			if(numbered_pos in self.ags5_nums):
				self.ags5_count+=1
			self.rm_count+=1
	
	def appearsToBeNumberedMut(self,nm):
		"""
		return False if it's not a numbered mutation, 
		but a 3-element array of the FROM, POSITION, and TO 
		information of the mutation if it is		
		"""
		#return self.codonCounterObj.appearsToBeNumberedMut(nm)
		import re
		nmRe=re.compile(r'^([A-Z\*])([0-9][0-9][ABC]?)([A-Z\*])$')
		nmmo=re.search(nmRe,nm)
		if(nmmo):
			mut_from=nmmo.group(1)
			mut_pos=nmmo.group(2)
			mut_to=nmmo.group(3)
			return [str(mut_from),str(mut_pos),str(mut_to)]
		else:
			return False
	
	
	
	def compute_ags(self,mode="AGS6",returnThreeZeroesIfDivByZero=False):
		"""
		Given the mutations whose countings have been
		accrued/aggregated, compute AGS scores based on them
		"""
		#print "mgr named ",self.name,"giving ags from counter obj named ",self.codonCounterObj.name
		if(mode=="AGS6"):
			sampTot=self.rm_count
			if(sampTot<=0):
				#avoid division by zero error
				if(returnThreeZeroesIfDivByZero):
					return [0,0,0]
				return ["NA",0,0]			
			sampAGSTot=float(self.ags6_count)
			pct_ags6=(sampAGSTot/sampTot)*100.0
			ags_numerator=pct_ags6-1.6*6.0
			ags_denominator=0.9
			ags6_score=ags_numerator/ags_denominator
			return [ags6_score,sampAGSTot,sampTot]
		elif(mode=="AGS5"):
			sampTot=self.rm_count
			if(sampTot<=0):
				#avoid division by zero error
				if(returnThreeZeroesIfDivByZero):
					return [0,0,0]
				return ["NA",0,0]	
			sampAGSTot=float(self.ags5_count)
			pct_ags5=(sampAGSTot/sampTot)*100.0
			ags_numerator=pct_ags5-1.6*5.0
			ags_denominator=0.9
			ags5_score=ags_numerator/ags_denominator
			return [ags5_score,sampAGSTot,sampTot]

		else:
			raise Exception("ERROR, UNHANDLED AGS MODE = "+mode)






if (__name__=="__main__"):
	#myCounter=codonCounter("/home/data/vdj_server/repertoire-summarization/codon_data/codon_pos_IGHV4")
	myAGSMgr=ags_manager("stdin")
	#print "Enter in RM data "
	for line in fileinput.input():
    		#print "got '"+line.strip()+"'"
		pm=line.strip()
		if(myAGSMgr.appearsToBeNumberedMut(pm)):
			#print "got a pm ",pm," ",myAGSMgr.appearsToBeNumberedMut(pm)
			myAGSMgr.receive_numbered_mut(pm)
		else:
			#print "not a pm"
			pass
	print "AGS5",myAGSMgr.compute_ags("AGS5")
	print "AGS6",myAGSMgr.compute_ags("AGS5")





