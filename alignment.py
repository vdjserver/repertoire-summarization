#!/usr/bin/env python




#class for two-seq alignment methods/tools/data
class alignment:

	#internal data
	s_aln=""
	q_aln=""
	s_start=(-1)
	s_end=(-1)
	q_start=(-1)
	q_end=(-1)
	

	#constructor
	def __init__(self,Cq_aln,Cs_aln,Cq_start,Cq_end,Cs_start,Cq_end):
		self.s_aln=self.Cs_aln
		self.q_aln=Cq_aln
		self.s_start=Cs_start
		self.s_end=Cs_end
		self.q_start=Cq_start
		self.q_end=Cq_end





	#given an alignment, and a position of a sequence IN the alignment
	#return the porition of the alignment whose bp are <= or >= that position in the alignment
	def getAlnAtAndCond(self,a_pos,st="query",cond="geq"):
		if(not(st=="query")):
			st="subject"
		aln=["",""] #Q, then S
		if(not(cond=="leq")):
			cond="geq"
		s_pos=self.s_start
		q_pos=self.q_start
		temp=0
		min_s_pos=float("inf")
		max_s_pos=float("-inf")
		min_q_pos=float("inf")
		max_q_pos=float("-inf")		
		while(temp<len(q_aln)):
			used_flag=False
			if(cond=="leq"):
				if(st=="query"):
					if(a_pos>=q_pos):
						used_flag=True
						aln[0]+=q_aln[temp]
						aln[1]+=s_aln[temp]
				else:
					if(a_pos>=s_pos):
						aln[0]+=q_aln[temp]
						aln[1]+=s_aln[temp]
						used_flag=True
			else:
				if(st=="query"):
					if(a_pos<=q_pos):
						aln[0]+=q_aln[temp]
						aln[1]+=s_aln[temp]
						used_flag=True
				else:
					if(a_pos<=s_pos):
						aln[0]+=q_aln[temp]
						aln[1]+=s_aln[temp]				
						used_flag=True
			if(used_flag):
				if(s_pos<=min_s_pos):
					min_s_pos=s_pos
				if(s_pos>=max_s_pos):
					max_s_pos=s_pos
				if(q_pos>=max_q_pos):
					max_q_pos
				if(q_pos<=min_q_pos):
					min_q_pos=q_pos
			if(s_aln[temp]!="-"):
				s_pos+=1
			if(q_aln[temp]!="-"):
				q_pos+=1
			temp+=1
		sub_aln=alignment(aln[0],aln[1],min_q_pos,max_q_pos,min_s_pos,max_s_pos)
		#return aln #Q then S
		return sub_aln
