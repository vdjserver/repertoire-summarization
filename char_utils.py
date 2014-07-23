#!/usr/bin/env python


import re

#find the number '-' in a btop
def getNumberIndelsFromBTOP(btop):
	num_indels=0
	for btop_char in btop:
		if(btop_char=="-"):
			num_indels+=1
	return num_indels




#given a position in inverted space and the sequence length in base pairs
#return the position in regular (un-inverted) space
def getRegPosFromInvertedPos(inv_pos, seq_len_in_bp):
    '''
    from a position in the inverted space, get the position in regular space
    '''
    fake_interval = [inv_pos, inv_pos]
    reg_fake_interval = getRevCompInterval(fake_interval, seq_len_in_bp)
    return reg_fake_interval[0]





#find the number of base substitutions (does NOT include indels) from a BTOP
def getNumberBaseSubsFromBTOP(btop):
	#remove all digits
	btop_orig=str(btop)
	btop=re.sub(r'[0-9]','',btop)
	#all that remains is indels and mutations (if there are any)
	if(len(btop)==0):
		return 0
	else:
		#found at least one mutation or indel
		#btops for these are pairs of S/Q data
		if(not(len(btop)%2==0)):
			#ummm....they should be in pairs!
			print "btop should be even for mutations and indels!"
			sys.exit(0)
		number_muts_and_indels=len(btop)/2
		num_indels=getNumberIndelsFromBTOP(btop_orig)
		num_muts=number_muts_and_indels-num_indels
		return num_muts


#given a btop, find the number of insertions in the query,
# and the number of deletions in the query
def getIndelMapFromBTOP(btop):
	#q=CAGTAGTAACTGGTGGAGTTGGGTCCGCCAGCCCCCAGG-AAGGGGCTGGAGTGGATTGGGGAAATCTATCATAGTGGGAGCACCAACTACAATCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACACGGCCAAGAACCAGATCTCCCTGAACCTGACCTCTGTGACCGCTGCGGACACGGCCGTGTATTTTTGTGCGAGAG
	#s=CAGTAGTAACTGGTGGAGTTGGGTCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGGAAATCTATCATAGTGGGAGCACCAACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCAGTAGACAAGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCGGACACGGCCGTGTATTACTGTGCGAGAG
	#btop=39-G53TC37CA1GT11AT10CG4CG12TC19TATC10
	#in the pairings of BTOP, the query character is first followed by the subject character
	indel_map=dict()
	indel_map['insertions']=0
	indel_map['deletions']=0
	btop=re.sub(r'[0-9]','',btop)
	while(len(btop)>2):
		q_char=btop[0]
		s_char=btop[1]
		if(q_char=="-"):
			indel_map['deletions']+=1
		if(s_char=="-"):
			indel_map['insertions']+=1
		btop=btop[2:]
	return indel_map
