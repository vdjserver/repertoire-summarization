#!/usr/bin/env python


from utils import readFileIntoArrayOneLinePerArrayElement




#read a file of CDR3s into an array
def readCDR3FileOnePerLineAsCDR3List(path,filterOutStop=True)
	import re
	cdr3_list=readFileIntoArrayOneLinePerArrayElement(path)
	if(filter_out_stop):
		new_list=list()
		for c in cdr3_list:
			if(re.search("\*",c)):
				#it has a stop, don't add it!
			else:
				new_list.append(c)
		cdr3_list=new_list
	return cdr3_list
	






#given an AA get the 10-d KIDERA vector
def getFactorsGivenResidue(aa):
	if(aa=="A"):
		ALA=[-1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.20,-0.08,0.21,-0.48]
		fact=AA
	elif(aa=="D"):
		ASP=[0.58,-0.22,-1.58,0.81,-0.92,0.15,-1.52,0.47,0.76,0.70]
		fact=ASP
	elif(aa=="C"):
		CYS=[0.12,-0.89,0.45,-1.05,-0.71,2.41,1.52,-0.69,1.13,1.10]
		fact=CYS
	elif(aa=="E"):
		GLU=[1.45,0.19,-1.61,1.17,-1.31,0.40,0.04,0.38,-0.35,-0.12]
		fact=GLU
	elif(aa=="F"):
		PHE=[-0.21,0.98,-0.36,-1.43,0.22,-0.81,0.67,1.10,1.71,-0.44]
		fact=PHE
	elif(aa=="G"):
		GLY=[1.46,-1.96,-0.23,-0.16,0.10,-0.11,1.32,2.36,-1.66,0.46]
		fact=GLY
	elif(aa=="H"):
		HIS=[-0.41,0.52,-0.28,0.28,1.61,1.01,-1.85,0.47,1.13,1.63]
		fact=HIS
	elif(aa="I"):
		ILE=[-0.73,-0.16,1.79,-0.77,-0.54,0.03,-0.83,0.51,0.66,-1.78]
		fact=ILE
	elif(aa="K"):
		LYS=[-0.34,0.82,-0.23,1.70,1.54,-1.62,1.15,-0.08,-0.48,0.60]
		fact=LYS
	elif(aa=="L"):
		LEU=[-1.04,0.00,-0.24,-1.10,-0.55,-2.05,0.96,-0.76,0.45,0.93]
		fact=LEU
	elif(aa=="M"):
		MET=[-1.40,0.18,-0.42,-0.73,2.00,1.52,0.26,0.11,-1.27,0.27]
		fact=MET
	elif(aa=="N"):
		ASN=[1.14,-0.07,-0.12,0.81,0.18,0.37,-0.09,1.23,1.10,-1.73]
		fact=ASN
	elif(aa=="P"):
		PRO=[2.06,-0.33,-1.15,-0.75,0.88,-0.45,0.30,-2.30,0.74,-0.28]
		fact=PRO
	elif(aa=="Q"):
		GLN=[-0.47,0.24,0.07,1.10,1.10,0.59,0.84,-0.71,-0.03,-2.33]
		fact=GLN
	elif(aa=="R"):
		ARG=[0.22,1.27,1.37,1.87,-1.70,0.46,0.92,-0.39,0.23,0.93]
		fact=ARG
	elif(aa=="S"):
		SER=[0.81,-1.08,0.16,0.42,-0.21,-0.43,-1.89,-1.15,-0.97,-0.23]
		fact=SER
	elif(aa=="T"):
		THR=[0.26,-0.70,1.21,0.63,-0.10,0.21,0.24,-1.15,-0.56,0.19]
		fact=THR
	elif(aa=="V"):
		VAL=[-0.74,-0.71,2.04,-0.40,0.50,-0.81,-1.07,0.06,-0.46,0.65]
		fact=VAL	
	elif(aa=="W"):
		TRP=[0.30,2.10,-0.72,-1.57,-1.16,0.57,-0.48,-0.40,-2.30,-0.60]
		fact=TRP
	elif(aa=="Y"):
		TYR=[1.38,1.48,0.80,-0.56,0.00,-0.68,-0.31,1.03,-0.05,0.53]
		fact=TYR
	else:
		raise Exception("Unknown amino acid ",aa)
	return fact



#generate an AA->KIDERA map
def createAAKideraMap():
	AAs="AVLIPFWMGSTCYNQDEKRH"
	AA_kidera_map=dict()
	for aa in AAs:
		AA_kidera_map[aa]=getFactorsGivenResidue(aa)
	return AA_kidera_map








#given a CDR3 AA string, compute the average KIDERA
def computeAverageKideraFromCDR3(cdr3):
	#Ten-dimensional Kidera Factor representations [18] (Table S1) for
	#each CDR3 sequence were computed by taking the average score
	#for each Kidera factor across all amino acids within the CDR3
	#sequence. This resulted in an average score for each of the Kidera
	#Factors represented as a ten dimensional vector and inserted into a
	#csv file.
	#Epstein M, Barenco M, Klein N, Hubank M, Callard RE (2014) Revealing 
	#Individual Signatures of Human T Cell CDR3 Sequence Repertoires with 
	#Kidera Factors. PLoS ONE 9(1): e86986. doi:10.1371/journal.pone.0086986
	aa_k_m=createAAKideraMap()
	kidera=[0,0,0,0,0,0,0,0,0,0]
	num_aminos=0
	for aa in cdr3:
		#get the kidera for the AA
		aa_kidera=aa_k_m(aa)
		#add it to the aggregate kidera
		for d in range(len(kidera)):
			kidera[d]+=aa_kidera[d]
	if(num_aminos!=0):
		for d in range(len(kidera)):
			kidera[d]=kidera[d]/num_aminos
	return kidera
		




#given one or more kideras, return the averaged kidera
def computeAverageKideraFromKideraList(kidera_list):
	kidera_avg=[0,0,0,0,0,0,0,0,0,0]
	for kidera in kidera_list:
		for d in range(10):
			kidera_avg[d]+=kidera[d]
	for d in range(10):
		kidera_avg[d]=float(kidera_avg[d])/float(len(kidera_list))
	return kidera_avg
	




#given two lists of kidera factors, return a distance matrix
def computeEuclideanDistMatrix(kidera_1,kidera_2):
	#init with kidera_1 as first index, then kidera_2 as second index
	dist_matrix = [[None for x in xrange(len(kidera_2))] for x in xrange(len(kidera_1))]
	for k1 in range(len(kidera_1)):
		for k2 in range(len(kidera_2)):
			dist=computeEuclidean(kidera_1[k1],kidera_2[k2])
			dist_matrix[k1][k2]=dist
	return dist_matrix





#given a pairwise distance matrix, convert it into a rank matrix
def constructRankMatrixGivenDistMatrix(dist_matrix):
	import scipy.stats as ss
	data_list=list()
	for r in range(len(dist_matrix)):
		for c in range(len(dist_matrix)):
			data_list.append(dist_matrix[r][c])
	data_list_ranking=ss.rankdata(data_list)
	data_to_ranking_dict=dict()
	for i in range(len(data_list)):
		data_to_ranking_dict[data_list[i]]=data_list_ranking[i]
	



#compute euclidean distance between two vectors
def computeEuclidean(v1,v2):
	#sum of squares of differences
	sos=0
	for d in v1:
		sos=(abs(v1[d]-v2[d))**2
	#ed gets computed as the square root of sos
	ed=sos**(0.5)
	return ed

#given a group assignment array compute the number of possible permutations
#it's a factorial
def computeNumAsnCombos(group_asn_arr):
	num_combos=math.factorial(len(group_asn_arr))
	return num_combos


#given a label assignment and rank 
#stat matrix, compute the R statistic!
def computeRStat(rank_sim_mat,group_asn_arr):
	in_grp_sum=0
	bt_grp_sum=0
	#num_in_sum=0
	total_num_samps_under_consideration=len(group_asn_arr)
	num_bt_sum=0
	num_gr_sum=0
	for r in range(len(group_asn_arr)):
		for c in range(len(group_asn_arr)):
			if(r<c):
				c_class=group_asn_arr[c]
				r_class=group_asn_arr[r]
				#num_in_sum+=1
				if(r_class==c_class):
					in_grp_sum+=group_asn_arr[r][c]
					num_gr_sum+=1
				else:
					bt_grp_sum+=group_asn_arr[r][c]
					num_bt_sum+=1
	in_grp_avg=float(in_grp_sum)/float(num_gr_sum)
	bt_grp_avg=float(bt_grp_sum)/float(num_bt_sum)
	n=float(total_num_samps_under_consideration)
	M=(n*(n-1.0))/2.0
	numerator=bt_grp_avg-bt_grp_avg
	denominator=M/2.0
	r_stat=numerator/denominator
	return r_stat




#return a string roundned to nearest 100th
def getPctStr(num,tot_poss):
	pct=float(num)/float(tot_poss)
	pct*=100.0
	pct_str=str(rount(pct,2))
	pct_str=pct_str[0:5]
	return pct_str







tstnegD2=readCDR3FileOnePerLineAsCDR3List("/home/data/Mei/TST-_Delta2.txt.CDR3.aa")
tstpozD2=readCDR3FileOnePerLineAsCDR3List("/home/data/Mei/TST+_Delta2.txt.CDR3.aa")








