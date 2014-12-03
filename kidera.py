#!/usr/bin/env python


from utils import readFileIntoArrayOneLinePerArrayElement




#read a file of CDR3s into an array
def readCDR3FileOnePerLineAsCDR3List(path,filterOutStop=True,filterOutX=True):
	import re
	cdr3_list=readFileIntoArrayOneLinePerArrayElement(path)
	if(filterOutStop):
		new_list=list()
		for c in cdr3_list:
			if(re.search("\*",c)):
				#it has a stop, don't add it!
				pass
			elif(re.search("X",c)):
				#filter out the X amino
				pass
			else:
				#print "adding ",c," from ",path
				new_list.append(c)
		cdr3_list=new_list
	return cdr3_list
	






#given an AA get the 10-d KIDERA vector
def getFactorsGivenResidue(aa):
	if(aa=="A"):
		ALA=[-1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.20,-0.08,0.21,-0.48]
		fact=ALA
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
	elif(aa=="I"):
		ILE=[-0.73,-0.16,1.79,-0.77,-0.54,0.03,-0.83,0.51,0.66,-1.78]
		fact=ILE
	elif(aa=="K"):
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



#from an array of CDR3 strings, return a list of as many kidera vectors
#or their average (if desired, default NOT)
def convertListOfCDR3sToListOfKideras(cdr3_list,getAvg=False):
	kidera_list=list()
	for cdr3 in cdr3_list:
		kidera_list.append(computeAverageKideraFromCDR3(cdr3))
	if(getAvg):
		return computeAverageKideraFromKideraList(kidera_list)
	else:
		return kidera_list



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
		#get the 10D kidera for the AA
		aa_kidera=aa_k_m[aa]
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
	



#compute euclidean distance between two vectors
def computeEuclidean(v1,v2):
	#sum of squares of differences
	sos=0
	for d in range(len(v1)):
		sos+=(abs(v1[d]-v2[d]))**2
	#ed gets computed as the square root of sos
	ed=sos**(0.5)
	return ed






#given two lists of kidera factors, return a distance matrix
def computeEuclideanDistMatrix(kidera_1,kidera_2):
	#init with kidera_1 as first index, then kidera_2 as second index
	dist_matrix = [[None for x in xrange(len(kidera_2))] for x in xrange(len(kidera_1))]
	for k1 in range(len(kidera_1)):
		for k2 in range(len(kidera_2)):
			#dist=computeEuclidean(kidera_1[k1],kidera_2[k2])
			#dist_matrix[k1][k2]=dist
			if(k1<=k2):
				dist=computeEuclidean(kidera_1[k1],kidera_2[k2])
				dist_matrix[k1][k2]=dist
			else:
				dist_matrix[k1][k2]=dist_matrix[k2][k1]
	return dist_matrix





#given a pairwise distance matrix, convert it into a rank matrix
def constructRankMatrixGivenDistMatrix(dist_matrix):
	import scipy.stats as ss
	data_list=list()
	for r in range(len(dist_matrix)):
		for c in range(len(dist_matrix[r])):
			if(r<c):
				#ignore diagonals
				data_list.append(dist_matrix[r][c])
	data_list_ranking=ss.rankdata(data_list)
	rank_mat=[[None for x in xrange(len(dist_matrix[0]))] for x in xrange(len(dist_matrix))]
	temp_i=0
	for r in range(len(dist_matrix)):
		for c in range(len(dist_matrix[r])):
			if(r<c):
				rank_mat[r][c]=data_list_ranking[temp_i]
				temp_i+=1
	for r in range(len(dist_matrix)):
		for c in range(len(dist_matrix[r])):
			if(r>c):
				rank_mat[r][c]=rank_mat[c][r]

	return rank_mat




#given a group assignment array compute the number of possible permutations
#it's a factorial
def computeNumAsnCombos(group_asn_arr,isList=True):
	import math
	if(isList):
		num_combos=math.factorial(len(group_asn_arr))
		return num_combos
	else:
		return math.factorial(int(group_asn_arr))	








#return a string roundned to nearest 100th
def getPctStr(num,tot_poss):
	pct=float(num)/float(tot_poss)
	pct*=100.0
	pct_str=str(rount(pct,2))
	pct_str=pct_str[0:5]
	return pct_str




#given a list of items
#return a random permutation of the list of items
#return it as a list
def randomPermutation(m_list):
	import random
	import numpy as np
	p_order=np.random.permutation(len(m_list))
	rand_permutation=list()
	for i in range(len(p_order)):
		rand_permutation.append(m_list[p_order[i]])
	return rand_permutation
	

#given a list of items have a generator
#that returns a given number (default 10)
#of permutations
def randomPermuter(m_list,numPermutations=10):
	t_list=m_list
	for p in range(numPermutations):
		t_list=randomPermutation(t_list)
		yield t_list


#get the dimenstions of a matrix (rows,columns)
def getMatrixDim(m):
	num_rows=len(m)
	if(num_rows==0):
		return [0,0]
	else:
		num_cols=len(m[0])
		return [num_rows,num_cols]


#create a merged list of kideras
#create a merged list of labels
#return them as a pair
def joinKideraListsAndMakeLabels(l1,l2,lab1="first",lab2="second"):
	merged=list()
	labels=list()
	for k in l1:
		merged.append(k)
		labels.append(lab1)
	for k in l2:
		merged.append(k)
		labels.append(lab2)
	return [merged,labels]


#pretty print a matrix
def printMatrixNice(m):
	print "The dim of a matrix is ",getMatrixDim(m)
	for r in range(len(m)):
		for c in range(len(m[r])):
			if(c!=len(m[r])-1):
				print str(round(m[r][c],2))+" , ",
			else:
				print str(round(m[r][c],2))
			


def kidera_analyze(f1,f2,lab1,lab2,num_permuts,fs_are_files=True):
	#load the input data first
	if(fs_are_files):
		#read two files of CDR3 polypeptides (one per line)
		cdr31=readCDR3FileOnePerLineAsCDR3List(f1)
		cdr32=readCDR3FileOnePerLineAsCDR3List(f2)
	else:
		#just load the data passed in
		cdr31=f1
		cdr32=f2
	#get corresponding kidera lists (one per polypeptide string)
	k1=convertListOfCDR3sToListOfKideras(cdr31)
	k2=convertListOfCDR3sToListOfKideras(cdr32)
	#merge them for subsequent use with a matrix
	merged_and_labels=joinKideraListsAndMakeLabels(k1,k2,lab1,lab2)
	merged=merged_and_labels[0]
	labels=merged_and_labels[1]
	#compute a euclidian distance matrix
	merged_dist_matrix=computeEuclideanDistMatrix(merged,merged)
	#rank it
	merged_rank_matrix=constructRankMatrixGivenDistMatrix(merged_dist_matrix)
	#compute the r statistic
	r_stat=computeRStat(merged_rank_matrix,labels)
	r_stats=list()
	num_pass=0
	permID=1
	for pl in randomPermuter(labels,num_permuts):
		if(permID%1000==0):
			print "On permutation #",permID
		#perm_merged_and_labeled=permute_merged_and_labeled(merged,labels)
		#perm_merged=perm_merged_and_labeled[0]
		#perm_merged_dist_matrix=computeEuclideanDistMatrix(perm_merged,perm_merged)
		#perm_merged_rank_matrix=constructRankMatrixGivenDistMatrix(perm_merged_dist_matrix)
		#temp_r=computeRStat(perm_merged_rank_matrix,pl)
		temp_r=computeRStat(merged_rank_matrix,pl)
		if(r_stat<0):
			if(r_stat<temp_r):
				num_pass+=1
		elif(r_stat>0):
			if(r_stat>temp_r):
				num_pass+=1
		permID+=1
	p_val=1.0-(float(num_pass)/float(num_permuts))
	#return the empirical r_stat, the p-value, then the permuted r_stats upon which the p-value is based
	ret_package=[r_stat,p_val,r_stats]
	return ret_package



def fix_matrix(m):
	for r in range(len(m)):
		for c in range(len(m[r])):
			if(m[r][c]==(-1)):
				m[r][c]=m[c][r]
	return m



def permute_merged_and_labeled(merged,labels):
	import numpy as np
	p_order=np.random.permutation(len(labels))
	new_labels=list(labels)
	new_merged=list()
	for i in range(len(p_order)):
		new_merged.append(merged[p_order[i]])
	return [new_merged,new_labels]
	




#given a label assignment and rank 
#stat matrix, compute the R statistic!
def computeRStat(rank_sim_mat,group_asn_arr):
	in_grp_sum=0
	bt_grp_sum=0
	total_num_samps_under_consideration=len(group_asn_arr)
	num_bt_sum=0
	num_gr_sum=0
	for r in range(len(group_asn_arr)):
		for c in range(len(group_asn_arr)):
			if(r<c):
				c_class=group_asn_arr[c]
				r_class=group_asn_arr[r]
				if(r_class==c_class):
					in_grp_sum+=rank_sim_mat[r][c]
					num_gr_sum+=1
				else:
					bt_grp_sum+=rank_sim_mat[r][c]
					num_bt_sum+=1
	#print "in g num",in_grp_sum," denom ",num_gr_sum
	in_grp_avg=float(in_grp_sum)/float(num_gr_sum)
	#print "bt g num",bt_grp_sum," denom ",num_bt_sum
	bt_grp_avg=float(bt_grp_sum)/float(num_bt_sum)
	n=float(total_num_samps_under_consideration)
	M=(n*(n-1.0))/2.0
	numerator=bt_grp_avg-in_grp_avg
	denominator=M/2.0
	#print "num=",numerator,"denom=",denominator
	r_stat=float(numerator)/float(denominator)
	return r_stat








def test_r_stat_computing():
	dm1 = [[None for x in xrange(6)] for x in xrange(6)]
	#this is a RANKED matrix
	#from FIG S1 in Epstein M, Barenco M, Klein N, Hubank M, Callard RE (2014) Revealing Individual Signatures of Human T 
	#Cell CDR3 Sequence Repertoires with Kidera Factors. PLoS ONE 9(1): e86986. doi:10.1371/journal.pone.0086986
	dm1=[
		[0,-1,-1,-1,-1,-1],
		[2,0,-1,-1,-1,-1],
		[3,6,0,-1,-1,-1],
		[1,4,5,0,-1,-1],
		[7,9,10,8,0,-1],
		[12,13,14,15,16,0]
		]
	labs1=["l1","l1","l1","l2","l2","l2"]
	#print "pre_fix d1 is "
	#printMatrixNice(dm1)
	#this "fix_matrix" fixes the -1 values
	dm1=fix_matrix(dm1)
	#print "post_fix d1 is "
	#printMatrixNice(dm1)
	r_stat_1=computeRStat(dm1,labs1)
	#print "it's ",r_stat_1
	if(r_stat_1==0.0):
		#should be zero
		pass
	else:
		raise Exception("Error, test r_stat 1 not equal to zero!")



def kidera_dist_print(f1,f2,fs_are_files=True):
	#load the input data first
	if(fs_are_files):
		#read two files of CDR3 polypeptides (one per line)
		cdr31=readCDR3FileOnePerLineAsCDR3List(f1)
		cdr32=readCDR3FileOnePerLineAsCDR3List(f2)
	else:
		#just load the data passed in
		cdr31=f1
		cdr32=f2
	#get corresponding kidera lists (one per polypeptide string)
	k1=convertListOfCDR3sToListOfKideras(cdr31)
	k2=convertListOfCDR3sToListOfKideras(cdr32)
	for i in k1:
		print i
	for i in k2:
		print i
	print "num k1 is ",len(k1)
	print "num k2 is ",len(k2)
	#merge them for subsequent use with a matrix
	merged_and_labels=joinKideraListsAndMakeLabels(k1,k2,lab1,lab2)
	merged=merged_and_labels[0]
	labels=merged_and_labels[1]
	#compute a euclidian distance matrix
	merged_dist_matrix=computeEuclideanDistMatrix(merged,merged)
	for r in range(len(merged_dist_matrix)):
		#print merged_dist_matrix[i]
		for c in range(len(merged_dist_matrix[r])):
			if(c<len(merged_dist_matrix[r])-1):
				#print str(merged_dist_matrix[r][c])+",",
				pass
			else:
				#print str(merged_dist_matrix[r][c])
				pass




#test_r_stat_computing()
f1="/home/data/Mei/TST-_Delta2.txt.CDR3.aa"
f2="/home/data/Mei/TST+_Delta2.txt.CDR3.aa"
lab1="neg"
lab2="poz"
num_permuts=10000


#kidera_dist_print(f1,f2)
ret_package=kidera_analyze(f1,f2,lab1,lab2,num_permuts)
r_stat=ret_package[0]
p_val=ret_package[1]
print "The r_stat is ",r_stat
print "The p-value is ",p_val
print "num permutations is ",num_permuts



