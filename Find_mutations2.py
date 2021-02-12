import scipy.stats as stats
import numpy as np
import sys
from scipy.stats import binom


def hist2vec(histogram):
	vec= histogram[0, 0] * np.ones(histogram[1, 0])
	for i in range(histogram.shape[1]-1):
		ve_temp= histogram[0, i + 1] * np.ones(histogram[1, i + 1])
		vec = np.concatenate((vec, ve_temp), axis=0)
	return vec


def log_likelihood(reads, A, f, P):
	L_k_log =0
	for k in range(int(reads.size/2)):
		L_k_log+=reads[1][k]*np.log(sum(f*P[A, reads[0][k]]))
	return L_k_log


def Check_Norm_valid(Norm_reads,Norm_allele,P ,p_equal):
	if Norm_allele.size>2:
		return -2
	elif Norm_allele.size==2:
		A1=Norm_allele[0]
		A2=Norm_allele[1]
		n1=Norm_reads[1][np.nonzero(Norm_reads[0][:]==A1)]
		n2=Norm_reads[1][np.nonzero(Norm_reads[0][:]==A2)]
		#print n1,n2,Norm_reads
		p=binom.cdf(min(n1,n2),n1+n2,0.5)
		#print p,p_equal
		if p<p_equal:
			return -1
		else:
			return 1
	else:
		return 1


def Check_Mutation(Norm_reads,Norm_allele,Norm_frac,Tum_reads,Tum_allele,Tum_frac,P,LOR_ratio,p_equal,KS_thresh):
	sN=Norm_allele.shape
	sN=sN[0]
	sT=Tum_allele.shape
	sT=sT[0]
	#print "Allel num",Norm_allele,sN,Tum_allele,sT,
	Thresh_hold=-LOR_ratio
	p=p_equal
	che=Check_Norm_valid(Norm_reads,Norm_allele,P,p)
	if che<1:
		return che
	else:
		#print ("Log data before sending",Norm_reads,Tum_allele,Tum_frac)
		L_Norm_Tum=log_likelihood(Norm_reads,Tum_allele,Tum_frac,P)
		L_Norm_Norm=log_likelihood(Norm_reads,Norm_allele,Norm_frac,P)
		L_Tum_Tum=log_likelihood(Tum_reads,Tum_allele,Tum_frac,P)
		L_Tum_Norm=log_likelihood(Tum_reads,Norm_allele,Norm_frac,P)

		AIC_Norm_Tum =2*sT-2*L_Norm_Tum
		AIC_Norm_Norm=2*sN-2*L_Norm_Norm
		AIC_Tum_Tum  =2*sT-2*L_Tum_Tum
		AIC_Tum_Norm =2*sN-2*L_Tum_Norm

		if AIC_Tum_Tum - AIC_Tum_Norm < Thresh_hold and AIC_Norm_Norm - AIC_Norm_Tum < Thresh_hold:

			#Checking the KS test
			vec_N=hist2vec(Norm_reads)
			vec_T=hist2vec(Tum_reads)
			[ks_vla,ks_p]=stats.ks_2samp(vec_N,vec_T)
			if (ks_p<KS_thresh):
				return 1
			else:
				return -3
		else:
			return 0


def main(cmd_args):
	probability_table = np.loadtxt(cmd_args[3], delimiter=',')
	tumor_file = open(cmd_args[1], "r")
	normal_file = open(cmd_args[2], "r")
	LOR_ratio = float(cmd_args[4])
	p_equal = float(cmd_args[5])
	KS_thresh = float(cmd_args[6])

	count, count_mut, c1 = 0, 0, 0
	for Norm_t in normal_file:
		Norm_t = Norm_t[0:len(Norm_t)-2]
		Norm_list=Norm_t.split(" -99999999 ")
		Norm_list[3].strip('\n')
		Tum_t = tumor_file.readline()
		Tum_list = Tum_t.split(" -99999999 ")
		chr_l = Norm_t.split("ZZZ")
		locus = chr_l[0]

		Norm_reads_list = Norm_list[0].split(" ")
		Norm_reads = np.empty(len(Norm_reads_list)-1, dtype=int)
		for i in range(1,len(Norm_reads_list)):
			Norm_reads[i-1]=int(Norm_reads_list[i])
		Norm_reads_mat = np.array([Norm_reads[0:Norm_reads.size:2],Norm_reads[1:Norm_reads.size:2]])
		fi=np.nonzero(Norm_reads_mat[0,:]<40)
		Norm_reads_mat=Norm_reads_mat[:,fi[0]]
		c1+=1

		Norm_allele = np.fromstring(Norm_list[2], dtype=int, sep=' ')

		Norm_frac= np.fromstring(Norm_list[3], dtype=float, sep=' ')
		if Norm_allele.size>1:
			Norm_allele_arg=Norm_allele.argsort()
			Norm_frac=Norm_frac[Norm_allele_arg]
			Norm_allele=Norm_allele[Norm_allele_arg]
			#print "Nprmal",Norm_allele,Norm_frac

		#print Tum_list[0]
		Tum_reads_list = Tum_list[0].split(" ")
		Tum_reads = np.empty(len(Tum_reads_list)-1, dtype=int)
		#print Tum_reads_list
		for i in range(1,len(Tum_reads_list)):
			Tum_reads[i-1]=int(Tum_reads_list[i])
		Tum_reads_mat=np.array([Tum_reads[0:Tum_reads.size:2],Tum_reads[1:Tum_reads.size:2]])
		fi=np.nonzero(Tum_reads_mat[0,:]<40)
		Tum_reads_mat=Tum_reads_mat[:,fi[0]]


		Tum_allele = np.fromstring(Tum_list[2], dtype=int, sep=' ')
		#print Tum_allele,"size is",Tum_allele.size
		Tum_frac= np.fromstring(Tum_list[3], dtype=float, sep=' ')
		#print "Before if",Norm_allele.shape," ", Tum_allele.shape
		if Tum_allele.size>1:
			Tum_allele_arg=Tum_allele.argsort()
			Tum_frac=Tum_frac[Tum_allele_arg]
			Tum_allele=Tum_allele[Tum_allele_arg]


		if np.array_equal(Norm_allele,Tum_allele):
			c=1
		else:
			count+=1
			q=Check_Mutation(Norm_reads_mat, Norm_allele, Norm_frac, Tum_reads_mat, Tum_allele, Tum_frac, probability_table, LOR_ratio, p_equal, KS_thresh)
			if q==1:
				count_mut+=1
			print (q,locus,Norm_reads,Norm_allele,Norm_frac,Tum_reads,Tum_allele,Tum_frac,"@")

	print ("All",count,"Mut",count_mut)


if __name__== '__main__':
	main(sys.argv)