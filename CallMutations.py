import scipy.stats as stats
import numpy as np
import sys
import math
from scipy.stats import binom
from collections import namedtuple
from icecream import ic
from typing import List

AlleleSet = namedtuple('Allele', ['alleles', 'fractions'])


class SUPPORT_LEVEL:
	# basically an enum for read support levels
	TOO_MANY_ALLELES = -2
	INSUFFICIENT = -1
	SUFFICIENT = 1


def check_read_support(Norm_reads, normal_alleles, p_equal) -> int:
	if normal_alleles.size < 2:
		return SUPPORT_LEVEL.SUFFICIENT
	elif normal_alleles.size==2:
		first_allele = normal_alleles[0]
		second_allele = normal_alleles[1]
		# reads supporting first allele
		first_allele_reads = Norm_reads[1][np.nonzero(Norm_reads[0][:]==first_allele)]
		second_allele_reads = Norm_reads[1][np.nonzero(Norm_reads[0][:]==second_allele)]
		p = binom.cdf(min(first_allele_reads, second_allele_reads), first_allele_reads+second_allele_reads, 0.5)
		if p < p_equal:
			return SUPPORT_LEVEL.INSUFFICIENT
		else:
			return SUPPORT_LEVEL.SUFFICIENT
	else: # Norm_allele.size>2:
		return SUPPORT_LEVEL.TOO_MANY_ALLELES


def log_likelihood(reads, alleles, fractions, probability_table):
	L_k_log = 0
	for k in range(int(reads.size/2)):
		L_k_log+=reads[1][k]*np.log(sum(fractions*probability_table[alleles, reads[0][k]]))
		"""
		if math.isinf(L_k_log):
			ic(fractions)
			ic(alleles)
			ic(reads[0][k])
			ic(probability_table[alleles, reads[0][k]])
			ic(sum(fractions*probability_table[alleles, reads[0][k]]))
			sys.exit(1)
		"""
	return L_k_log


def hist2vec(histogram):
	vector = histogram[0, 0] * np.ones(histogram[1, 0])
	for i in range(histogram.shape[1]-1):
		ve_temp= histogram[0, i + 1] * np.ones(histogram[1, i + 1])
		vector = np.concatenate((vector, ve_temp), axis=0)
	return vector


def check_mutation(normal_reads, normal_alleles, tumor_reads, tumor_allele, probability_table, LOR_ratio, p_equal, KS_threshold):
	num_normal_alleles = normal_alleles.alleles.shape[0]
	num_tumor_alleles = tumor_allele.alleles.shape[0]
	threshold = -LOR_ratio
	read_support_level = check_read_support(normal_reads, normal_alleles.alleles, p_equal)
	if read_support_level != SUPPORT_LEVEL.SUFFICIENT:
		return read_support_level
	else:
		L_Norm_Tum = log_likelihood(normal_reads, tumor_allele.alleles, tumor_allele.fractions, probability_table)
		L_Norm_Norm =  log_likelihood(normal_reads, normal_alleles.alleles, normal_alleles.fractions, probability_table)
		L_Tum_Tum = log_likelihood(tumor_reads, tumor_allele.alleles, tumor_allele.fractions, probability_table)
		L_Tum_Norm = log_likelihood(tumor_reads, normal_alleles.alleles, normal_alleles.fractions, probability_table)

		AIC_Norm_Tum = 2*num_tumor_alleles-2*L_Norm_Tum
		AIC_Norm_Norm = 2*num_normal_alleles-2*L_Norm_Norm
		AIC_Tum_Tum = 2*num_tumor_alleles-2*L_Tum_Tum
		AIC_Tum_Norm = 2*num_normal_alleles-2*L_Tum_Norm

		if AIC_Tum_Tum - AIC_Tum_Norm < threshold and AIC_Norm_Norm - AIC_Norm_Tum < threshold:
			# KS test
			vec_N = hist2vec(normal_reads)
			vec_T = hist2vec(tumor_reads)
			_, ks_p = stats.ks_2samp(vec_N, vec_T)
			if (ks_p < KS_threshold):
				return 1
			else:
				return -3
		else:
			return 0


def split_line(line: str) -> List[str]:
	line = line[:-2].strip()
	line_split = line.split((" -99999999 "))
	return line_split

def strip_brackets(text: str):
	# returns given text with brackets stripped off
	return (text.replace("[", "")).replace("]", "")

def get_reads(split_line: List[str])-> np.array:
	reads_list = split_line[0].split(" ")
	reads = np.empty(len(reads_list) - 1, dtype=int)
	for i in range(1, len(reads_list)):
		reads[i - 1] = int(reads_list[i])
	reads_mat = np.array([reads[0:reads.size:2], reads[1:reads.size:2]])
	repeat_filter = np.nonzero(reads_mat[0, :] < 40)
	reads_mat = reads_mat[:, repeat_filter[0]]
	return reads_mat


def get_alleles(split_line: List[str]) -> AlleleSet:
	alleles = np.fromstring(strip_brackets(split_line[2]), dtype=int, sep=' ')
	allelic_fractions = np.fromstring(split_line[3], dtype=float, sep=' ')
	if alleles.size>1:
		allele_sorted = alleles.argsort()
		allelic_fractions=allelic_fractions[allele_sorted]
		alleles = alleles[allele_sorted]
	return AlleleSet(alleles=alleles, fractions=allelic_fractions)

def write_output(output_lines, output_file: str):
	with open(output_file, 'w') as results:
		results.write("\n".join(output_lines))
	return 1

def main(cmd_args):
	probability_table = np.loadtxt(cmd_args[3], delimiter=',')
	tumor_file = open(cmd_args[1], "r")
	normal_file = open(cmd_args[2], "r")
	LOR_ratio = float(cmd_args[4])
	p_equal = float(cmd_args[5])
	KS_thresh = float(cmd_args[6])

	output_lines = []
	total_allelic_differentiated, mutated_loci = 0, 0
	normal_file_lines = normal_file.readlines()
	for i in range(len(normal_file_lines)):
		norm_line = normal_file_lines[i]
		locus = norm_line.split("ZZZ")[0]
		Norm_list = split_line(norm_line)
		Tum_list = split_line(tumor_file.readline())
		Norm_reads_mat = get_reads(Norm_list)
		normal_alleles = get_alleles(Norm_list)
		tumor_reads_mat = get_reads(Tum_list)
		tumor_alleles = get_alleles(Tum_list)
		if not np.array_equal(normal_alleles.alleles, tumor_alleles.alleles):
			total_allelic_differentiated+=1
			decision = check_mutation(Norm_reads_mat, normal_alleles, tumor_reads_mat, tumor_alleles, probability_table, LOR_ratio, p_equal, KS_thresh)
			if decision==1:
				mutated_loci+=1
			output_lines.append(f"{decision} {locus} {Norm_reads_mat} {normal_alleles.alleles} {normal_alleles.fractions} {tumor_reads_mat} {tumor_alleles.alleles} {tumor_alleles.fractions} @")
	output_lines.append(f"Total Loci with Differing Called Alleles: {total_allelic_differentiated}\n Total Mutated Loci: {mutated_loci}")
	write_output(output_lines, cmd_args(7))

if __name__== '__main__':
	main(sys.argv)
