import sys
import scipy.stats as stats
import numpy as np
from collections import namedtuple


AlleleSet = namedtuple("AlleleSet", "log_likelihood repeat_lengths frequencies")


def allele_maximum_likelihood(histogram: np.array, num_allele: int, probability_table: np.array) -> AlleleSet:
	#The search for the alleles will be only from lenghts that are larger than 5 reads.
	#This mean that I need to have two lists of availbile lenghts. All of them, and only those that are larger than 5 reads.
	assert num_allele >= 1, "NUMBER OF ALLELES MUST BE EQUAL TO OR GREATER THAN ONE"
	repeat_lengths = histogram[0,:]
	num_reads = histogram[1,:]

	supported_repeat_lengths = repeat_lengths[np.where(num_reads>5)[0]]
	max_log_likelihood=-1e9

	for _ in range(10):
		# randomly select num_alles repeat lengths
		random_repeat_lengths = supported_repeat_lengths[np.random.permutation(supported_repeat_lengths.size)[0:num_allele]]

		new_frequencies = np.zeros(num_allele)
		frequencies = np.ones(num_allele)/num_allele
		Z_i_j = np.zeros([44, num_allele])
		new_theta = np.zeros(num_allele)

		change = 1e6
		prev_log_likelihood = 1e6
		while change > 1e-5:

			#Step 1: make the Z_i_j matrix
			for length in repeat_lengths:
				for j in range(num_allele):
					Z_i_j[length, j] = probability_table[random_repeat_lengths[j], length]*frequencies[j]/np.sum(probability_table[random_repeat_lengths[:], length]*frequencies[:] + 1e-10)

			#Step 2: From the Z_i_j's estimate the new frequencies.
			for j in range(num_allele):
				new_frequencies[j] = np.sum(Z_i_j[repeat_lengths, j]*num_reads)/np.sum(num_reads)

			#Step number 2. Maximize the new Thetas
			Theta_new_temp = np.zeros(supported_repeat_lengths.size)
			for j in range(num_allele):
				for k in range(supported_repeat_lengths.size):
					Test_theta = supported_repeat_lengths[k]
					Theta_new_temp[k] = sum(Z_i_j[repeat_lengths, j] * np.log(probability_table[Test_theta, repeat_lengths] + 1e-10) * num_reads)
				new_theta[j] = supported_repeat_lengths[Theta_new_temp.argmax()]

			for j in range(num_allele):
				random_repeat_lengths[j] = new_theta[j]
				frequencies[j] = new_frequencies[j]

			#Calcualte the likelihood
			#L(Theta|D)=P(D|Theta)=PI_i(d_i|Theta)=PI_i(SUM_j(f[j]*p(d_i|theta_j))).
			#In our case we can combine all the reads with the same number of repeats thus:
			log_likelihood = 0
			for k in np.arange(repeat_lengths.size):
				log_likelihood += num_reads[k] * np.log(sum(frequencies * probability_table[random_repeat_lengths, repeat_lengths[k]]) + 1e-10 )
			if log_likelihood > max_log_likelihood:
				max_log_likelihood = log_likelihood
				best_alleles = random_repeat_lengths
				best_frequencies = frequencies
			change = np.abs(prev_log_likelihood - log_likelihood)
			prev_log_likelihood = log_likelihood
	return AlleleSet(log_likelihood = max_log_likelihood , repeat_lengths=best_alleles, frequencies=best_frequencies)


def find_alleles(histogram, probability_table) -> AlleleSet:
	repeat_lengths = histogram[0,:]
	num_reads = histogram[1,:]
	supported_repeat_lengths = repeat_lengths[np.where(num_reads>5)[0]]

	first_allele_set = allele_maximum_likelihood(histogram, 1, probability_table)
	second_allele_set = allele_maximum_likelihood(histogram, 2, probability_table)
	distribution_2_alleles = 2*(second_allele_set.log_likelihood - first_allele_set.log_likelihood)
	if distribution_2_alleles > 0:
		p_value_2_alleles = stats.chi2.pdf(distribution_2_alleles, 2)
		if p_value_2_alleles > 0.05:
			return first_allele_set
		elif supported_repeat_lengths.size == 2:
			return second_allele_set
		else:
			third_allele_set = allele_maximum_likelihood(histogram, 3, probability_table)
			distribution_3_alleles = 2*(third_allele_set.log_likelihood - second_allele_set.log_likelihood)
			p_value_3_alleles = stats.chi2.pdf(distribution_3_alleles, 2)
			if  p_value_3_alleles > 0.05:
				return second_allele_set
			elif supported_repeat_lengths.size == 3:
				return third_allele_set
			else:
				fourth_allele_set = allele_maximum_likelihood(histogram, 4, probability_table)
				distribution_4_alleles = 2*(fourth_allele_set.log_likelihood - third_allele_set.log_likelihood)
				p_value_4_alleles = stats.chi2.pdf(distribution_4_alleles,2)
				if p_value_4_alleles > 0.05:
					return third_allele_set
				else:
					return fourth_allele_set
	return first_allele_set


def get_repeat_threshold(ms_length):
	if ms_length == 1:
		return 5
	elif ms_length == 2:
		return 4
	elif ms_length >= 3:
		return 3


def filter_repeats(histogram, ms_length):
	""" Filters loci with very long or short repeats"""
	filtered = np.nonzero(histogram[0, :] < 40)  # insertion lengths of under 40
	histogram = histogram[:, filtered[0]]
	filtered = np.nonzero(histogram[0, :] > get_repeat_threshold(ms_length))
	return histogram[:, filtered[0]]


def write_results(filename, alleles: AlleleSet, split_histogram):
	with open(filename + ".all", "w") as results_file:
		for orig_substr in split_histogram: # writes original histogram
			results_file.write("{0} ".format(orig_substr))
		results_file.write("-999 ")
		results_file.write("{0} ".format(round(alleles.log_likelihood, 6)))
		results_file.write("-999 ")
		for repeat_length in alleles.repeat_lengths:
			results_file.write("{0} ".format(repeat_length))
		results_file.write("-999 ")
		for likelihood in alleles.frequencies:
			results_file.write("{0} ".format(round(likelihood, 6)))
		results_file.write("\n")
		results_file.close()


def main(cmd_arguments):
	probability_table = np.loadtxt(cmd_arguments[2], delimiter=',') # probability file
	histogram_file = open(cmd_arguments[1], "r") # histogram file
	for histogram_str in histogram_file:
		split_histogram=histogram_str.strip().split(", ")
		for i in range(1, len(split_histogram)):
			split_histogram[i] = int(split_histogram[i])
		histogram = np.array([split_histogram[1:len(split_histogram):2], split_histogram[2:len(split_histogram):2]])
		ms_length = float(split_histogram[0][split_histogram[0].rfind(":")+1:]) # length of microsatellite
		histogram = filter_repeats(histogram, ms_length)
		repeat_lengths = histogram[0,:]
		num_reads=histogram[1,:]
		supported_repeat_lengths = repeat_lengths[np.where(num_reads > 5)[0]]  # filters low coverage MSI loci

		if supported_repeat_lengths.size==1:
			allele_set = AlleleSet(log_likelihood=0, repeat_lengths=[supported_repeat_lengths], frequencies=[1])

		elif supported_repeat_lengths.size>1:
			allele_set = find_alleles(histogram, probability_table)

		write_results(cmd_arguments[1], allele_set, split_histogram)


if __name__ == '__main__':
	main(sys.argv)