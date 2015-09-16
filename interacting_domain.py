def interacting_domains(chrom_pro_coord, chrom_ER_coord, chrom, state, matrix_version = True): # state 'left' or 'right'		
	import numpy as np
	import config_variables
	data_folder = config_variables.data_folder

	interacting_domains = np.loadtxt(data_folder+'report_hESC_Combined_converted.csv.gz', dtype = str, usecols = (4, 12, 13), delimiter = ',')
	interacting_domains_chrom = interacting_domains[interacting_domains[:, 0] == chrom][:, 1:].astype(int)

	size_domains = len(interacting_domains_chrom)

	interacting_domains_chrom[:, 0] = interacting_domains_chrom[:, 0] - 1

	bins = np.unique(interacting_domains_chrom)
	bins = np.r_[0, bins, 2*bins[-1]]

	features = np.r_[chrom_pro_coord, chrom_ER_coord]

	first_col_allocat = np.digitize(features[:, 0], bins) - 1
	second_col_allocat = np.digitize(features[:, 1], bins) - 1
	mean_allocat = np.digitize(features.mean(1), bins) - 1

	selective = False # cares about whether a domain exists or not
	
	existing_domains = np.zeros(bins.size, bool)

	for index in np.arange(bins.size-1):
		if [bins[index], bins[index+1]] in map(list, interacting_domains_chrom):
			existing_domains[index] = True
	
	fall_into_domains_1st_col = np.in1d(first_col_allocat, np.where(existing_domains))
	fall_into_domains_2nd_col = np.in1d(second_col_allocat, np.where(existing_domains))
	fall_into_domains_mean = np.in1d(mean_allocat, np.where(existing_domains))

	falling_into_a_domain_mask = fall_into_domains_1st_col + fall_into_domains_2nd_col + fall_into_domains_mean

	consensus_allocation = -1*np.ones(len(features))

	for index, falls_into_domain, el in zip(np.arange(len(features)), falling_into_a_domain_mask, np.c_[first_col_allocat[:, None], second_col_allocat[:, None], mean_allocat[:, None]]):
		if selective:
			if falls_into_domain:
				unique = np.unique(el)
				if len(unique) == 1: consensus_allocation[index] = unique
				elif len(unique) > 1:
					if state == 'left': consensus_allocation[index] = min(unique)
					elif state == 'right': consensus_allocation[index] = max(unique)
		else:
			unique = np.unique(el)
			if len(unique) == 1: consensus_allocation[index] = unique
			elif len(unique) > 1:
				if state == 'left': consensus_allocation[index] = min(unique)
				elif state == 'right': consensus_allocation[index] = max(unique)

	length_chr = len(chrom_pro_coord) + len(chrom_ER_coord)
	if matrix_version == True:
		state_matrix = np.zeros((length_chr, length_chr), bool)

		for index, el in enumerate(consensus_allocation): #consensus_allocation stores a bin allocations of each feature 
			if el <> -1: 								  #if selective is off then it will be always true, eitherwards it may be false and then the state matrix will be false
				state_matrix[index] = consensus_allocation == el

		return state_matrix

	else: 
		return consensus_allocation, size_domains

