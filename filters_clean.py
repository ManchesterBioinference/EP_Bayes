import numpy as np
import config_variables
full_list_promoters = config_variables.full_list_promoters
option = config_variables.option
dataset_time_series_dict = config_variables.dataset_time_series_dict
time_points = config_variables.time_points
datasets_names = config_variables.datasets_names
distant_enh_only = config_variables.distant_enh_only
full_list_enhancers = config_variables.full_list_enhancers
ER_pro_filtered_ = config_variables.ER_pro_filtered_
chrom_names = config_variables.chrom_names
link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name


def extract_TSS_coordinates(upstream):

	data = np.loadtxt(name_of_time_series_promoter_file_for_TSS_start, dtype = str,  delimiter = '\t')	
	plus_strand = data[:, 4] == '+'
	TSS_coordinates = np.zeros(len(plus_strand), int)
	TSS_coordinates[plus_strand] = data[plus_strand, 1].astype(int) + upstream
	TSS_coordinates[np.invert(plus_strand)] = data[np.invert(plus_strand), 2].astype(int) + upstream

	return TSS_coordinates


def single_element_cleaner(list_of_datasets, remove_single_domain_elements, filtered_elements):
	import itertools
	import interacting_domain

	mode = config_variables.mode
	upstream = config_variables.upstream
	TSS_or_intra_genic_for_domain_filter = config_variables.TSS_or_intra_genic_for_domain_filter
	chrom_mask_non_single_domain_elements = {}

	chroms, coordinates = dataset_time_series_dict[list_of_datasets[0]][0], dataset_time_series_dict[list_of_datasets[0]][1]

	interacting_domains = np.loadtxt('report_hESC_Combined_converted.csv', dtype = str, usecols = (4, 12, 13), delimiter = ',')

	if mode == "promoter_enhancer_interactions":
				
		if TSS_or_intra_genic_for_domain_filter == "TSS_only": 
			TSS_coordinates = extract_TSS_coordinates(upstream)
			coordinates = np.column_stack((TSS_coordinates-1, TSS_coordinates+1))

	for chrom_ in np.unique(chroms):
		filtered_elements_chrom = filtered_elements[chroms == chrom_]
		chrom_coordinates = coordinates[chroms == chrom_]

		if len(chrom_coordinates) and sum(interacting_domains[:, 0] == chrom_):
			matrix_left = interacting_domain.interacting_domains(chrom_coordinates, np.array([]).reshape(0,2), chrom_, state = "left", matrix_version = True)
			matrix_right = interacting_domain.interacting_domains(chrom_coordinates, np.array([]).reshape(0,2), chrom_, state = "right", matrix_version = True)
			mask = np.ones_like(matrix_left)
			mask[range(len(mask)), range(len(mask))] = False
			mask[:, np.invert(filtered_elements_chrom)] = False
			matrix_allocations_joint = matrix_left*mask + matrix_right*mask

			#number_of_elemenets_in_the_same_domain = matrix_allocations.sum(1)

			matrix_allocations_unique_to_left = matrix_allocations_joint - matrix_right*mask
			matrix_allocations_unique_to_right = matrix_allocations_joint - matrix_left*mask

			shared_allocations = matrix_allocations_joint - (matrix_allocations_unique_to_left + matrix_allocations_unique_to_right)

			promoter_is_shared = shared_allocations.sum(1).astype(bool) # if any is true

			survived = np.zeros(len(chrom_coordinates), bool)
			survived[promoter_is_shared] = ((matrix_allocations_unique_to_left.sum(1) > 0) + (matrix_allocations_unique_to_right.sum(1) > 0))[promoter_is_shared]

			survived[np.invert(promoter_is_shared)] = ((matrix_left*mask).sum(1) > 0)[np.invert(promoter_is_shared)]


			chrom_mask_non_single_domain_elements[chrom_] = survived

		else: chrom_mask_non_single_domain_elements[chrom_] = np.ones(sum(chroms == chrom_), bool)
	
	chrom_mask_non_single_domain_elements_total = np.array(list(itertools.chain.from_iterable([chrom_mask_non_single_domain_elements[chrom__] for chrom__ in np.unique(chroms)])))

	return chrom_mask_non_single_domain_elements_total


def features_filtered(filter_, count_f, list_of_datasets, options, name_of_overlap_file, add_overl, remove_single_domain_elements = False):

	#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


	rep_1_ts, rep_2_ts = dataset_time_series_dict[list_of_datasets[0]][2][:,:time_points], dataset_time_series_dict[list_of_datasets[1]][2][:,:time_points]

	filter_ = float(filter_)

	correlations = np.zeros(len(rep_1_ts))

	for ind, el in enumerate(np.c_[rep_1_ts, rep_2_ts]): correlations[ind] = np.corrcoef(el[:time_points], el[time_points:])[0,1]

	correl_filt = np.invert(np.isnan(correlations))
	correlations[np.isnan(correlations)] = -10.

	correl_filt = (correlations > filter_)*correl_filt
	
	filt = correl_filt

	for name_of_dataset_to_open, norm_name in zip(list_of_datasets[option], datasets_names[option]):
		              
		dataset_time_series = dataset_time_series_dict[name_of_dataset_to_open][2][:,:time_points]
		
		counts_k = (dataset_time_series).sum(1)#*norm_k[0]
		
		counts_filt_k = counts_k > count_f

		filt = filt*counts_filt_k

	#creates a mask for overlaps---------------------------		
	if add_overl:
		overl_features_ = np.unique(np.loadtxt(name_of_overlap_file, dtype = int, delimiter = '\t', usecols = (3, )))
		overl = np.zeros(len(rep_1_ts), bool)
		overl[overl_features_] = True	
		filt = filt + overl
	#------------------------------------------------------
	#if remove_single_domain_elements:
	#	chrom_mask_non_single_domain_elements_total = single_element_cleaner(list_of_datasets, remove_single_domain_elements, filt)
	#else: chrom_mask_non_single_domain_elements_total = True

	#filt = filt * chrom_mask_non_single_domain_elements_total
	
	dict_chrom_not_survived = {}
	dict_chrom_survived = {}

	for i in np.concatenate((np.array(range(1, 23), dtype='S2'), ['X'], ['Y'])):
		chrom_mask = dataset_time_series_dict[list_of_datasets[0]][0] == 'chr{0}'.format(i)

		dict_chrom_survived['chr{0}'.format(i)] = np.where(filt*chrom_mask)[0]
		dict_chrom_not_survived['chr{0}'.format(i)] = np.where(np.invert(filt)*chrom_mask)[0]


	return dict_chrom_survived, dict_chrom_not_survived, filt, correl_filt

def distant_enh_only_filter(name_of_overlap_file):
	chroms = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"]["ER"]][0]
	overlapping_peaks = np.loadtxt(name_of_overlap_file, usecols = (3,), dtype = int)
	filt = np.zeros(len(chroms), bool)
	filt[overlapping_peaks] = True
	dict_chrom_distant = {}
	dict_chrom_proximal = {}
	for chrom in chrom_names:
		chrom_mask = chroms == chrom
		dict_chrom_proximal[chrom] = np.where(filt*chrom_mask)[0]
		dict_chrom_distant[chrom] = np.where(np.invert(filt)*chrom_mask)[0]
	return dict_chrom_distant, dict_chrom_proximal, filt


