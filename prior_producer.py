def prior_producer():

	import config_variables
	import correl_distance_extractor_clean
	np = config_variables.np
	log_distances = config_variables.log_distances
	chroms_in_prior	= config_variables.chroms_in_prior
	dataset_names_option = config_variables.dataset_names_option
	mode = config_variables.mode
	domain = config_variables.domain
	negative_interactions = config_variables.negative_interactions
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided
	domain_like_chromosome_correction = config_variables.domain_like_chromosome_correction

	prior_elements = {}


	if config_variables.MoG_classificator: 
		prior_mode = False
	else: prior_mode = True


	for type_of_interaction in ["promoter_enhancer_interactions", "enhancer_enhancer_interactions"]:

		prior_elements[type_of_interaction] = {}

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			prior_elements[type_of_interaction][classification_of_interactions] = {}

			for attribute_of_interaction in ["distance", "correlation"]:
				prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction] = {}

				if attribute_of_interaction == "correlation":

					for data_set_name in dataset_names_option:
						prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name] = {}
						prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name]["attribute_values"] = {}
						prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name]["prior_bins"] = np.array([])
						prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name]["prior_frequencies"] = np.array([])

						for chrom_ in chroms_in_prior:
							prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name]["attribute_values"][chrom_] = np.array([])

				else:
					prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction]["attribute_values"] = {}
					prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction]["prior_bins"] = np.array([])
					prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction]["prior_frequencies"] = np.array([])

					for chrom_ in chroms_in_prior:
						prior_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction]["attribute_values"][chrom_] = np.array([])

	if not(domain) and not(domain_like_chromosome_correction):
		if one_sided_or_two_sided == "double_sided":
			if log_distances:
				prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = False, logged = True, prior_mode = prior_mode)
			else:
				prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = False, logged = False, prior_mode = prior_mode)

		elif one_sided_or_two_sided == "single_sided":

			if log_distances:
				prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = True, logged = True, prior_mode = prior_mode)
			else: 
				prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = True, logged = False, prior_mode = prior_mode)


	else:

		if one_sided_or_two_sided == "double_sided":
			prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = False, logged = False, prior_mode = prior_mode)
		elif one_sided_or_two_sided == "single_sided":
			prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = True, logged = False, prior_mode = prior_mode)

		import interacting_domain_clean as interacting_domain 
		from prepare_interactions_clean import un_string
		dataset_time_series_dict = config_variables.dataset_time_series_dict
		link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
		chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
		chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh
		
	
		def interacting_enhancers_coordinates_function():

			interacting_enhancers_coord = {}

			enh_coordinates = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"]['ER']][1]

			for chrom in  chroms_in_prior:		
				if mode == "promoter_enhancer_interactions":
					chrom_interacting_enhancers_pro = np.unique(un_string(chr_interactions_dict_pro_enh[chrom])[:,1])
					interacting_enhancers_coord[chrom] = enh_coordinates[chrom_interacting_enhancers_pro]
				#if config_variables.alternative_classificator_outside_enhancers: 
					#interacting_enhancers_coord[chrom] = enh_coordinates[config_variables.chrom_interacting_enhancers_pro[chrom]]
					
				elif mode == "enhancer_enhancer_interactions": 

					chrom_interacting_enhancers_enh = np.unique(un_string(chr_interactions_dict_enh_enh[chrom])[:,0])
					chrom_interacting_enhancers_enh = np.unique(np.r_[chrom_interacting_enhancers_enh, np.unique(un_string(chr_interactions_dict_enh_enh[chrom])[:,1])])
					interacting_enhancers_coord[chrom] = enh_coordinates[chrom_interacting_enhancers_enh]

			return interacting_enhancers_coord

		prior_elements[mode]["interacting_enhancers_coordinates"] = interacting_enhancers_coordinates_function()

		#chrom, interacting_enhancer_coordinates, distances, possible_distances_counts, length = chrom, interacting_enhancer_coordinates_chrom, total_array, possible_distances_counts, 200
		def adaptive_domain_dummy_inter_matrix_version(chrom, interacting_enhancer_coordinates, distances, possible_distances_counts, length):


			possible_interaction_centres = interacting_enhancer_coordinates.mean(1)[:, None] + distances #should be if dist - then one end if dist + then other end
			reshape_size = possible_interaction_centres.size
			possible_interaction_coordinates = np.c_[(possible_interaction_centres.reshape(reshape_size) - length)[:, None], (possible_interaction_centres.reshape(reshape_size) + length)[:, None]]


			def state_dependent_allocation(state, chrom):
				allocation_of_coordinates, size_domains = interacting_domain.interacting_domains(interacting_enhancer_coordinates, possible_interaction_coordinates, chrom, state, matrix_version = False, chromosome_version = domain_like_chromosome_correction) # matrix_mode = False

				allocations_of_interacting_enhancers = allocation_of_coordinates[:len(interacting_enhancer_coordinates)]

				allocations_of_interacting_enhancer_possible_interactions = allocation_of_coordinates[len(interacting_enhancer_coordinates):]

				allocations_of_interacting_enhancer_possible_interactions = allocations_of_interacting_enhancer_possible_interactions.reshape(len(interacting_enhancer_coordinates), len(distances))	

				dummy_interactions_in_the_same_domain_as_enhancer = allocations_of_interacting_enhancers[:,None] == allocations_of_interacting_enhancer_possible_interactions
	
				return dummy_interactions_in_the_same_domain_as_enhancer, size_domains

			dummy_interactions_in_the_same_domain_as_enhancer_left, size_domains = state_dependent_allocation("left", chrom)

			dummy_interactions_in_the_same_domain_as_enhancer_right, size_domains = state_dependent_allocation("right", chrom)

			possible_distances_counts_chrom = ((dummy_interactions_in_the_same_domain_as_enhancer_left + dummy_interactions_in_the_same_domain_as_enhancer_right).astype(int)).sum(0)

			possible_distances_counts = possible_distances_counts + possible_distances_counts_chrom

			return possible_distances_counts, size_domains

	
		chrom_domain_sizes = {}	

		import itertools

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			total_array = [prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
			total_array = list(itertools.chain.from_iterable(total_array))
			possible_distances_counts = np.zeros_like(total_array)	
			for chrom in  chroms_in_prior:
				interacting_enhancer_coordinates_chrom = prior_elements[mode]["interacting_enhancers_coordinates"][chrom]
				possible_distances_counts, size_domains_chrom = adaptive_domain_dummy_inter_matrix_version(chrom, interacting_enhancer_coordinates_chrom, total_array, possible_distances_counts, 200)
		
				chrom_domain_sizes[chrom] = size_domains_chrom

			prior_elements[mode][classification_of_interactions]["distance"]["possible_distances_counts"] = possible_distances_counts	

			if any(possible_distances_counts == 0):
				#interacting_enhancer_coordinates = prior_elements[mode]["interacting_enhancers_coordinates"]
				#possible_distances_counts, total_array = corrector(interacting_enhancer_coordinates, total_array, possible_distances_counts, 200)

				if mode == "promoter_enhancer_interactions" and classification_of_interactions == "positive_interactions":
					print "possible_distances_p_e_counts_true should not have 0 count distances which could be corrected by corrector() i.e something wrong"

		if one_sided_or_two_sided == "double_sided":
			if log_distances:
				prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = False, logged = True, prior_mode = prior_mode)

		elif one_sided_or_two_sided == "single_sided":
			if log_distances:
				prior_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, prior_elements, chroms_in_prior, absolute = True, logged = True, prior_mode = prior_mode)


	return prior_elements

