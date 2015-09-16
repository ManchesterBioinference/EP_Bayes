def allocator(infered_elements, prior_elements):

	import config_variables
	import correl_distance_extractor_clean
	mode = config_variables.mode
	domain = config_variables.domain
	chroms_to_infer = config_variables.chroms_to_infer
	np = config_variables.np
	dataset_names_option = config_variables.dataset_names_option
	negative_interactions = config_variables.negative_interactions
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided

	def allocate_to_empirical_probabilities_in_histogram(array_to_allocate, adequate_histogram_bins, probabilities_of_a_bin):

		init_shape = array_to_allocate.shape
		
		array_to_allocate = array_to_allocate.reshape(array_to_allocate.size)
		
		indexes_of_bins_it_falls_into = np.digitize(array_to_allocate, adequate_histogram_bins)#, right=False)
		indexes_of_nan = np.isnan(array_to_allocate)
		indexes_of_bins_it_falls_into[indexes_of_nan] = 1

		if np.min(indexes_of_bins_it_falls_into) == 0 or np.max(indexes_of_bins_it_falls_into) == len(adequate_histogram_bins):
			print 'allocation_problem', np.min(indexes_of_bins_it_falls_into) == 0, np.max(indexes_of_bins_it_falls_into) == len(adequate_histogram_bins)
			print np.where(indexes_of_bins_it_falls_into == 0)[0], array_to_allocate[np.where(indexes_of_bins_it_falls_into == 0)[0]]
					
			output_probabilities = probabilities_of_a_bin[indexes_of_bins_it_falls_into - 1]
			output_probabilities[indexes_of_nan] = np.min(probabilities_of_a_bin)*0.9
			
			return indexes_of_bins_it_falls_into.reshape(init_shape), output_probabilities.reshape(init_shape)#, indexes_of_bins_it_falls_into.reshape(init_shape)

		else:
		
			output_probabilities = probabilities_of_a_bin[indexes_of_bins_it_falls_into - 1]
			output_probabilities[indexes_of_nan] = np.min(probabilities_of_a_bin)*0.9
			return indexes_of_bins_it_falls_into.reshape(init_shape), output_probabilities.reshape(init_shape)#, indexes_of_bins_it_falls_into.reshape(init_shape)



		print 'enhancer-enhancers allocation'
	

#	def bin_extender(bins, prob):
#		mass_left, mass_right = np.diff(bins[:2])*prob[0], np.diff(bins[-2:])*prob[-1]
#		bins[0], bins[-1] = low_dist, up_dist	
#		prob[0], prob[-1] = mass_left/np.diff(bins[:2]), mass_right/np.diff(bins[-2:])

#		return bins, prob

	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		for attribute_of_interaction in ["distance", "correlation"]:
			if attribute_of_interaction == "correlation":
				for data_set_name in dataset_names_option:
					for chrom_ in chroms_to_infer:				
						array_to_allocate = infered_elements[mode][classification_of_interactions][attribute_of_interaction][data_set_name]["attribute_values"][chrom_]
						for classification_of_interactions_2 in ["positive_interactions", "negative_interactions"]:

							probabilities_of_a_bin = prior_elements[mode][classification_of_interactions_2][attribute_of_interaction][data_set_name]["prior_frequencies"]
							adequate_histogram_bins = prior_elements[mode][classification_of_interactions_2][attribute_of_interaction][data_set_name]["prior_bins"]
	
							bin_allocation_probabilities = allocate_to_empirical_probabilities_in_histogram(array_to_allocate, adequate_histogram_bins, probabilities_of_a_bin)[1]
							infered_elements[mode][classification_of_interactions][attribute_of_interaction][data_set_name]["probabilities_of_being_{0}".format(classification_of_interactions_2)][chrom_] = bin_allocation_probabilities


			else:
				for chrom_ in chroms_to_infer:
					array_to_allocate = infered_elements[mode][classification_of_interactions][attribute_of_interaction]["attribute_values"][chrom_]
					for classification_of_interactions_2 in ["positive_interactions", "negative_interactions"]:	

						probabilities_of_a_bin = prior_elements[mode][classification_of_interactions_2][attribute_of_interaction]["prior_frequencies"]
						adequate_histogram_bins = prior_elements[mode][classification_of_interactions_2][attribute_of_interaction]["prior_bins"]

						#if  low_dist < adequate_histogram_bins[0] or adequate_histogram_bins[-1] < up_dist: adequate_histogram_bins, probabilities_of_a_bin = bin_extender(adequate_histogram_bins, probabilities_of_a_bin)
						
						bin_allocation_probabilities = allocate_to_empirical_probabilities_in_histogram(array_to_allocate, adequate_histogram_bins, probabilities_of_a_bin)[1]
						infered_elements[mode][classification_of_interactions][attribute_of_interaction]["probabilities_of_being_{0}".format(classification_of_interactions_2)][chrom_] = bin_allocation_probabilities


	return infered_elements
