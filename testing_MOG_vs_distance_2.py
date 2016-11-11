	def calculate_distances_between_promoters_and_non_promoters():
		"""
		as the name suggests 
		"""
		from chrom_specific_negative_interactions import extract_TSS_coordinates, calculate_distances

		TSS_coordinates = extract_TSS_coordinates(config_variables.upstream)

		indexes_p = total_p + np.arange(len(pro_coords[pro_chroms == chrom]))#dict_chrom_pro_survived[chrom]
		indexes_e = chrom_enh_survived

		point_coordinates_promoter, point_coordinates_enhancer = TSS_coordinates[indexes_p], np.mean(enh_coords[indexes_e], axis = 1)

		distances_matrix = calculate_distances(False, point_coordinates_promoter, point_coordinates_enhancer)

		dist_from_enhancer_to_promoters = distances_matrix[len(indexes_p):, 0: len(indexes_p)]

		return dist_from_enhancer_to_promoters

		#dist=np.array([abs(np.mean(peak) - np.mean(active_promoters, axis = 1)) for peak in peaks])


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

			dist = calculate_distances_between_promoters_and_non_promoters()

			neg_dist = dist < 0

			if log_distances:
				negative_distances = -1*dist[neg_dist]
				positive_distances = dist[np.invert(neg_dist)]
				negative_distances = -1*np.log10(negative_distances)
				positive_distances = np.log10(positive_distances)
			else: 
				negative_distances = dist[neg_dist]
				positive_distances = dist[np.invert(neg_dist)]

			dist[neg_dist] = negative_distances
			dist[np.invert(neg_dist)] = positive_distances 

			probabilities_of_a_bin = prior_elements[mode]["negative_interactions"]["distance"]["prior_frequencies"]
			adequate_histogram_bins = prior_elements[mode]["negative_interactions"]["distance"]["prior_bins"]

			distance_prob_component_neg_neg = allocate_to_empirical_probabilities_in_histogram(negative_distances, adequate_histogram_bins, probabilities_of_a_bin)[1]
			distance_prob_component_pos_neg = allocate_to_empirical_probabilities_in_histogram(positive_distances, adequate_histogram_bins, probabilities_of_a_bin)[1]

			distance_prob_component_neg = np.zeros_like(dist)
			distance_prob_component_neg[neg_dist] = distance_prob_component_neg_neg
			distance_prob_component_neg[np.invert(neg_dist)] = distance_prob_component_pos_neg

			dist = calculate_distances_between_promoters_and_non_promoters()

			neg_dist = dist < 0

			if log_distances:
				negative_distances = -1*dist[neg_dist]
				positive_distances = dist[np.invert(neg_dist)]
				negative_distances = -1*np.log10(negative_distances)
				positive_distances = np.log10(positive_distances)
			else: 
				negative_distances = dist[neg_dist]
				positive_distances = dist[np.invert(neg_dist)]

			dist[neg_dist] = negative_distances
			dist[np.invert(neg_dist)] = positive_distances 

			probabilities_of_a_bin = prior_elements[mode]["positive_interactions"]["distance"]["prior_frequencies"]
			adequate_histogram_bins = prior_elements[mode]["positive_interactions"]["distance"]["prior_bins"]

			distance_prob_component_neg_pos = allocate_to_empirical_probabilities_in_histogram(negative_distances, adequate_histogram_bins, probabilities_of_a_bin)[1]
			distance_prob_component_pos_pos = allocate_to_empirical_probabilities_in_histogram(positive_distances, adequate_histogram_bins, probabilities_of_a_bin)[1]

			distance_prob_component_pos = np.zeros_like(dist)
			distance_prob_component_pos[neg_dist] = distance_prob_component_neg_pos
			distance_prob_component_pos[np.invert(neg_dist)] = distance_prob_component_pos_pos













