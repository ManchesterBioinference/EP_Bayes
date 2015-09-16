def correl_distance_extractor(mode, dictionary_elements, chroms, absolute = False, logged = False, prior_mode = False):
	import config_variables
	negative_interactions = config_variables.negative_interactions
	dataset_names_option = config_variables.dataset_names_option
	np = config_variables.np


	low_dist = up_dist = []
	for chrom in chroms:

		type_of_interaction_false_interactions = negative_interactions.chrom_specific_negative_interactions(chrom, mode, prior_mode)
		distance_of_positive_type_of_interaction, distance_of_negative_type_of_interaction = negative_interactions.distance_and_correl_of_interactions(chrom, type_of_interaction_false_interactions, mode, attribute_distance = True)

		if absolute:

			distance_of_positive_type_of_interaction = abs(distance_of_positive_type_of_interaction)
			distance_of_negative_type_of_interaction = abs(distance_of_negative_type_of_interaction)


		if logged:

			mask_positive_distances_pos_inter = np.zeros(len(distance_of_positive_type_of_interaction), bool)
			mask_positive_distances_pos_inter[distance_of_positive_type_of_interaction > 0] = True
			mask_negative_distances_pos_inter = np.invert(mask_positive_distances_pos_inter)

			distance_of_positive_type_of_interaction[mask_positive_distances_pos_inter] = np.log10(distance_of_positive_type_of_interaction[mask_positive_distances_pos_inter])
			distance_of_positive_type_of_interaction[mask_negative_distances_pos_inter] = -1*(np.log10(-1*distance_of_positive_type_of_interaction[mask_negative_distances_pos_inter]))


			mask_positive_distances_neg_inter = np.zeros(len(distance_of_negative_type_of_interaction), bool)
			mask_positive_distances_neg_inter[distance_of_negative_type_of_interaction > 0] = True
			mask_negative_distances_neg_inter = np.invert(mask_positive_distances_neg_inter)

			distance_of_negative_type_of_interaction[mask_positive_distances_neg_inter] = np.log10(distance_of_negative_type_of_interaction[mask_positive_distances_neg_inter])
			distance_of_negative_type_of_interaction[mask_negative_distances_neg_inter] = -1*(np.log10(-1*distance_of_negative_type_of_interaction[mask_negative_distances_neg_inter]))


		dictionary_elements[mode]["positive_interactions"]["distance"]["attribute_values"][chrom] = distance_of_positive_type_of_interaction
		dictionary_elements[mode]["negative_interactions"]["distance"]["attribute_values"][chrom] = distance_of_negative_type_of_interaction

		
		for data_set_name in dataset_names_option:
	
			correlation_of_positive_type_of_interaction, correlation_of_negative_type_of_interaction = negative_interactions.distance_and_correl_of_interactions(chrom, type_of_interaction_false_interactions, mode, data_set_name, attribute_correlation = True)

			dictionary_elements[mode]["positive_interactions"]["correlation"][data_set_name]["attribute_values"][chrom] = correlation_of_positive_type_of_interaction
			dictionary_elements[mode]["negative_interactions"]["correlation"][data_set_name]["attribute_values"][chrom] = correlation_of_negative_type_of_interaction

		mean_correl_pos, mean_correl_neg = np.zeros_like(correlation_of_positive_type_of_interaction), np.zeros_like(correlation_of_negative_type_of_interaction)
		for data_set_name in dataset_names_option[:2]:
			mean_correl_pos = mean_correl_pos + dictionary_elements[mode]["positive_interactions"]["correlation"][data_set_name]["attribute_values"][chrom]
			mean_correl_neg = mean_correl_neg + dictionary_elements[mode]["negative_interactions"]["correlation"][data_set_name]["attribute_values"][chrom]

		for data_set_name in dataset_names_option[:2]:
			dictionary_elements[mode]["positive_interactions"]["correlation"][data_set_name]["attribute_values"][chrom] = mean_correl_pos/2.
			dictionary_elements[mode]["negative_interactions"]["correlation"][data_set_name]["attribute_values"][chrom] = mean_correl_neg/2.


	return dictionary_elements

