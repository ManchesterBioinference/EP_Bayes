def prepare_upper_and_lower_bounds_for_priors(prior_elements, infered_elements):

	import config_variables
	np = config_variables.np
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided
	mode = config_variables.mode
	chroms_in_prior = config_variables.chroms_in_prior
	chroms_to_infer = config_variables.chroms_to_infer

	max_total_prior = {}
	min_total_prior = {}
	max_total_infer = {}
	min_total_infer = {}

	if one_sided_or_two_sided == "double_sided":

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			max_total_prior[classification_of_interactions] = {}
			min_total_prior[classification_of_interactions] = {}
			max_total_infer[classification_of_interactions] = {}
			min_total_infer[classification_of_interactions] = {}

			for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):

				distance_values_prior = prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"]
				distance_values_infer = infered_elements[mode][classification_of_interactions]["distance"]["attribute_values"]

				max_total_prior[classification_of_interactions][positive_or_negative_side] = max(np.array([max(sign*distance_values_prior[chrom_][sign*distance_values_prior[chrom_] > 0]) for chrom_ in chroms_in_prior]))
				min_total_prior[classification_of_interactions][positive_or_negative_side] = min(np.array([min(sign*distance_values_prior[chrom_][sign*distance_values_prior[chrom_] > 0]) for chrom_ in chroms_in_prior]))

				max_total_infer[classification_of_interactions][positive_or_negative_side] = max(np.array([max(sign*distance_values_infer[chrom_][sign*distance_values_infer[chrom_] > 0]) for chrom_ in chroms_to_infer]))
				min_total_infer[classification_of_interactions][positive_or_negative_side] = min(np.array([min(sign*distance_values_infer[chrom_][sign*distance_values_infer[chrom_] > 0]) for chrom_ in chroms_to_infer]))

		up_dist = {}
		low_dist = {}
		for sign, positive_or_negative_side in zip([1, -1],["positive_side", "negative_side"]):
			up_dist[positive_or_negative_side] = []
			low_dist[positive_or_negative_side] = []
			for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
				up_dist[positive_or_negative_side] = up_dist[positive_or_negative_side] + [max([max_total_infer[classification_of_interactions][positive_or_negative_side], max_total_prior[classification_of_interactions][positive_or_negative_side]])]
				low_dist[positive_or_negative_side] = low_dist[positive_or_negative_side] + [min([min_total_infer[classification_of_interactions][positive_or_negative_side], min_total_prior[classification_of_interactions][positive_or_negative_side]])]
			
			up_dist[positive_or_negative_side] = (1. + 10.**-15)* max(up_dist[positive_or_negative_side])
			low_dist[positive_or_negative_side] = (1. - 10.**-15)* min(low_dist[positive_or_negative_side])


		swap_limits_upper = up_dist['negative_side']
		swap_limits_lower = low_dist['negative_side']
		up_dist['negative_side'] = -1*swap_limits_lower
		low_dist['negative_side'] = -1*swap_limits_upper

		

		return low_dist, up_dist

	elif one_sided_or_two_sided == "single_sided":

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			max_total_prior[classification_of_interactions] = {}
			min_total_prior[classification_of_interactions] = {}
			max_total_infer[classification_of_interactions] = {}
			min_total_infer[classification_of_interactions] = {}

			distance_values_prior = prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"]
			distance_values_infer = infered_elements[mode][classification_of_interactions]["distance"]["attribute_values"]

			max_total_prior[classification_of_interactions] = max([max(distance_values_prior[chrom_]) for chrom_ in chroms_in_prior])
			min_total_prior[classification_of_interactions] = min([min(distance_values_prior[chrom_]) for chrom_ in chroms_in_prior])

			max_total_infer[classification_of_interactions] = max([max(distance_values_infer[chrom_]) for chrom_ in chroms_to_infer])
			min_total_infer[classification_of_interactions] = min([min(distance_values_infer[chrom_]) for chrom_ in chroms_to_infer])

		up_dist = []
		low_dist = []

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			up_dist = up_dist + [max([max_total_infer[classification_of_interactions], max_total_prior[classification_of_interactions]])]
			low_dist = low_dist + [min([min_total_infer[classification_of_interactions], min_total_prior[classification_of_interactions]])]
			
		up_dist = (1. + 10.**-15)*max(up_dist)
		low_dist = (1. - 10.**-15)*min(low_dist)

		return low_dist, up_dist



