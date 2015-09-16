def infered_elements_filler():

	import config_variables
	import correl_distance_extractor_clean
	mode = config_variables.mode
	domain = config_variables.domain
	chroms_to_infer = config_variables.chroms_to_infer
	np = config_variables.np
	dataset_names_option = config_variables.dataset_names_option
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided
	log_distances = config_variables.log_distances

	infered_elements = {} 

	for type_of_interaction in ["promoter_enhancer_interactions", "enhancer_enhancer_interactions"]:

		infered_elements[type_of_interaction] = {}

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			infered_elements[type_of_interaction][classification_of_interactions] = {}

			for attribute_of_interaction in ["distance", "correlation"]:
				infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction] = {}

				if attribute_of_interaction == "correlation":

					for data_set_name in dataset_names_option:
						infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name] = {}
						infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name]["attribute_values"] = {}
						for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
							infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name][probability_of_being_positive_or_negative] = {}

						for chrom_ in chroms_to_infer:
							infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name]["attribute_values"][chrom_] = np.array([])
							for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
								infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][data_set_name][probability_of_being_positive_or_negative][chrom_] = np.array([])
	
				else:
					infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction]["attribute_values"] = {}
					for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
						infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative] = {}

					for chrom_ in chroms_to_infer:
						infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction]["attribute_values"][chrom_] = np.array([])
						for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
							infered_elements[type_of_interaction][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][chrom_] = np.array([])

	if one_sided_or_two_sided == "double_sided":
		if log_distances: 
			infered_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, infered_elements, chroms_to_infer, absolute = False, logged = True)
		else:
			infered_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, infered_elements, chroms_to_infer, absolute = False, logged = False)

	elif one_sided_or_two_sided == "single_sided":

		if log_distances: 
			infered_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, infered_elements, chroms_to_infer, absolute = True, logged = True)
		else: 
			infered_elements = correl_distance_extractor_clean.correl_distance_extractor(mode, infered_elements, chroms_to_infer, absolute = True, logged = False)

	return infered_elements

