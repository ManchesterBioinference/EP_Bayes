def executor(posterior, selected_combinations, type_of_models):

	import config_variables 
	from prepare_interactions_clean import un_string
	import MAP_function_clean as MAP_
	import classifiers_clean
	import itertools
	mode = config_variables.mode
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided
	domain = config_variables.domain
	np = config_variables.np
	chroms_to_infer = config_variables.chroms_to_infer
	chroms_in_prior = config_variables.chroms_in_prior
	results_folder = config_variables.results_folder
	MoG_classificator = config_variables.MoG_classificator
	dict_option_ = config_variables.dict_option
	
	#f1 = open(results_folder + "MAP_performance_{0}_{1}_p_{2}_i_{3}_d{4}".format(mode, one_sided_or_two_sided, "-".join(chroms_in_prior[[0,-1]]), ",".join(chroms_to_infer), domain), "w")


	def get_MAPS_for_domain_non_domain_interacting_enhancers(domain_atr):

		chrom_interacting_enhancers_pro = config_variables.chrom_interacting_enhancers_pro

		from  prepare_interactions_clean import un_string
		interacting_non_intracting_mask = {}

		#total_number_of_interacting_enhancers = 0

		if config_variables.disentagled_features_validation: 
			chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS
		else:
			chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh

		for chrom___ in chroms_to_infer:

			mask = np.in1d(np.unique(un_string(chr_interactions_dict_pro_enh[chrom___])[:,1]), chrom_interacting_enhancers_pro[chrom___])
			if domain_atr == "outside_domain":
				interacting_non_intracting_mask[chrom___] = np.invert(mask)
			elif domain_atr == "within_domain":	
				interacting_non_intracting_mask[chrom___] = mask

			#total_number_of_interacting_enhancers += np.sum(interacting_non_intracting_mask[chrom___])

		return interacting_non_intracting_mask


	sensitivity_match_MAP = {}
	match_MAP = {}
	infered_elements = {}
	MAP_probabilites = {}
	#total_number_of_interacting_enhancers = {}

	probabilities_for_promoters_of_interacting_enhancers = {}	
	for type_of_model in type_of_models:
		sensitivity_match_MAP[type_of_model] = {}
		match_MAP[type_of_model] = {}
		infered_elements[type_of_model] = {}
		MAP_probabilites[type_of_model] = {}
		probabilities_for_promoters_of_interacting_enhancers[type_of_model] = {}

		for domain_atr in ["outside_domain", "within_domain", "unsplit"]:
			sensitivity_match_MAP[type_of_model][domain_atr] = {}


	for type_of_model in type_of_models:

		for chrom in chroms_to_infer:
		
			# classifiers_clean takes value of mode from config_variables on it's own

			i_s_f, j_s_f = MAP_.positive_negative_interactions_for_MAP(chrom)

			if type_of_model in ["dist", "MOG_dist"]:

				match_MAP[type_of_model][chrom], infered_elements[type_of_model][chrom], MAP_probabilites[type_of_model][chrom], probabilities_for_promoters_of_interacting_enhancers[type_of_model][chrom] = MAP_.MAP(posterior["positive_interactions"][type_of_model], posterior["negative_interactions"][type_of_model], chrom, i_s_f, j_s_f)
			
			elif type_of_model in ["correl", "correl_dist", "MOG_correl_dist"]:
				match_MAP[type_of_model][chrom], infered_elements[type_of_model][chrom], MAP_probabilites[type_of_model][chrom] = {}, {}, {}
				for option_ in selected_combinations:
					comb = "_".join([dict_option_[el] for el in option_])
					print comb
					
					match_MAP[type_of_model][chrom][comb], infered_elements[type_of_model][chrom][comb], MAP_probabilites[type_of_model][chrom][comb], probabilities_for_promoters_of_interacting_enhancers[type_of_model][chrom] = MAP_.MAP(posterior["positive_interactions"][type_of_model][comb], posterior["negative_interactions"][type_of_model][comb], chrom, i_s_f, j_s_f)


	for type_of_model in type_of_models:

		if type_of_model in ["dist", "MOG_dist"]:
			 
			match_of_total_number_of_predictions = np.sum([np.sum(match_MAP[type_of_model][chrom_]) for chrom_ in chroms_to_infer])
			total_number_of_predictions = float(np.sum([len(match_MAP[type_of_model][chrom_]) for chrom_ in chroms_to_infer]))
			sensitivity_match_MAP[type_of_model]["unsplit"] = match_of_total_number_of_predictions/total_number_of_predictions
			for domain_atr in ["outside_domain", "within_domain"]:
				interacting_non_intracting_mask = get_MAPS_for_domain_non_domain_interacting_enhancers(domain_atr)

				match_of_total_number_of_predictions = np.sum([np.sum(match_MAP[type_of_model][chrom_][interacting_non_intracting_mask[chrom_]]) for chrom_ in chroms_to_infer])
				total_number_of_predictions = float(np.sum([len(match_MAP[type_of_model][chrom_][interacting_non_intracting_mask[chrom_]]) for chrom_ in chroms_to_infer]))	
				sensitivity_match_MAP[type_of_model][domain_atr] = match_of_total_number_of_predictions/total_number_of_predictions

		elif type_of_model in ["correl", "correl_dist", "MOG_correl_dist"]:
			for option_ in selected_combinations:
				comb = "_".join([dict_option_[el] for el in option_])
				 
				match_of_total_number_of_predictions = np.sum([np.sum(match_MAP[type_of_model][chrom_][comb]) for chrom_ in chroms_to_infer])
				total_number_of_predictions = float(np.sum([len(match_MAP[type_of_model][chrom_][comb]) for chrom_ in chroms_to_infer]))	
				sensitivity_match_MAP[type_of_model]["unsplit"][comb] = match_of_total_number_of_predictions/total_number_of_predictions

				for domain_atr in ["outside_domain", "within_domain"]:
					interacting_non_intracting_mask = get_MAPS_for_domain_non_domain_interacting_enhancers(domain_atr)

					match_of_total_number_of_predictions = np.sum([np.sum(match_MAP[type_of_model][chrom_][comb][interacting_non_intracting_mask[chrom_]]) for chrom_ in chroms_to_infer])
					total_number_of_predictions = float(np.sum([len(match_MAP[type_of_model][chrom_][comb][interacting_non_intracting_mask[chrom_]]) for chrom_ in chroms_to_infer]))	
					sensitivity_match_MAP[type_of_model][domain_atr][comb] = match_of_total_number_of_predictions/total_number_of_predictions


	#f1.close()
	return match_MAP, sensitivity_match_MAP, MAP_probabilites, infered_elements, probabilities_for_promoters_of_interacting_enhancers
	#return MAP_probabilites, infered_elements, match_MAP, sensitivity_match_MAP

