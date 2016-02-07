def executor():

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
	
	f1 = open(results_folder + "MAP_performance_{0}_{1}_p_{2}_i_{3}_d{4}".format(mode, one_sided_or_two_sided, "-".join(chroms_in_prior[[0,-1]]), ",".join(chroms_to_infer), domain), "w")


	stuff = [1,2,3,4]
	dict_option_ = {0: 'Pol2_2012-03', 1: 'Pol2',  2: 'H2AZ', 3: 'ER', 4: 'H3K4me3', 5: '2012-03_RNA', 6: 'RNA'}
	combinations = []
	for L in range(0, len(stuff)+1):
	  for subset in itertools.combinations(stuff, L):
		if len(subset): combinations += [list(subset)]

	combinations = np.array(combinations)

	sensitivity_match_MAP = {}
	sensitivity_match_MAP["correl_dist"] = {}
	sensitivity_match_MAP["correl"] = {}
	sensitivity_match_MAP["dist"] = {}

	for domain_atr in ["outside_domain", "within_domain"]:
		sensitivity_match_MAP["correl_dist"][domain_atr] = {}
		sensitivity_match_MAP["correl"][domain_atr] = {}
		sensitivity_match_MAP["dist"][domain_atr] = {}

	total_number_of_interacting_enhancers = {}

	for option_ in combinations:
	
		comb = ",".join([dict_option_[el] for el in option_])
		if option_ == combinations[-1]: comb = "All"
		print comb
		f1.write(comb + "\n"   )

		match_MAP_correl_dist = {}
		infered_elements_correl_dist = {}
		MAP_probabilites_correl_dist = {}

		match_MAP_correl = {}
		infered_elements_correl = {}
		MAP_probabilites_correl = {}	

		match_MAP_dist = {}
		infered_elements_dist = {}
		MAP_probabilites_dist = {}

		posterior_correl_dist_true, posterior_correl_dist_false = classifiers_clean.posterior_producer([0], option_)		
		posterior_correl_true, posterior_correl_false = classifiers_clean.posterior_producer([], option_)
		posterior_dist_true, posterior_dist_false = classifiers_clean.posterior_producer([0], [])

		for chrom in chroms_to_infer:
		
			# classifiers_clean takes value of mode from config_variables on it's own

			i_s_f, j_s_f = MAP_.positive_negative_interactions_for_MAP(chrom)	
			match_MAP_correl_dist[chrom], infered_elements_correl_dist[chrom], MAP_probabilites_correl_dist[chrom] = MAP_.MAP(posterior_correl_dist_true, posterior_correl_dist_false, chrom, i_s_f, j_s_f)
			match_MAP_correl[chrom], infered_elements_correl[chrom], MAP_probabilites_correl[chrom] = MAP_.MAP(posterior_correl_true, posterior_correl_false, chrom, i_s_f, j_s_f)
			match_MAP_dist[chrom], infered_elements_dist[chrom], MAP_probabilites_dist[chrom] = MAP_.MAP(posterior_dist_true, posterior_dist_false, chrom, i_s_f, j_s_f)

		
		match_MAP_correl_dist_total = list(itertools.chain.from_iterable([match_MAP_correl_dist[chrom_] for chrom_ in chroms_to_infer]))
		match_MAP_correl_total = list(itertools.chain.from_iterable([match_MAP_correl[chrom_] for chrom_ in chroms_to_infer]))
		match_MAP_dist_total = list(itertools.chain.from_iterable([match_MAP_dist[chrom_] for chrom_ in chroms_to_infer]))

		f1.write("correl_dist:" + str(np.sum(match_MAP_correl_dist_total)/float(len(match_MAP_correl_dist_total))) + "\n"+ "\n")
		f1.write("dist:" + str(np.sum(match_MAP_dist_total)/float(len(match_MAP_dist_total))) + "\n"+ "\n")
		f1.write("correl:" + str(np.sum(match_MAP_correl_total)/float(len(match_MAP_correl_total))) + "\n" + "\n")
		print "correl_dist:", np.sum(match_MAP_correl_dist_total)/float(len(match_MAP_correl_dist_total))
		print "dist:", np.sum(match_MAP_dist_total)/float(len(match_MAP_dist_total))
		print "correl:", np.sum(match_MAP_correl_total)/float(len(match_MAP_correl_total))
		sensitivity_match_MAP["correl_dist"][",".join(np.array(option_, str))] = np.sum(match_MAP_correl_dist_total)/float(len(match_MAP_correl_dist_total))
		sensitivity_match_MAP["correl"][",".join(np.array(option_, str))] = np.sum(match_MAP_correl_total)/float(len(match_MAP_correl_total))
		sensitivity_match_MAP["dist"][",".join(np.array(option_, str))] = np.sum(match_MAP_dist_total)/float(len(match_MAP_dist_total))

		def get_MAPS_for_domain_non_domain_interacting_enhancers(domain_atr):

			chrom_interacting_enhancers_pro = config_variables.chrom_interacting_enhancers_pro

			from  prepare_interactions_clean import un_string
			interacting_non_intracting_mask = {}

			total_number_of_interacting_enhancers = 0

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

				total_number_of_interacting_enhancers += np.sum(interacting_non_intracting_mask[chrom___])
			

			match_MAP_correl_dist_total = list(itertools.chain.from_iterable([match_MAP_correl_dist[chrom_][interacting_non_intracting_mask[chrom_]] for chrom_ in chroms_to_infer]))
			match_MAP_correl_total = list(itertools.chain.from_iterable([match_MAP_correl[chrom_][interacting_non_intracting_mask[chrom_]] for chrom_ in chroms_to_infer]))
			match_MAP_dist_total = list(itertools.chain.from_iterable([match_MAP_dist[chrom_][interacting_non_intracting_mask[chrom_]] for chrom_ in chroms_to_infer]))

			sensitivity_match_MAP_correl_dist = np.sum(match_MAP_correl_dist_total)/float(len(match_MAP_correl_dist_total))
			sensitivity_match_MAP_correl = np.sum(match_MAP_correl_total)/float(len(match_MAP_correl_total))
			sensitivity_match_MAP_dist = np.sum(match_MAP_dist_total)/float(len(match_MAP_dist_total))

			return sensitivity_match_MAP_correl_dist, sensitivity_match_MAP_correl, sensitivity_match_MAP_dist, total_number_of_interacting_enhancers
	
		for domain_atr in ["outside_domain", "within_domain"]:
			[sensitivity_match_MAP["correl_dist"][domain_atr][",".join(np.array(option_, str))], 
			sensitivity_match_MAP["correl"][domain_atr][",".join(np.array(option_, str))],	
			sensitivity_match_MAP["dist"][domain_atr][",".join(np.array(option_, str))],
			total_number_of_interacting_enhancers[domain_atr]] = get_MAPS_for_domain_non_domain_interacting_enhancers(domain_atr)

	config_variables.total_number_of_interacting_enhancers = total_number_of_interacting_enhancers

	f1.close()
	return MAP_probabilites_correl_dist, infered_elements_correl_dist, match_MAP_correl_dist, sensitivity_match_MAP

