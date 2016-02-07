def interactions_producer_filter(generator_mode, domain, max_path, TSS_or_intra_genic_for_domain_filter, gene_or_TSS_mode):

	
	import os.path
	import config_variables
	np = config_variables.np
	temp_output = config_variables.temp_output

	if generator_mode == "filter_independent_generator":
		filtered_promoters = config_variables.filtered_promoters
		filtered_enhancers = config_variables.filtered_enhancers
		pro_survived, enh_survived = np.where(np.ones_like(filtered_promoters))[0], np.where(np.ones_like(filtered_enhancers))[0]
		f_name = temp_output + 'clean_true_interactions_all_without_prom_prom_filter_independent_d{0}_{1}_max_path_{2}'.format(domain, TSS_or_intra_genic_for_domain_filter,max_path)


	elif generator_mode == "filter_correl_dependent_generator":
		Pol_2_correl_filtered_promoters = config_variables.Pol_2_correl_filtered_promoters
		Pol_2_correl_filtered_enhancers = config_variables.Pol_2_correl_filtered_enhancers
		filter_pro_ = config_variables.filter_value
		filter_enh_ = config_variables.filter_enh

		pro_survived, enh_survived = np.where(Pol_2_correl_filtered_promoters)[0], np.where(Pol_2_correl_filtered_enhancers)[0]
		f_name = temp_output + 'clean_true_interactions_all_without_prom_prom_{0}_{1}_d{2}_{3}_max_path_{4}'.format(filter_pro_, filter_enh_, domain, TSS_or_intra_genic_for_domain_filter, max_path)
		
	elif generator_mode == "filter_dependent_generator":
		filter_pro_ = config_variables.filter_value
		filter_enh_ = config_variables.filter_enh
		count_f_p_ = config_variables.count_f_p
		count_f_e_ = config_variables.count_f_e
		filter_option = config_variables.filt_option 
		ER_pro_filtered_ = config_variables.ER_pro_filtered_
		pro_survived, enh_survived, f_name = np.where(filtered_promoters)[0], np.where(filtered_enhancers)[0]
		f_name = temp_output + 'clean_true_interactions_all_without_prom_prom_{0}_{1}_cop{2}_coe{3}_opt{4}_ov{5}_d{6}_{7}_max_path_{8}'.format(filter_pro_, filter_enh_, count_f_p_, count_f_e_, filter_option, ER_pro_filtered_, domain, TSS_or_intra_genic_for_domain_filter, max_path)

	upstream_validation = config_variables.upstream
	downstream_validation = config_variables.downstream
	upstream = config_variables.upstream

	
	if gene_or_TSS_mode == "TSS_MODE":
		f_name = f_name + "_{0}_{1}_TSS".format(upstream_validation, downstream_validation)
	elif gene_or_TSS_mode == "GENE_MODE":
		f_name = f_name + "_{0}_{1}_Gene".format(upstream_validation, downstream_validation)

	generate = not(os.path.isfile(f_name))

	if generate:
		import interaction_generator_clean
		interactions_of_path = interaction_generator_clean.generator(pro_survived, enh_survived, domain, max_path) # generator also depends on the TSS_or_intra_genic_for_domain_filter if domain == True takes the value from config_varables
		np.savetxt(f_name, interactions_of_path, fmt = '%s', delimiter = '\t')

	return f_name




