	name_of_saved_priors_file = config_variables.temp_output + 'priors_average_{0}_{1}'.format(one_sided_or_two_sided, config_variables.mode_of_code)

	name_of_saved_priors_file = name_of_saved_priors_file + "{0}_{1}_{2}_{3}_{4}".format(config_variables.number_of_bins[0], config_variables.number_of_bins[1], config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation: name_of_saved_priors_file += "_TSS"
	else: name_of_saved_priors_file += "_GENE"

	def does_file_exists(extra_string):
		name_of_saved_priors_file_ = "{0}_{1}".format(name_of_saved_priors_file, extra_string)

		name_of_saved_priors_file_exist = os.path.exists(name_of_saved_priors_file_)
			
		return pair_wise_prob_exist, name_of_saved_priors_file_


	for data_set_name in dataset_names_option:
		name_of_saved_priors_file_ = name_of_saved_priors_file + data_set_name
		does_file_exists(name_of_saved_priors_file_)

 in ["distance", "correlation"]
	def saver(name_of_saved_priors_file_, attribute_of_interaction, data_set):

		name_of_saved_priors_file_ = name_of_saved_priors_file_ + attribute_of_interaction
		name_of_saved_priors_file_exist = os.path.exists(name_of_saved_priors_file_)

		if attribute_of_interaction == "correlation":
			np.save(name_of_saved_priors_file_ + "{0}_freq".format(data_set), prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_frequencies"]
			np.save(name_of_saved_priors_file_ + "{0}_bins".format(data_set), prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_bins"]
		else:
			np.save(name_of_saved_priors_file_ + "_freq", prior_elements[mode][classification_of_interactions]["distance"]["prior_frequencies"]
			np.save(name_of_saved_priors_file_ + "_bins", prior_elements[mode][classification_of_interactions]["distance"]["prior_bins"]


	for data_set_name in dataset_names_option:

		saver(name_of_saved_priors_file_, attribute_of_interaction, data_set)

	open_atr = [('distance', 'dist')] + [("correlation", dataset) for dataset in dataset_names_option]

	for data_set_name in open_atr:
		name_of_saved_freq = name_of_saved_priors_file_ + "{0}_freq".format(data_set)
		if os.path.exists(name_of_saved_freq): prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_frequencies"] = np.load(name_of_saved_freq)
