def execute(sensitivity_match_MAP, number_of_interacting_enhancers_, option_to_plot = "ALL", type_of_models = [], posterior_MOG = None, kappa_0=1.0, mu_0=1.0, alpha_0=1.0, Beta_0=1.0, number_of_samples=10, burn_in=0):

	#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

	from matplotlib.backends.backend_pdf import PdfPages
	import config_variables
	import numpy as np
	results_folder = config_variables.results_folder

	name_of_output_FDR_file = results_folder + 'ATEST_FDR_file_{0}_{1}_{4}_{2}_{3}_average_PolII'.format(config_variables.chroms_in_prior[0], config_variables.chroms_to_infer[0], config_variables.one_sided_or_two_sided, config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_FDR_file	+= "_{0}".format("_".join(type_of_models))

	name_of_output_FDR_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	#name_of_output_FDR_file += "_{0}_{1}".format("_".join(np.array(number_of_samples,str)), "_".join(np.array(burn_in, str)))

	if config_variables.disentagled_features_validation: 

		name_of_output_FDR_file += "_TSS" 
	else:
		name_of_output_FDR_file += "_GENE"
	if option_to_plot == "ALL": name_of_output_FDR_file += "_ALL"

	if "MOG_correl_dist" in type_of_models or "MOG_dist" in type_of_models: name_of_output_FDR_file += "_MOG_{0}_{1}_{2}_{3}_{4}_{5}__test".format(kappa_0, mu_0, alpha_0, Beta_0, "_".join(np.array(number_of_samples, str)), "_".join(np.array(burn_in, str)))

	pdf = PdfPages(name_of_output_FDR_file + ".pdf")



	np = config_variables.np
	negative_interactions = config_variables.negative_interactions

	if config_variables.FDR_mode:
		
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes = results_folder + "file_with_FDRs_{0}_{1}_smo_{2}_{3}".format(config_variables.chroms_in_prior[0], config_variables.chroms_to_infer[0], config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

		if config_variables.disentagled_features_validation:
			name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_TSS"
		else:
			name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_GENE"

		if "MOG_correl_dist" in type_of_models or "MOG_dist" in type_of_models: 
			name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_{0}".format("_".join(type_of_models))
			name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_MOG_{0}_{1}_{2}_{3}_{4}_{5}__test".format(kappa_0, mu_0, alpha_0, Beta_0, "_".join(np.array(number_of_samples, str)), "_".join(np.array(burn_in, str)))


		output = open(name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes, "w")
		output.write("#" + "\t".join(["Data_set", "FDR", "Precision", "True_links_above_FDR", "Threshold"]) + "\n")	




	save_file_with_cutoffs = True
	if save_file_with_cutoffs:
		
		name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes = results_folder + "file_with_cutoffs_{0}_{1}_smo_{2}_{3}".format(config_variables.chroms_in_prior[0], config_variables.chroms_to_infer[0], config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

		name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

		if config_variables.disentagled_features_validation:
			name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes += "_TSS"
		else:
			name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes += "_GENE"

		if "MOG_correl_dist" in type_of_models or "MOG_dist" in type_of_models: 
			name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes	+= "_{0}".format("_".join(type_of_models))
			name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes += "_MOG_{0}_{1}_{2}_{3}_{4}_{5}__test".format(kappa_0, mu_0, alpha_0, Beta_0, "_".join(np.array(number_of_samples, str)), "_".join(np.array(burn_in, str)))


		output_2 = open(name_of_output_file_with_cutoffs_estimated_on_odd_even_chromosomes, "w")
		output_2.write("#" + "\t".join(["Data_set", "Mode", "Cutoff10%", "Cutoff20%", "Cutoff30%"]) + "\n")	



	TSS_coordinates = config_variables.negative_interactions.extract_TSS_coordinates(config_variables.upstream)

	def positive_negative_interactions_for_MAP(chrom):
		
		indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

		if mode == "promoter_enhancer_interactions":

			false_inter_pro = negative_interactions.chrom_specific_negative_interactions(chrom, mode)

			i_s_f, j_s_f = false_inter_pro[:,0] + total_p, false_inter_pro[:,1] + total_e

		if mode == "enhancer_enhancer_interactions":

			false_inter_enh = negative_interactions.chrom_specific_negative_interactions(chrom, mode)
	
			i_s_f, j_s_f = false_inter_enh[:,0] + total_e, false_inter_enh[:,1] + total_e

		return i_s_f, j_s_f


	def filter_interactions_in_domain(posterior_t, posterior_f, chrom, domain, invert_domain):
		enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)

		i_s_f, j_s_f = positive_negative_interactions_for_MAP(chrom)
		
		length_chr = len(indexes_p) + len(indexes_e)
		interaction_matrix = np.zeros((length_chr, length_chr))
		posterior_t, posterior_f = posterior_t[chrom], posterior_f[chrom]
		
		if domain:
			if config_variables.TSS_or_intra_genic_for_domain_filter == "Intra_genic": coords_pro_domain = pro_coordinates[indexes_p]
			elif config_variables.TSS_or_intra_genic_for_domain_filter == "TSS_only": coords_pro_domain = np.column_stack((TSS_coordinates[indexes_p]-1, TSS_coordinates[indexes_p]+1))
			domain_matrix = interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'left', True)
			domain_matrix = domain_matrix + interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'right', True)
			if invert_domain: domain_matrix = np.invert(domain_matrix)

		else:
			domain_matrix = True

			
		if mode == "promoter_enhancer_interactions":

			#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
			#true_inter_pro = un_string(chr_interactions_dict_pro_enh[chrom][:, :2]).astype(int)
			#i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

			if config_variables.disentagled_features_validation: 
				chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
			else:
				chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

			true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)
			i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

			interaction_matrix[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = posterior_t
			interaction_matrix[i_s_f - total_p, j_s_f + len(indexes_p) - total_e] = posterior_f

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = True

			true_pro_enh_inter_filtered = interacting_mask * domain_matrix

			print np.sum(true_pro_enh_inter_filtered)

			chrom_posterior_t_filtered = interaction_matrix[true_pro_enh_inter_filtered]

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_f - total_p, j_s_f + len(indexes_p) - total_e] = True

			false_pro_enh_inter_filtered = interacting_mask * domain_matrix

			print np.sum(false_pro_enh_inter_filtered)

			chrom_posterior_f_filtered = interaction_matrix[false_pro_enh_inter_filtered]

			return	chrom_posterior_t_filtered, chrom_posterior_f_filtered

		if mode == "enhancer_enhancer_interactions":

			chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh
			true_inter_enh = un_string(chr_interactions_dict_enh_enh[chrom][:, :2]).astype(int)
			i_s_t, j_s_t = true_inter_enh[:,0], true_inter_enh[:,1]

			
			interaction_matrix[i_s_t + len(indexes_p) - total_e, j_s_t + len(indexes_p) - total_e] = posterior_t
			interaction_matrix[i_s_f + len(indexes_p) - total_e, j_s_f + len(indexes_p) - total_e] = posterior_f
			interaction_matrix[j_s_t + len(indexes_p) - total_e, i_s_t + len(indexes_p) - total_e] = posterior_t # transpose to create a full matrix
			interaction_matrix[j_s_f + len(indexes_p) - total_e, i_s_f + len(indexes_p) - total_e] = posterior_f # transpose to create a full matrix


			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_t + len(indexes_p) - total_e, j_s_t + len(indexes_p) - total_e] = True

			true_enh_enh_inter_filtered = interacting_mask * domain_matrix
			chrom_posterior_t_filtered = interaction_matrix[true_enh_enh_inter_filtered]

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_f + len(indexes_p) - total_e, j_s_f + len(indexes_p) - total_e] = True

			false_enh_enh_inter_filtered = interacting_mask * domain_matrix
			chrom_posterior_f_filtered = interaction_matrix[false_enh_enh_inter_filtered]
			
			return chrom_posterior_t_filtered, chrom_posterior_f_filtered


	

	from  prepare_interactions_clean import un_string
	normalised = False
	import interacting_domain
	import itertools

	def domain_filter(inpu, domain, invert_domain):

		posterior_t, posterior_f = inpu 
		chrom_posterior_t_filtered, chrom_posterior_f_filtered = {}, {}
		for chrom__ in chroms_to_infer:
			chrom_posterior_t_filtered[chrom__], chrom_posterior_f_filtered[chrom__] = filter_interactions_in_domain(posterior_t, posterior_f, chrom__, domain, invert_domain)

		posterior_t_filtered = np.array(list(itertools.chain.from_iterable([chrom_posterior_t_filtered[chrom_] for chrom_ in chroms_to_infer])))
		posterior_f_filtered = np.array(list(itertools.chain.from_iterable([chrom_posterior_f_filtered[chrom_] for chrom_ in chroms_to_infer])))
	
		return posterior_t_filtered, posterior_f_filtered
	#----------------------------------------------------------------------------------------------------------------------------------------------------------------------	


	import itertools
	import matplotlib.pyplot as plt
	import config_variables
	np = config_variables.np	
	classificator_elements = config_variables.classificator_elements
	classifiers_clean = config_variables.classifiers_clean
	filter_values = config_variables.filter_values
	datasets_names = config_variables.datasets_names
	chroms_to_infer = config_variables.chroms_to_infer
	mode = 	config_variables.mode

	
	#dict_option = {0: 'Pol2_2012-03', 1: 'Pol2',  2: 'H2AZ', 3: 'ER', 4: 'H3K4me3', 5: '2012-03_RNA', 6: 'RNA'}
	dict_option = dict(zip(range(len(datasets_names)), datasets_names))

	def calculate_single_ROC_best_True_sensitivity(probabilities_true, probabilities_false, length_of_positives, length_of_negatives, percent_1, percent_2, percent_3, thresh = False, give_indexes_for_thresholds = False):

		_True_positives_of_threshold = []
		_False_positives_of_threshold = []

		sorted_prob_true = np.sort(probabilities_true)
		sorted_prob_false = np.sort(probabilities_false)

		sorted_thresholds = np.sort(np.unique(np.r_[probabilities_true, probabilities_false]))
		sorted_thresholds = np.unique(np.r_[sorted_thresholds, np.max(sorted_thresholds)*1.01])

		len_prob_true = len(probabilities_true)
		len_prob_false = len(probabilities_false)

		print 'len prob: ', len_prob_true, len_prob_false

		_True_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_true, sorted_thresholds)[0][::-1])

		_False_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_false, sorted_thresholds)[0][::-1])

		Precision = np.array(_True_positives_of_threshold, dtype = float)/(np.array(_True_positives_of_threshold, dtype = float) + np.array(_False_positives_of_threshold, dtype = float))
	
		True_positive_Rate = np.array(_True_positives_of_threshold)/float(len_prob_true)#float(length_of_positives)#
		
		False_positive_Rate = np.array(_False_positives_of_threshold)/float(len_prob_false)#float(length_of_negatives)#


		if give_indexes_for_thresholds:

			threshold_1, threshold_2, threshold_3 = percent_1, percent_2, percent_3
			index_100_first_occurance =	np.where(sorted_thresholds[::-1] <= threshold_1)[0][0] - 1
			index_200_first_occurance = np.where(sorted_thresholds[::-1] <= threshold_2)[0][0] - 1
			index_300_first_occurance = np.where(sorted_thresholds[::-1] <= threshold_3)[0][0] - 1



			threshold_1, threshold_2, threshold_3 = [], [], []
			return True_positive_Rate, False_positive_Rate, Precision, index_100_first_occurance, index_200_first_occurance, index_300_first_occurance, threshold_1, threshold_2, threshold_3

		else:

			index_100_first_occurance = np.where(True_positive_Rate >= percent_1)[0][0]
			index_200_first_occurance = np.where(True_positive_Rate >= percent_2)[0][0]
			index_300_first_occurance = np.where(True_positive_Rate >= percent_3)[0][0]

			if config_variables.FDR_mode and thresh:

				for FDR in config_variables.FDR:

					print "Precision", Precision[:40]
					print "TPR", True_positive_Rate[:40]
					try: 
						left_most = np.where(Precision >= 1-FDR)[0][-1]
						to_save = "\t".join(np.array([config_variables.comb, config_variables.dist_or_correl_attribute, FDR, Precision[left_most], _True_positives_of_threshold[left_most], sorted_thresholds[::-1][left_most + 2]]).astype("|S30")) # should be sorted_thresholds[::-1][left_most + 1] but +2 corrects rounding errors caused by saving into a txt file. 
						output.write(to_save + "\n")
						#var = raw_input("Please enter something (pause): ")
					except:
						to_save = "\t".join([config_variables.comb, config_variables.dist_or_correl_attribute, str(FDR), "N\A", "N\A", "N\A"])
						output.write(to_save + "\n")


		if thresh:
			threshold_1 = sorted_thresholds[::-1][index_100_first_occurance + 1]
			threshold_2 = sorted_thresholds[::-1][index_200_first_occurance + 1]
			threshold_3 = sorted_thresholds[::-1][index_300_first_occurance + 1]
			return True_positive_Rate, False_positive_Rate, Precision, index_100_first_occurance, index_200_first_occurance, index_300_first_occurance, threshold_1, threshold_2, threshold_3
	
		#print 'number of thresholds', len(True_positive_Rate), len(False_positive_Rate)

		#return True_positive_Rate, False_positive_Rate, Precision, index_100_first_occurance, index_200_first_occurance, index_300_first_occurance


	from pylab import rcParams
	#rcParams['figure.figsize'] = 20, 8

	#stuff = [0, 1, 2, 3, 4]

	import selected_combinations as sel
	combinations, selected_combinations = sel.selected_combinations(option_to_plot)

	filter_values_ = filter_values[[0]]


	#http://stackoverflow.com/questions/14270391/python-matplotlib-multiple-bars - alternatively
	#http://matplotlib.org/examples/pylab_examples/subplots_demo.html


	percent_1_, percent_2_, percent_3_ = 0.1, 0.2, 0.3

	dict_option_ = dict_option
	datasets_names_ = datasets_names

	if "correl_dist" in type_of_models: threshold_1_dist_correl, threshold_2_dist_correl, threshold_3_dist_correl = {}, {}, {}
	if "correl" in type_of_models: threshold_1_correl, threshold_2_correl, threshold_3_correl = {}, {}, {}
	if "dist" in type_of_models: threshold_1_dist, threshold_2_dist, threshold_3_dist = {}, {}, {}
	if "MOG_correl_dist" in type_of_models: threshold_1_dist_correl_MOG, threshold_2_dist_correl_MOG, threshold_3_dist_correl_MOG = {}, {}, {}
	if "MOG_dist" in type_of_models: threshold_1_dist_MOG, threshold_2_dist_MOG, threshold_3_dist_MOG = {}, {}, {}
	

	#total_number_of_interacting_enhancers = config_variables.total_number_of_interacting_enhancers

	posterior_correl_dist_true_unsplit, posterior_correl_dist_false_unsplit = {}, {}
	posterior_correl_true_unsplit, posterior_correl_false_unsplit = {}, {}
	posterior_correl_dist_true_unsplit_MOG, posterior_correl_dist_false_unsplit_MOG = {}, {}



	for domain_atr, domain, invert_domain, thresh, give_indexes_for_thresholds in np.array([[None, False, False, True, False], ["within_domain", True, False, False, True], ["outside_domain", True, True, False, True]])[0:3]:

		if option_to_plot == "ALL":
			plt.rcParams['xtick.labelsize'] = 28
			plt.rc('ytick', labelsize = 28)
			f, ax = plt.subplots(1, len(selected_combinations), sharex=True, sharey=True, figsize=(50,10))
			f.subplots_adjust(left=0.035, bottom=0.1, right=0.975, top=0.925, hspace=0.1, wspace=0.1)
			red_blue_yellow_cyan_marker_size = 12
			red_blue_yellow_cyan_marker_size_legend_box = 19
			legend_box_names_font_size = 20
			size_of_combination_name = 23
			size_of_y_label = 30

			ax[0].set_ylabel('Precision', fontsize = size_of_y_label)



		elif option_to_plot == "SELECTIVE":
			plt.rcParams['xtick.labelsize'] = 28
			plt.rc('ytick', labelsize = 28)
			f, ax = plt.subplots(1, len(selected_combinations), sharex=True, sharey=True, figsize=(20,10))
			f.subplots_adjust(left=0.085, bottom=0.15, right=0.965, top=0.925, hspace=0.1, wspace=0.05)
			red_blue_yellow_cyan_marker_size = 12
			red_blue_yellow_cyan_marker_size_legend_box = 19
			legend_box_names_font_size = 22
			size_of_combination_name = 35
			size_of_y_label = 35

			
			ax[0].set_ylabel('Precision', fontsize = size_of_y_label)

		import matplotlib.lines as mlines
		blue_line = mlines.Line2D([], [], color='blue', marker='^', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data+prior')#'NB')#
		yellow_line = mlines.Line2D([], [], color='darkviolet', marker='s', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='prior')
		red_line = mlines.Line2D([], [], color='red', marker='o', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data')
		pink_line = mlines.Line2D([], [], color='cyan', marker='*', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='LVA data+prior')
		green_line = mlines.Line2D([], [], color='green', marker="v", markersize = red_blue_yellow_cyan_marker_size_legend_box, label='LVA prior')
		#handles, labels = ax[0].get_legend_handles_labels()
		#ax[0].legend(handles, labels, fontsize = legend_box_names_font_size)





		for index_filt_val, filter_value in enumerate(filter_values_):


			if "dist" in type_of_models: 
				if not(domain_atr): posterior_dist_true_unsplit, posterior_dist_false_unsplit = posterior_MOG["positive_interactions"]["dist"], posterior_MOG["negative_interactions"]["dist"]#classifiers_clean.posterior_producer([0], [], total_posterior = False)

				posterior_dist_true, posterior_dist_false = domain_filter((posterior_dist_true_unsplit, posterior_dist_false_unsplit), domain, invert_domain)

			if "MOG_dist" in type_of_models: 
				if not(domain_atr): posterior_dist_true_unsplit_MOG, posterior_dist_false_unsplit_MOG = posterior_MOG["positive_interactions"]["MOG_dist"], posterior_MOG["negative_interactions"]["MOG_dist"]

				posterior_dist_true_MOG, posterior_dist_false_MOG = domain_filter((posterior_dist_true_unsplit_MOG, posterior_dist_false_unsplit_MOG), domain, invert_domain)
			
			#posterior_dist_true, posterior_dist_false = np.array(posterior_dist_true), np.array(posterior_dist_false)
			#posterior_dist_true, posterior_dist_false = posterior_dist_true[posterior_dist_true <> 1.], posterior_dist_false[posterior_dist_false <> 1.]

			for index_opt, option_ in zip(range(len(selected_combinations)), selected_combinations):#[0:5]:#[1:3]

					
				comb = ",".join([dict_option_[el] for el in option_])
				comb_MOG = "_".join([dict_option[el] for el in option_])


				if option_ == combinations[-1]: comb = "All"
				if config_variables.FDR_mode: config_variables.comb = comb
				print comb

				if domain_atr <> None:
					#number_of_interacting_enhancers = total_number_of_interacting_enhancers[domain_atr]
					sensitivity_dist_correl = sensitivity_match_MAP["correl_dist"][domain_atr][comb_MOG]
					sensitivity_correl = sensitivity_match_MAP["correl"][domain_atr][comb_MOG]
					sensitivity_dist = sensitivity_match_MAP["dist"][domain_atr]#[",".join(np.array(option_, str))]
					
					if "MOG_dist" in type_of_models: sensitivity_dist_MOG = sensitivity_match_MAP["MOG_dist"][domain_atr]
					if "MOG_correl_dist" in type_of_models: sensitivity_dist_correl_MOG = sensitivity_match_MAP["MOG_correl_dist"][domain_atr][comb_MOG]

					# sensitivity_MOG

				else:
					#number_of_interacting_enhancers = number_of_interacting_enhancers_
					sensitivity_dist_correl = sensitivity_match_MAP["correl_dist"]["unsplit"][comb_MOG]
					sensitivity_correl = sensitivity_match_MAP["correl"]["unsplit"][comb_MOG]
					sensitivity_dist = sensitivity_match_MAP["dist"]["unsplit"]#[",".join(np.array(option_, str))]
					
					if "MOG_dist" in type_of_models: sensitivity_dist_MOG = sensitivity_match_MAP["MOG_dist"]["unsplit"]
					if "MOG_correl_dist" in type_of_models: sensitivity_dist_correl_MOG = sensitivity_match_MAP["MOG_correl_dist"]["unsplit"][comb_MOG]

				full_len = sum([len(classificator_elements[-1.][mode]["positive_interactions"]["distance"]["probabilities_of_being_positive_interactions"]["posterior_component_values"][chrom_]) for chrom_ in chroms_to_infer])
				new_len = sum([len(classificator_elements[filter_value][mode]["positive_interactions"]["distance"]["probabilities_of_being_positive_interactions"]["posterior_component_values"][chrom_]) for chrom_ in chroms_to_infer])

				length_of_positives_pro = full_len
				length_of_negatives_pro = sum([len(classificator_elements[-1.][mode]["negative_interactions"]["distance"]["probabilities_of_being_positive_interactions"]["posterior_component_values"][chrom_]) for chrom_ in chroms_to_infer])
				#full_len = len(true_interactions_dist_correl_pairwise_prob_pro_filter_comb[-1.])
				#new_len = len(true_interactions_dist_correl_pairwise_prob_pro_filter_comb[filter_value])

				per = float(full_len)/float(new_len)
				percent_1 = percent_1_ * per
				percent_2 = percent_2_ * per
				percent_3 = percent_3_ * per


				
				if not(domain_atr): posterior_correl_dist_true_unsplit[comb], posterior_correl_dist_false_unsplit[comb] = classifiers_clean.posterior_producer([0], option_, total_posterior = False)

				posterior_correl_dist_true, posterior_correl_dist_false = domain_filter((posterior_correl_dist_true_unsplit[comb], posterior_correl_dist_false_unsplit[comb]), domain, invert_domain)

				length_of_positives_pro, length_of_negatives_pro = len(posterior_correl_dist_true), len(posterior_correl_dist_false) #domain adjusted
				#length_of_positives_pro, length_of_negatives_pro = len(posterior_correl_dist_true_unsplit), len(posterior_correl_dist_false_unsplit) #domain unadjusted

				if "correl" in type_of_models: 
					if not(domain_atr): posterior_correl_true_unsplit[comb], posterior_correl_false_unsplit[comb] = classifiers_clean.posterior_producer([], option_, total_posterior = False)

					posterior_correl_true, posterior_correl_false = domain_filter((posterior_correl_true_unsplit[comb], posterior_correl_false_unsplit[comb]), domain, invert_domain)


				if "MOG_correl_dist" in type_of_models: 
					if not(domain_atr): posterior_correl_dist_true_unsplit_MOG[comb], posterior_correl_dist_false_unsplit_MOG[comb] = posterior_MOG["positive_interactions"]["MOG_correl_dist"][comb_MOG], posterior_MOG["negative_interactions"]["MOG_correl_dist"][comb_MOG]
					posterior_correl_dist_true_MOG, posterior_correl_dist_false_MOG = domain_filter((posterior_correl_dist_true_unsplit_MOG[comb], posterior_correl_dist_false_unsplit_MOG[comb]), domain, invert_domain)



				if "correl_dist" in type_of_models: 
					if give_indexes_for_thresholds: 
						percent_1, percent_2, percent_3 = threshold_1_dist_correl[comb], threshold_2_dist_correl[comb], threshold_3_dist_correl[comb]

					config_variables.dist_or_correl_attribute = "distance_correl"
					True_positive_Rate_dist_correl_pro, False_positive_Rate_dist_correl_pro, precision_dist_correl_pro, index_100_dist_correl_pro, index_200_dist_correl_pro, index_300_dist_correl_pro, threshold_1_dist_correl_, threshold_2_dist_correl_, threshold_3_dist_correl_ = calculate_single_ROC_best_True_sensitivity(posterior_correl_dist_true, posterior_correl_dist_false, length_of_positives_pro, length_of_negatives_pro, percent_1, percent_2, percent_3, thresh = thresh, give_indexes_for_thresholds = give_indexes_for_thresholds)

				if "dist" in type_of_models: 
					if give_indexes_for_thresholds: 
						percent_1, percent_2, percent_3 = threshold_1_dist[comb], threshold_2_dist[comb], threshold_3_dist[comb]
					config_variables.dist_or_correl_attribute = "distance"
					True_positive_Rate_dist_pro, False_positive_Rate_dist_pro, precision_dist_pro, index_100_dist_pro, index_200_dist_pro, index_300_dist_pro, threshold_1_dist_, threshold_2_dist_, threshold_3_dist_ = calculate_single_ROC_best_True_sensitivity(posterior_dist_true, posterior_dist_false, length_of_positives_pro, length_of_negatives_pro, percent_1, percent_2, percent_3, thresh = thresh, give_indexes_for_thresholds = give_indexes_for_thresholds)


				if "correl" in type_of_models: 
					if give_indexes_for_thresholds: 
						percent_1, percent_2, percent_3 = threshold_1_correl[comb], threshold_2_correl[comb], threshold_3_correl[comb]
					config_variables.dist_or_correl_attribute = "correl"
					True_positive_Rate_correl_pro, False_positive_Rate_correl_pro, precision_correl_pro, index_100_correl_pro, index_200_correl_pro, index_300_correl_pro, threshold_1_correl_, threshold_2_correl_, threshold_3_correl_ = calculate_single_ROC_best_True_sensitivity(posterior_correl_true, posterior_correl_false, length_of_positives_pro, length_of_negatives_pro, percent_1, percent_2, percent_3, thresh = thresh, give_indexes_for_thresholds = give_indexes_for_thresholds)


					#MOG--------------------------
				if "MOG_correl_dist" in type_of_models: 
					if give_indexes_for_thresholds: 
						percent_1, percent_2, percent_3 = threshold_1_dist_correl_MOG[comb], threshold_2_dist_correl_MOG[comb], threshold_3_dist_correl_MOG[comb]
					config_variables.dist_or_correl_attribute = "MOG_correl_dist"
					True_positive_Rate_dist_correl_pro_MOG, False_positive_Rate_dist_correl_pro_MOG, precision_dist_correl_pro_MOG, index_100_dist_correl_pro_MOG, index_200_dist_correl_pro_MOG, index_300_dist_correl_pro_MOG, threshold_1_dist_correl_MOG_, threshold_2_dist_correl_MOG_, threshold_3_dist_correl_MOG_ = calculate_single_ROC_best_True_sensitivity(posterior_correl_dist_true_MOG, posterior_correl_dist_false_MOG, length_of_positives_pro, length_of_negatives_pro, percent_1, percent_2, percent_3, thresh = thresh, give_indexes_for_thresholds = give_indexes_for_thresholds)

				if "MOG_dist" in type_of_models: 
					if give_indexes_for_thresholds: 
						percent_1, percent_2, percent_3 = threshold_1_dist_MOG[comb], threshold_2_dist_MOG[comb], threshold_3_dist_MOG[comb]
					config_variables.dist_or_correl_attribute = "MOG_dist"
					True_positive_Rate_dist_pro_MOG, False_positive_Rate_dist_pro_MOG, precision_dist_pro_MOG, index_100_dist_pro_MOG, index_200_dist_pro_MOG, index_300_dist_pro_MOG, threshold_1_dist_MOG_, threshold_2_dist_MOG_, threshold_3_dist_MOG_ = calculate_single_ROC_best_True_sensitivity(posterior_dist_true_MOG, posterior_dist_false_MOG, length_of_positives_pro, length_of_negatives_pro, percent_1, percent_2, percent_3, thresh = thresh, give_indexes_for_thresholds = give_indexes_for_thresholds)

				
				if thresh: 

					#"Data_set", "Mode", "Cutoff10%", "Cutoff20%", "Cutoff30%"

					if "correl_dist" in type_of_models: 	
						threshold_1_dist_correl[comb], threshold_2_dist_correl[comb], threshold_3_dist_correl[comb] = threshold_1_dist_correl_, threshold_2_dist_correl_, threshold_3_dist_correl_
						 
						to_save_2 = "\t".join(np.array([config_variables.comb, "correl_dist", threshold_1_dist_correl_, threshold_2_dist_correl_, threshold_3_dist_correl_]).astype("|S30"))
						output_2.write(to_save_2 + "\n")

					if "dist" in type_of_models:
						threshold_1_dist[comb], threshold_2_dist[comb], threshold_3_dist[comb] = threshold_1_dist_, threshold_2_dist_, threshold_3_dist_

						to_save_2 = "\t".join(np.array([config_variables.comb, "dist", threshold_1_dist_, threshold_2_dist_, threshold_3_dist_]).astype("|S30"))
						output_2.write(to_save_2 + "\n")

					if "correl" in type_of_models: 
						threshold_1_correl[comb], threshold_2_correl[comb], threshold_3_correl[comb] = threshold_1_correl_, threshold_2_correl_, threshold_3_correl_
						
						to_save_2 = "\t".join(np.array([config_variables.comb, "correl", threshold_1_correl_, threshold_2_correl_, threshold_3_correl_]).astype("|S30"))
						output_2.write(to_save_2 + "\n")

					if "MOG_correl_dist" in type_of_models: 
						threshold_1_dist_correl_MOG[comb], threshold_2_dist_correl_MOG[comb], threshold_3_dist_correl_MOG[comb] = threshold_1_dist_correl_MOG_, threshold_2_dist_correl_MOG_, threshold_3_dist_correl_MOG_
						
						to_save_2 = "\t".join(np.array([config_variables.comb, "MOG_correl_dist", threshold_1_dist_correl_MOG_, threshold_2_dist_correl_MOG_, threshold_3_dist_correl_MOG_]).astype("|S30"))
						output_2.write(to_save_2 + "\n")

					if "MOG_dist" in type_of_models:
						threshold_1_dist_MOG[comb], threshold_2_dist_MOG[comb], threshold_3_dist_MOG[comb] = threshold_1_dist_MOG_, threshold_2_dist_MOG_, threshold_3_dist_MOG_
						

						to_save_2 = "\t".join(np.array([config_variables.comb, "MOG_dist", threshold_1_dist_MOG_, threshold_2_dist_MOG_, threshold_3_dist_MOG_]).astype("|S30"))
						output_2.write(to_save_2 + "\n")


				centres_of_ticks = np.arange(4) + 0.5
				ind = centres_of_ticks

				OX = [0.1,"0.2\n TPR ",0.3, "\nMAP"]
				#ax[index_opt].set_title(comb, fontsize = size_of_combination_name)

				ax[index_opt].vlines(3., 0, 1, colors=u'SlateGray', linestyles=u'dashed')
				ax[index_opt].set_xlim([0., 4.])
				ax[index_opt].set_ylim([0., 1.])

				ax[index_opt].set_xticks(centres_of_ticks)
				ax[index_opt].set_xticklabels(np.array(OX, str))

				ax[index_opt].set_title(comb, fontsize = size_of_combination_name)

				if "correl_dist" in type_of_models: 				
					probabilities_dist_correl = np.r_[precision_dist_correl_pro[[index_100_dist_correl_pro, index_200_dist_correl_pro, index_300_dist_correl_pro]], sensitivity_dist_correl]
					n = np.r_[length_of_positives_pro*True_positive_Rate_dist_correl_pro[[index_100_dist_correl_pro, index_200_dist_correl_pro, index_300_dist_correl_pro]] + False_positive_Rate_dist_correl_pro[[index_100_dist_correl_pro, index_200_dist_correl_pro, index_300_dist_correl_pro]]*length_of_negatives_pro, number_of_interacting_enhancers_]
					yerr = (probabilities_dist_correl*(1-probabilities_dist_correl)/n)**0.5
					ax[index_opt].errorbar(ind, probabilities_dist_correl, yerr= yerr, fmt='^', color="b", alpha=0.5, linewidth=3., markersize=red_blue_yellow_cyan_marker_size)
					ax[index_opt].plot(ind, probabilities_dist_correl, alpha=1.0, color="b", marker= "^", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "NB")#data+distance
				if "dist" in type_of_models:
					probabilities_dist = np.r_[precision_dist_pro[[index_100_dist_pro, index_200_dist_pro, index_300_dist_pro]], sensitivity_dist]
						#n = np.r_[np.array([True_positive_Rate_dist_pro[index_100_dist_pro], True_positive_Rate_dist_pro[index_200_dist_pro], True_positive_Rate_dist_pro[index_300_dist_pro]])*number_of_interacting_enhancers_, number_of_interacting_enhancers_]
					n = np.r_[length_of_positives_pro*True_positive_Rate_dist_pro[[index_100_dist_pro, index_200_dist_pro, index_300_dist_pro]] + False_positive_Rate_dist_pro[[index_100_dist_pro, index_200_dist_pro, index_300_dist_pro]]*length_of_negatives_pro, number_of_interacting_enhancers_]

					yerr = (probabilities_dist*(1-probabilities_dist)/n)**0.5
					ax[index_opt].errorbar(ind, probabilities_dist, yerr = yerr, fmt='s', color="darkviolet", alpha=0.5, linewidth=3., markersize=red_blue_yellow_cyan_marker_size)
					ax[index_opt].plot(ind, probabilities_dist, alpha=1.0, color="darkviolet", marker= "s", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "distance")
				if "correl" in type_of_models: 
					probabilities_correl = np.r_[precision_correl_pro[[index_100_correl_pro, index_200_correl_pro, index_300_correl_pro]], sensitivity_correl]
						#n = np.r_[np.array([True_positive_Rate_correl_pro[index_100_correl_pro], True_positive_Rate_correl_pro[index_200_correl_pro], True_positive_Rate_correl_pro[index_300_correl_pro]])*number_of_interacting_enhancers_, number_of_interacting_enhancers_]
					n = np.r_[length_of_positives_pro*True_positive_Rate_correl_pro[[index_100_correl_pro, index_200_correl_pro, index_300_correl_pro]] + False_positive_Rate_correl_pro[[index_100_correl_pro, index_200_correl_pro, index_300_correl_pro]]*length_of_negatives_pro, number_of_interacting_enhancers_]
					yerr = (probabilities_correl*(1-probabilities_correl)/n)**0.5
					ax[index_opt].errorbar(ind, probabilities_correl, yerr= yerr, fmt='o', color="red", alpha=0.5, linewidth=3., markersize=red_blue_yellow_cyan_marker_size)
					ax[index_opt].plot(ind, probabilities_correl, alpha=1.0, color="red", marker= "o", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "data")

					#-------------------------MOG
				if "MOG_correl_dist" in type_of_models: 
					probabilities_dist_correl_MOG = np.r_[precision_dist_correl_pro_MOG[[index_100_dist_correl_pro_MOG, index_200_dist_correl_pro_MOG, index_300_dist_correl_pro_MOG]], sensitivity_dist_correl_MOG]
					#n = np.array([True_positive_Rate_dist_correl_pro_MOG[index_100_dist_correl_pro_MOG], True_positive_Rate_dist_correl_pro_MOG[index_200_dist_correl_pro_MOG], True_positive_Rate_dist_correl_pro_MOG[index_300_dist_correl_pro_MOG]])*number_of_interacting_enhancers_
					n = np.r_[length_of_positives_pro*True_positive_Rate_dist_correl_pro_MOG[[index_100_dist_correl_pro_MOG, index_200_dist_correl_pro_MOG, index_300_dist_correl_pro_MOG]] + False_positive_Rate_dist_correl_pro_MOG[[index_100_dist_correl_pro_MOG, index_200_dist_correl_pro_MOG, index_300_dist_correl_pro_MOG]]*length_of_negatives_pro, number_of_interacting_enhancers_]
					yerr = (probabilities_dist_correl_MOG*(1-probabilities_dist_correl_MOG)/n)**0.5
					ax[index_opt].errorbar(ind, probabilities_dist_correl_MOG, yerr= yerr, fmt='*', color="cyan", alpha=0.5, linewidth=3., markersize=red_blue_yellow_cyan_marker_size)
					ax[index_opt].plot(ind, probabilities_dist_correl_MOG, alpha=1.0, color="cyan", marker= "*", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "LVA data+prior")

				if "MOG_dist" in type_of_models:
					probabilities_dist_MOG = np.r_[precision_dist_pro_MOG[[index_100_dist_pro_MOG, index_200_dist_pro_MOG, index_300_dist_pro_MOG]], sensitivity_dist_MOG]
					#n = np.array([True_positive_Rate_dist_pro_MOG[index_100_dist_pro_MOG], True_positive_Rate_dist_pro_MOG[index_200_dist_pro_MOG], True_positive_Rate_dist_pro_MOG[index_300_dist_pro_MOG]])*number_of_interacting_enhancers_
					n = np.r_[length_of_positives_pro*True_positive_Rate_dist_pro_MOG[[index_100_dist_pro_MOG, index_200_dist_pro_MOG, index_300_dist_pro_MOG]] + False_positive_Rate_dist_pro_MOG[[index_100_dist_pro_MOG, index_200_dist_pro_MOG, index_300_dist_pro_MOG]]*length_of_negatives_pro, number_of_interacting_enhancers_]
					yerr = (probabilities_dist_MOG*(1-probabilities_dist_MOG)/n)**0.5
					ax[index_opt].errorbar(ind, probabilities_dist_MOG, yerr= yerr, fmt="v", color="green", alpha=0.5, linewidth=3., markersize=red_blue_yellow_cyan_marker_size)
					ax[index_opt].plot(ind, probabilities_dist_MOG, alpha=1.0, color="green", marker= "v", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "LVA prior")






				handles, labels = ax[0].get_legend_handles_labels()
				ax[0].legend(handles, labels, fontsize = legend_box_names_font_size, numpoints=1, handletextpad = 0.0, borderpad=0.2, labelspacing=0.2)
				


				#[index_opt].plot(ind, np.r_[precision_correl_pro[[index_100_correl_pro, index_200_correl_pro, index_300_correl_pro]], sensitivity_match_MAP["correl"][",".join(np.array(option_, str))]], alpha=1.0, color="red", marker= "o", linewidth=0.0)



			pdf.savefig()
	if config_variables.FDR_mode: output.close(); 
	output_2.close()
	pdf.close()	
	#plt.show()
	plt.close("all")



