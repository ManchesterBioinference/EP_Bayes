def executor(PR_CURVES = "SELECTIVE", mode_of_code = "EVEN", mode_of_features_and_interactions = "FEATURES_AND_INTERACTIONS_TOGETHER", GENE_OR_PROMOTER_MODE = "GENE_MODE", redo_raw_CHIA_PET_interactions = True, rising_of_falling_POL2_or_ER = ["ER", "Pol2"][1], plot_TF_enrichments_in_cluster = False, upstream = 300, downstream = 0, upstream_t_s = 300, downstream_t_s = 0, do_clustering = False, re_do_clustering = False, cluster_figure_selection = None, DB_version = False):

	import numpy as np
	import re
	from sys import argv
	import matplotlib.pyplot as plt
	import itertools
	import bisect as bis
	import random as random
	import time
	import kern_density_est
	import smooth_priors
	import smooth_priors_non_domain
	import smooth_correl
	import smooth_priors_domain
	from matplotlib.backends.backend_pdf import PdfPages

	#np.seterr(all=None, divide='raise', over=None, under=None, invalid=None)

	if mode_of_features_and_interactions == "FEATURES_AND_INTERACTIONS_TOGETHER":
		disentagled_features_validation = False # it just mean that it's either gene or TSS mode
		upstream_t_s = upstream
		downstream_t_s = downstream

	elif mode_of_features_and_interactions == "FEATURES_AND_INTERACTIONS_SEPERATE":
		if GENE_OR_PROMOTER_MODE == "GENE_MODE": disentagled_features_validation = False
		if GENE_OR_PROMOTER_MODE == "TSS_MODE": disentagled_features_validation = True
			


	filter_value, filter_enh, count_f_p, count_f_e, ER_pro_filtered_, path = '-1.', '-1.', 30, 30, 'True', 1 # change to 2 if you want path 2 interactions
	alternative_classificator = True
	alternative_classificator_outside_enhancers = False#True # that option is for domain. Generator
	domain = False
	domain_like_chromosome_correction = False
	interacting_negatives = False
	log_distances = True
	plot_atr, plot_atr_kernel = False, True
	use_smooth_prior_for_estimation = True
	likelihood_cross_validation = True
	distant_enh_only = True # matters for enhancer-enhancer interactions and MAPS for all enhancers not only the interacting ones
	filter_values = np.array([-1., -0.6, -0.2])
	filter_value = filter_values[0]
	number_of_bins = 4000, 4000 # to implement-easy

	FDR = np.array([0.15, 0.2 ,  0.25,  0.3 ,  0.35,  0.4])

	#scripts:

	if mode_of_code == "FULL":
		chroms_in_prior = np.arange(0,23,1)#+1#np.arange(0,13,1)#np.arange(0,13,1)
		chroms_to_infer = np.arange(0,23,1)#np.arange(0,23,2)#np.arange(0,13,1)#np.arange(0,23,2)#np.arange(0,13,1)
		FDR_mode = False
		interacting_enhancers_only = False
		TOP_PR_interaction_plotter_clean_chrom_to_plot = chroms_to_infer[1]	
		option_for_predictive_FULL_mode = 2
		genes_predicted_with_FDR_for_GRO_seq_validation = 0.25
		TOP_PR_interaction_plotter_FDR_thresholds_to_plot = FDR[:3]

	elif mode_of_code == "ODD":
		chroms_in_prior = np.arange(0,23,2)#np.arange(0,13,1)#np.arange(0,13,1)
		chroms_to_infer = np.arange(0,23,2)#np.arange(0,23,2)#np.arange(0,13,1)#np.arange(0,23,2)#np.arange(0,13,1)
		FDR_mode = True # that fuction apply for odd-odd and odd-even only. 
		interacting_enhancers_only = True

	elif mode_of_code == "EVEN":
		chroms_in_prior = np.arange(0,23,2)#np.arange(0,13,1)#np.arange(0,13,1)
		chroms_to_infer = np.arange(0,22,2)+1#np.arange(0,23,2)#np.arange(0,13,1)#np.arange(0,23,2)#np.arange(0,13,1)
		FDR_mode = True
		interacting_enhancers_only = True


	

	mode = ["promoter_enhancer_interactions", "enhancer_enhancer_interactions"][0]
	one_sided_or_two_sided = ["single_sided", "double_sided"][1]
	TSS_or_intra_genic_for_domain_filter = ["Intra_genic", "TSS_only"][0]
	generator_mode = ["filter_independent_generator", "filter_correl_dependent_generator", "filter_dependent_generator"][1]

	promoter_overlaps_enhancer_file = "intersect_with_full_genes_l_{0}_r_{1}".format(upstream, downstream)
	name_of_promoter_file_for_overlap = "Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed"
	name_of_enhancer_file_for_overlap = "common_region_peaks_extended_less_time_points_corrected_0_indexed"#"common_region_peaks_extended_less_time_points_sorted"
	name_of_time_series_promoter_file_for_TSS_start = 'Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered'
	name_of_overlap_file_pro = 'ER_promoters_{0}_{1}'.format(upstream, downstream)
	name_of_overlap_file_enh = 'ER_peaks_overlapping_promoters_{0}_{1}'.format(upstream, downstream)

	# you can now make every feature to behave differentely. However how to count a signal would depend on where TSS is.

	#----------------------------------------------------------
	#MoG

	MoG_classificator, Run_MoG_classificator = False, False
	kappa_0, mu_0, alpha_0, Beta_0 = 1.0, 0.0, 5.0, 2.0
	number_of_samples = 100001#30

	#----------------------------------------------------------

	chrom_names = np.array(map(lambda x: "chr{0}".format(x), np.r_[np.arange(1, 23).astype(dtype='S2'), ['X'], ['Y']]))
	chroms_to_infer = chrom_names[chroms_to_infer]
	chroms_in_prior = chrom_names[chroms_in_prior]
	option = [0,1,2,3,4]
	filt_option = option
	time_points = 8



	datasets_names = np.array(['PolII_2012-03', 'PolII', 'H2AZ', 'ER', 'H3K4me3'])#, '2012-03_RNA', 'RNA'])
	dataset_names_option = datasets_names[option]
	dict_option = dict(zip(range(len(datasets_names)), datasets_names))

	link_data_set_name_to_file_name = {}

	name_of_time_series_file = {}
	name_of_time_series_file["enhancers"] = name_of_enhancer_file_for_overlap + "_unfiltered_count"

	if upstream_t_s <> 0: name_of_time_series_file["promoters"] = "Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_{0}_unfiltered_count".format(upstream_t_s)
	else: name_of_time_series_file["promoters"] = "Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_unfiltered_count"

	full_list_enhancers = np.array([name_of_time_series_file["enhancers"] + "_{0}".format(name_of_TF) for name_of_TF in datasets_names])
	full_list_promoters = np.array([name_of_time_series_file["promoters"] + "_{0}".format(name_of_TF) for name_of_TF in datasets_names])

	link_data_set_name_to_file_name["enhancers"] = dict(zip(datasets_names, full_list_enhancers))
	link_data_set_name_to_file_name["promoters"] = dict(zip(datasets_names, full_list_promoters))



	#------------------------------------------------------------------------------------------------------------
	import time_series_prepare_filter as initiate_time_series
	initiate_time_series.datasets_names = datasets_names
	initiate_time_series.time_points = time_points
	dataset_time_series_dict = initiate_time_series.time_series_prepare(full_list_promoters[option], full_list_enhancers[option])


	#------------------------------------------------------------------------------------------------------------------

	classificator_elements = {}

	for filter_value_ in filter_values:
		classificator_elements[filter_value_] = {}
		for mode_ in ["promoter_enhancer_interactions", "enhancer_enhancer_interactions"]:
			classificator_elements[filter_value_][mode_] = {}
			for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
				classificator_elements[filter_value_][mode_][classification_of_interactions] = {}
				for attribute_of_interaction in ["distance", "correlation"]:
					classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction] = {}
					for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:

						classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative] = {}
						classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative]["prior_bins"] = np.array([])
						classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative]["prior_frequencies"] = np.array([])

						if attribute_of_interaction == "correlation":

							for data_set_name in dataset_names_option:
								classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][data_set_name] = {}
								classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][data_set_name]["prior_bins"] = np.array([])
								classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][data_set_name]["prior_frequencies"] = np.array([])
								classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][data_set_name]["posterior_component_values"] = {}
								for chrom_ in chroms_to_infer:
									classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][data_set_name]["posterior_component_values"][chrom_] = np.array([])

						else:
							for chrom_ in chroms_to_infer:
								classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative]["posterior_component_values"] = {}
								classificator_elements[filter_value_][mode_][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative]["posterior_component_values"][chrom_] = np.array([])




	import config_variables
	config_variables.np = np
	config_variables.link_data_set_name_to_file_name = link_data_set_name_to_file_name
	config_variables.chroms_in_prior = chroms_in_prior
	config_variables.mode = mode
	config_variables.dataset_names_option = dataset_names_option
	config_variables.count_f_p = count_f_p
	config_variables.count_f_e = count_f_e
	config_variables.filter_enh = filter_enh
	config_variables.domain = domain
	config_variables.dataset_time_series_dict = dataset_time_series_dict
	config_variables.re = re
	config_variables.path = path
	config_variables.upstream = upstream
	config_variables.interacting_negatives = interacting_negatives
	config_variables.interacting_enhancers_only = interacting_enhancers_only
	config_variables.chroms_to_infer = chroms_to_infer
	config_variables.filter_value = filter_value
	config_variables.filter_values = filter_values
	config_variables.datasets_names = datasets_names
	config_variables.full_list_promoters = full_list_promoters
	config_variables.option = option
	config_variables.time_points = time_points
	config_variables.distant_enh_only = distant_enh_only
	config_variables.full_list_enhancers = full_list_enhancers
	config_variables.ER_pro_filtered_ = ER_pro_filtered_
	config_variables.TSS_or_intra_genic_for_domain_filter = TSS_or_intra_genic_for_domain_filter
	config_variables.one_sided_or_two_sided = one_sided_or_two_sided
	config_variables.chrom_names = chrom_names
	config_variables.promoter_overlaps_enhancer_file = promoter_overlaps_enhancer_file
	config_variables.name_of_time_series_promoter_file_for_TSS_start = name_of_time_series_promoter_file_for_TSS_start
	config_variables.upstream = upstream
	config_variables.downstream = downstream
	config_variables.upstream_t_s = upstream_t_s
	config_variables.name_of_promoter_file_for_overlap = name_of_promoter_file_for_overlap
	config_variables.name_of_enhancer_file_for_overlap = name_of_enhancer_file_for_overlap
	config_variables.name_of_overlap_file_pro = name_of_overlap_file_pro
	config_variables.name_of_overlap_file_enh = name_of_overlap_file_enh
	config_variables.filt_option = filt_option
	config_variables.log_distances = log_distances
	config_variables.domain_like_chromosome_correction = domain_like_chromosome_correction
	config_variables.alternative_classificator = alternative_classificator
	config_variables.likelihood_cross_validation = likelihood_cross_validation
	config_variables.alternative_classificator_outside_enhancers = alternative_classificator_outside_enhancers
	config_variables.dict_option = dict_option
	config_variables.kappa_0, config_variables.mu_0, config_variables.alpha_0, config_variables.Beta_0 = kappa_0, mu_0, alpha_0, Beta_0
	config_variables.MoG_classificator, config_variables.Run_MoG_classificator = MoG_classificator, Run_MoG_classificator
	config_variables.number_of_samples = number_of_samples
	config_variables.use_smooth_prior_for_estimation = use_smooth_prior_for_estimation
	config_variables.FDR = FDR
	config_variables.FDR_mode = FDR_mode
	config_variables.number_of_bins = number_of_bins
	config_variables.disentagled_features_validation = disentagled_features_validation
	config_variables.mode_of_code = mode_of_code

	#-----------------------------------------------

	#prepares variables and calculates model for a filter_value





	#run twice to get TSS and GENE mode.
	if redo_raw_CHIA_PET_interactions: import interaction_finder_wrapper


	import filters_clean

	#if not(domain) or alternative_classificator:
	#	dict_chrom_pro_survived, dict_chrom_pro_not_survived, filtered_promoters, Pol_2_correl_filtered_promoters = filters_clean.features_filtered(filter_value, count_f_p, full_list_promoters, filt_option, name_of_overlap_file_pro, add_overl = False)
	#	dict_chrom_enh_survived, dict_chrom_enh_not_survived, filtered_enhancers, Pol_2_correl_filtered_enhancers = filters_clean.features_filtered(filter_enh, count_f_e, full_list_enhancers, filt_option, name_of_overlap_file_enh, add_overl = False)
	#else:
	#	if mode == "promoter_enhancer_interactions":
	#		dict_chrom_pro_survived, dict_chrom_pro_not_survived, filtered_promoters, Pol_2_correl_filtered_promoters = filters_clean.features_filtered(filter_value, count_f_p, full_list_promoters, filt_option, name_of_overlap_file_pro, add_overl = False, remove_single_domain_elements = True)

	#		filter_, count_f, list_of_datasets, options, name_of_overlap_file, add_overl, remove_single_domain_elements = filter_value, count_f_p, full_list_promoters, filt_option, name_of_overlap_file_pro, False, True

	#		dict_chrom_enh_survived, dict_chrom_enh_not_survived, filtered_enhancers, Pol_2_correl_filtered_enhancers = filters_clean.features_filtered(filter_enh, count_f_e, full_list_enhancers, filt_option, name_of_overlap_file_enh, add_overl = False)

	#	else:
	#		dict_chrom_pro_survived, dict_chrom_pro_not_survived, filtered_promoters, Pol_2_correl_filtered_promoters = filters_clean.features_filtered(filter_value, count_f_p, full_list_promoters, filt_option, name_of_overlap_file_pro, add_overl = False)
	#		dict_chrom_enh_survived, dict_chrom_enh_not_survived, filtered_enhancers, Pol_2_correl_filtered_enhancers = filters_clean.features_filtered(filter_enh, count_f_e, full_list_enhancers, filt_option, name_of_overlap_file_enh, add_overl = False, remove_single_domain_elements = True)

	dict_chrom_pro_survived, dict_chrom_pro_not_survived, filtered_promoters, Pol_2_correl_filtered_promoters = filters_clean.features_filtered(filter_value, count_f_p, full_list_promoters, filt_option, name_of_overlap_file_pro, add_overl = False)
	dict_chrom_enh_survived, dict_chrom_enh_not_survived, filtered_enhancers, Pol_2_correl_filtered_enhancers = filters_clean.features_filtered(filter_enh, count_f_e, full_list_enhancers, filt_option, name_of_overlap_file_enh, add_overl = False)

	config_variables.Pol_2_correl_filtered_promoters = Pol_2_correl_filtered_promoters
	config_variables.Pol_2_correl_filtered_enhancers = Pol_2_correl_filtered_enhancers

	config_variables.dict_chrom_distant, config_variables.dict_chrom_proximal, config_variables.proximal_enhancers_mask = filters_clean.distant_enh_only_filter(name_of_overlap_file_enh)




	if do_clustering: 

		correl_value_filter = False
		distant_enh_only_log = False

		if cluster_figure_selection == "cluster_ER_enhancer": cluster_mode_setting = ["promoters", "enhancers"][1]; datasets_to_concat = datasets_names[[3]]; filter_each_dataset = 100; correl_value_filter = False; distant_enh_only_log = True;

		elif cluster_figure_selection == "cluster_Pol2s_enhancer": cluster_mode_setting = ["promoters", "enhancers"][1]; datasets_to_concat = datasets_names[[0, 1]]; filter_each_dataset = 30; correl_value_filter = 0.2; distant_enh_only_log = True;

		elif cluster_figure_selection == "cluster_Pol2s_promoter": cluster_mode_setting = ["promoters", "enhancers"][0]; datasets_to_concat = datasets_names[[0,1]]; filter_each_dataset = 30; correl_value_filter = 0.2;

		elif cluster_figure_selection == "cluster_ER_promoter": cluster_mode_setting = ["promoters", "enhancers"][0]; datasets_to_concat = datasets_names[[3]]; filter_each_dataset = 100;

		elif cluster_figure_selection == "cluster_Pol2s_ER_enhancer": cluster_mode_setting = ["promoters", "enhancers"][1]; datasets_to_concat = datasets_names[[1, 3]]; filter_each_dataset = 200;

		elif cluster_figure_selection == "cluster_Pol2s_ER_enhancer_test": cluster_mode_setting = ["promoters", "enhancers"][1]; datasets_to_concat = datasets_names[[1, 3]]; filter_each_dataset = 300;

		config_variables.dataset_time_series_dict_mean_std = initiate_time_series.time_series_prepare_mean_std(full_list_promoters[option], full_list_enhancers[option])		
		config_variables.name_of_time_series_file = name_of_time_series_file
		import AP_clustering	
		config_variables.name_of_overlap_file_dict = dict(zip(["promoters", "enhancers"], [name_of_overlap_file_pro, name_of_overlap_file_enh]))
		merged_time_series_to_cluster = AP_clustering.concatenator(
		cluster_mode = cluster_mode_setting, 
		merge_time_series_option = datasets_to_concat, 
		count_filter_each_data_set = filter_each_dataset, 
		pol2_rep_correl_filt = correl_value_filter, 
		distant_enh_only = distant_enh_only_log)

		if re_do_clustering:

			config_variables.merged_time_series_to_cluster = merged_time_series_to_cluster	
			AP_clustering.AP_clustering(merged_time_series_to_cluster, number_of_clusters = 40)
			config_variables.labels = np.loadtxt(merged_time_series_to_cluster + "_labels", dtype = str)

		import os as os

		path_to_R = os.getcwd() + "/R_scripts/"
		
		if cluster_figure_selection == "cluster_ER_enhancer": os.system("Rscript " + path_to_R + "ER_enhancer.R") 

		elif cluster_figure_selection == "cluster_Pol2s_enhancer": os.system("Rscript " + path_to_R + "PolIIs_enhancer.R")

		elif cluster_figure_selection == "cluster_Pol2s_promoter": os.system("Rscript " + path_to_R + "PolIIs_promoter.R")

		elif cluster_figure_selection == "cluster_ER_promoter": os.system("Rscript " + path_to_R + "ER_promoter.R")

		elif cluster_figure_selection == "cluster_Pol2s_ER_enhancer": os.system("Rscript " + path_to_R + "Pol2_ER.R")

		import sys
		sys.exit("exiting code")

	if plot_TF_enrichments_in_cluster:

		merged_time_series_to_cluster = "common_region_peaks_extended_less_time_points_corrected_0_indexed_unfiltered_count_concat_PolII_ER_200"
		import overlapper_hg19_clean
		overlapper_hg19_clean.executor(merged_time_series_to_cluster, diff_bind_version = DB_version, mode_atr = rising_of_falling_POL2_or_ER) # mode attribute specifies whether it should use ER mean or Pol2 mean of a cluster to assess raising or falling tendencies.
		import sys
		sys.exit("exiting code")

	import generator_executor
	f_name = generator_executor.interactions_producer_filter(generator_mode, domain, 2, TSS_or_intra_genic_for_domain_filter, "GENE_MODE") #in order to get path 2 interactions change to 3

	config_variables.dict_chrom_pro_survived = dict_chrom_pro_survived
	config_variables.dict_chrom_pro_not_survived = dict_chrom_pro_not_survived
	config_variables.f_name = f_name
	config_variables.filtered_promoters = filtered_promoters
	config_variables.filtered_enhancers = filtered_enhancers
	config_variables.dict_chrom_enh_survived = dict_chrom_enh_survived
	config_variables.dict_chrom_enh_not_survived = dict_chrom_enh_not_survived

	import prepare_interactions_clean

	alternative_classificator_outside_enhancers = True # had something to do with enhancers outside domains
	if alternative_classificator_outside_enhancers:
		f_name_2 = generator_executor.interactions_producer_filter(generator_mode, True, 2, TSS_or_intra_genic_for_domain_filter, "GENE_MODE")
		chr_interactions_dict_pro_enh, chr_interactions_dict_enh_enh, dict_total_enh, dict_total_pro = prepare_interactions_clean.filter_true_interactions_of_promoters_and_enhancers_which_didnt_survive_filtering(f_name_2)
		from  prepare_interactions_clean import un_string
		chrom_interacting_enhancers_pro = {}	
		for chrom__ in chrom_names: chrom_interacting_enhancers_pro[chrom__] = np.unique(un_string(chr_interactions_dict_pro_enh[chrom__])[:,1])
		config_variables.chrom_interacting_enhancers_pro = chrom_interacting_enhancers_pro

	chr_interactions_dict_pro_enh, chr_interactions_dict_enh_enh, dict_total_enh, dict_total_pro = prepare_interactions_clean.filter_true_interactions_of_promoters_and_enhancers_which_didnt_survive_filtering(f_name)


	if disentagled_features_validation: 
		f_name_TSS = generator_executor.interactions_producer_filter(generator_mode, domain, 2, TSS_or_intra_genic_for_domain_filter, "TSS_MODE")
		config_variables.chr_interactions_dict_pro_enh_TSS, config_variables.chr_interactions_dict_enh_enh_TSS, config_variables.dict_total_enh_TSS, config_variables.dict_total_pro_TSS = prepare_interactions_clean.filter_true_interactions_of_promoters_and_enhancers_which_didnt_survive_filtering(f_name_TSS)

	config_variables.dict_total_enh = dict_total_enh
	config_variables.dict_total_pro = dict_total_pro
	config_variables.chr_interactions_dict_pro_enh = chr_interactions_dict_pro_enh
	config_variables.chr_interactions_dict_enh_enh = chr_interactions_dict_enh_enh

	import chrom_specific_negative_interactions as negative_interactions

	config_variables.negative_interactions = negative_interactions

	import prior_producer
	import classificator_clean
	import prepare_upper_and_lower_bounds_for_priors as prior_bounds
	import prior_histograms_cl
	import allocator

	config_variables.alternative_classificator_outside_enhancers = False
	prior_elements = prior_producer.prior_producer()
	config_variables.alternative_classificator_outside_enhancers = False#True
	infered_elements = classificator_clean.infered_elements_filler()

	low_dist, up_dist = prior_bounds.prepare_upper_and_lower_bounds_for_priors(prior_elements, infered_elements)
	prior_elements = prior_histograms_cl.prior_bins_prob_and_plotter(prior_elements, low_dist, up_dist, use_smooth_prior_for_estimation, plot_atr, plot_atr_kernel)
	infered_elements = allocator.allocator(infered_elements, prior_elements)

	#for mode in modes:
	#	for filter_value in filter_values:


	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		for attribute_of_interaction in ["distance", "correlation"]:
			for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
				if attribute_of_interaction == "correlation":
					for data_set_name in dataset_names_option:
						for chrom_ in chroms_to_infer:
							update = infered_elements[mode][classification_of_interactions][attribute_of_interaction][data_set_name][probability_of_being_positive_or_negative][chrom_]
							classificator_elements[filter_value][mode][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][data_set_name]["posterior_component_values"][chrom_] = update

				elif attribute_of_interaction == "distance":
					for chrom_ in chroms_to_infer:
						update = infered_elements[mode][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative][chrom_]
						classificator_elements[filter_value][mode][classification_of_interactions][attribute_of_interaction][probability_of_being_positive_or_negative]["posterior_component_values"][chrom_] = update


	def func_star(args): return MOG.executor(*args)

	if config_variables.Run_MoG_classificator:
		from multiprocessing import Pool
		config_variables.probabilities_of_a_bin = prior_elements[mode]["positive_interactions"]["distance"]["prior_frequencies"]
		config_variables.adequate_histogram_bins = prior_elements[mode]["positive_interactions"]["distance"]["prior_bins"]

		import finite_MOG_object_orientated_1d_times_n_case_log_calc_prob_visited_float64_distance_low_distances_active_promoters_clean as MOG
		stuff = [1, 2, 3, 4]
		combinations = []
		for L in range(0, len(stuff)+1):
			for subset in itertools.combinations(stuff, L):
				if len(subset): combinations += [list(subset)]

		selected_combinations = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()
		p = Pool(5)


		#for option_correl__ in selected_combinations:
		#	for chrom_ in chroms_to_infer:
		#		MOG.executor(number_of_samples, option_correl__, chrom_)



		arguments = [(number_of_samples, option_correl__, chrom_) for chrom_ in chroms_to_infer for option_correl__ in selected_combinations] #[(number_of_samples, option_correl__, chrom_) for option_correl__ in selected_combinations]

		p.map_async(func_star, arguments)
	





	config_variables.classificator_elements = classificator_elements
	import classifiers_clean
	config_variables.classifiers_clean = classifiers_clean


	if mode_of_code == "ODD" or mode_of_code == "EVEN":
		#import PR_top
		#PR_top.execute()
		#import PR_top_MAP_dots
		import MAP_invoker
		MAP_probabilites_correl_dist, infered_elements_correl_dist, match_MAP_correl_dist, sensitivity_match_MAP = MAP_invoker.executor()
		import PR_top_MAP_dots_alternative_domain
		PR_top_MAP_dots_alternative_domain.execute(sensitivity_match_MAP, np.sum([len(match_MAP_correl_dist[chrom_]) for chrom_ in chroms_to_infer]), option_to_plot = PR_CURVES)



	if mode_of_code == "FULL":
		#import MAP_clustering_labels_clean
		#MAP_clustering_labels_clean.executor(MAP_probabilites_correl_dist, infered_elements_correl_dist)

		import TOP_FDR_PR_gene_list_clean
		TOP_FDR_PR_gene_list_clean.executor(selection_option = option_for_predictive_FULL_mode)
		import TOP_FDR_PR_table_clean
		TOP_FDR_PR_table_clean.executor(selection_option = option_for_predictive_FULL_mode)
	
		import script_python_analys_PR
		script_python_analys_PR.executor(selection_option = option_for_predictive_FULL_mode, FDR_level = genes_predicted_with_FDR_for_GRO_seq_validation)

		#import MAP_interaction_plotter_clean
		#MAP_interaction_plotter_clean.executor(MAP_probabilites_correl_dist, infered_elements_correl_dist, match_MAP_correl_dist)
		import TOP_PR_interaction_plotter_clean
		TOP_PR_interaction_plotter_clean.executor(selection_option = option_for_predictive_FULL_mode, chrom_to_plot = TOP_PR_interaction_plotter_clean_chrom_to_plot, FDR_thresholds_to_plot = TOP_PR_interaction_plotter_FDR_thresholds_to_plot)


