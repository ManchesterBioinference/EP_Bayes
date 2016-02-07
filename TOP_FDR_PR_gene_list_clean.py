def executor(selection_option):

	import config_variables
	import itertools
	import numpy as np
	from  prepare_interactions_clean import un_string
	classifiers_clean = config_variables.classifiers_clean

	cluster_specific = False
	normalised = False
	#filter_value = config_variables.filter_value
	#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
	#chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh

	#dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived
	dict_chrom_enh_survived = config_variables.dict_chrom_enh_survived
	dict_chrom_distant = config_variables.dict_chrom_distant

	#dataset_time_series_dict = config_variables.dataset_time_series_dict
	negative_interactions = config_variables.negative_interactions

	dict_option = config_variables.dict_option
	#link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
	#name_of_time_series_promoter_file_for_TSS_start = config_variables.name_of_time_series_promoter_file_for_TSS_start

	mode = config_variables.mode
	chroms_to_infer = config_variables.chroms_to_infer
	results_folder = config_variables.results_folder

	gene_names = np.loadtxt(config_variables.name_of_time_series_promoter_file_for_TSS_start, dtype = str, usecols = (3,))	

	stuff = [1, 2, 3, 4]
	combinations = []

	for L in range(0, len(stuff)+1):
		for subset in itertools.combinations(stuff, L):
			if len(subset): combinations += [list(subset)]

	selected_combinations = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()


	def positive_negative_interactions_for_MAP(chrom):
		

		indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

		if mode == "promoter_enhancer_interactions":

			false_inter_pro = negative_interactions.chrom_specific_negative_interactions(chrom, mode)

			i_s_f, j_s_f = false_inter_pro[:,0] + total_p, false_inter_pro[:,1] + total_e

		if mode == "enhancer_enhancer_interactions":

			false_inter_enh = negative_interactions.chrom_specific_negative_interactions(chrom, mode)
	
			i_s_f, j_s_f = false_inter_enh[:,0] + total_e, false_inter_enh[:,1] + total_e

		return i_s_f, j_s_f

	#------------------------------------------------------------------------------
	#for paolo
	def interactions_above_threshold(posterior_t, posterior_f, threshold_up, threshold_low, chrom, label, domain = False):
		enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)
		
		length_chr = len(indexes_p) + len(indexes_e)
		probability_matrix = np.zeros((length_chr, length_chr))
		posterior_t_chrom, posterior_f_chrom = posterior_t[chrom], posterior_f[chrom]

		i_s_f, j_s_f = positive_negative_interactions_for_MAP(chrom)

		if domain:
			if config_variables.TSS_or_intra_genic_for_domain_filter == "Intra_genic": coords_pro_domain = pro_coordinates[indexes_p]
			elif config_variables.TSS_or_intra_genic_for_domain_filter == "TSS_only": coords_pro_domain = np.column_stack((TSS_coordinates[indexes_p]-1, TSS_coordinates[indexes_p]+1))
			domain_matrix = interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'left', True)
			domain_matrix = domain_matrix + interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'right', True)
		else:
			domain_matrix = True


		if cluster_specific:
			
			if config_variables.distant_enh_only: 
				survived_mask = np.zeros(len(indexes_e))
				survived_mask[dict_chrom_enh_survived[chrom]] = True
				distance_mask = np.zeros(len(indexes_e))
				distance_mask[dict_chrom_distant[chrom]] = True
				survived = np.where(survived_mask*distance_mask)[0]
			else:
				survived = dict_chrom_enh_survived[chrom]
	
			clustering_matrix = np.zeros_like(probability_matrix).astype(int)		
			clustering_matrix[:len(indexes_p), len(indexes_p) + survived] = labels[survived]
			clustering_matrix += clustering_matrix.T
			clustering_matrix_mask = clustering_matrix == label

		else:
			#_matrix = np.zeros_like(probability_matrix).astype(int)	
			clustering_matrix_mask = True


		#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
		#true_inter_pro = un_string(chr_interactions_dict_pro_enh[chrom][:, :2]).astype(int)


		if config_variables.disentagled_features_validation: 
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
		else:
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

		#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
		#true_inter_pro = un_string(chr_interactions_dict_pro_enh[chrom][:, :2]).astype(int)
		true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

		i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

		#probability_matrix[:,:] = np.min([np.min(posterior_t), np.min(posterior_f)])*0.999
		probability_matrix[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = posterior_t_chrom
		probability_matrix[i_s_f - total_p, j_s_f + len(indexes_p) - total_e] = posterior_f_chrom

		probability_matrix += probability_matrix.T

			#--------------------------------------------------------------------------the part to look at for Paolo !
			
		if normalised:
			norm_factors = np.sum(probability_matrix, axis = 0)
			probability_matrix = probability_matrix/norm_factors
			probability_matrix[np.isnan(probability_matrix)] = 0.


		#-------------------------------------------------------------------- that's where cluster_specifcic can play part

		#---------------------- distal
		mask_for_counts = clustering_matrix_mask * (probability_matrix >= threshold_low) * (probability_matrix < threshold_up) #* interacting_mask

		number_of_enhancer_target_promoters_within_a_cluster_distal = mask_for_counts.sum(0)[:len(indexes_p)] # affected by cluster labels

		filtered_probability_matrix = np.zeros_like(probability_matrix) # events with probability one have got probability 0. ease the calculations below

		filtered_probability_matrix[mask_for_counts] = probability_matrix[mask_for_counts]

		#probability_of_enhancer_target_promoters_within_a_cluster = 1. - (1. - filtered_probability_matrix).prod(0)[:len(indexes_p)]  # affected by cluster labels
		probability_of_enhancer_target_promoters_within_a_cluster = -1*np.log(1. - filtered_probability_matrix).sum(0)[:len(indexes_p)]  # affected by cluster labels
		#---------------------- distal

		#---------------------- proximal

		promoter_overlaps_enhancer = np.loadtxt(config_variables.promoter_overlaps_enhancer_file, usecols = (0, 4, 8), dtype = str)

		promoter_overlaps_enhancer_chrom = promoter_overlaps_enhancer[promoter_overlaps_enhancer[:,0] == chrom]

		interacting_mask = np.zeros_like(probability_matrix).astype(bool)

		interacting_mask[promoter_overlaps_enhancer_chrom[:, 1].astype(int) - total_p, promoter_overlaps_enhancer_chrom[:, 2].astype(int) - total_e + len(indexes_p)] = True
		interacting_mask += interacting_mask.T

		interacting_mask = interacting_mask*clustering_matrix_mask
		
		number_of_enhancer_target_promoters_within_a_cluster_proximal = interacting_mask.sum(0)[:len(indexes_p)]

		#---------------------- distances
		TSS_coordinates = negative_interactions.extract_TSS_coordinates(config_variables.upstream)
		point_coordinates_promoter, point_coordinates_enhancer = TSS_coordinates[indexes_p], np.mean(enh_coordinates[indexes_e], axis = 1)
		distances_matrix = negative_interactions.calculate_distances(domain, point_coordinates_promoter, point_coordinates_enhancer)

		#distances_of_proximal = distances_matrix[[promoter_overlaps_enhancer_chrom[:, 1].astype(int) - total_p, promoter_overlaps_enhancer_chrom[:, 2].astype(int) - total_e + len(indexes_p)]]

		interacting_mask = np.zeros_like(probability_matrix).astype(bool)
		interacting_mask[len(indexes_p) + dict_chrom_enh_survived[chrom] - total_e] = True
		#interacting_mask[config_variables.dict_chrom_pro_survived[chrom] - total_p] = False
		#interacting_mask[len(indexes_p) + config_variables.dict_chrom_proximal[chrom] - total_e] = True


		interacting_mask += interacting_mask.T

		#def dist_filter_(distances_matrix, low_lim, up_lim):

		#	prox_distances_matrix = (distances_matrix <= up_lim) * (distances_matrix > low_lim) * interacting_mask  
		#	number_of_enhancer_target_promoters_within_a_cluster_proximal_within_distance_interval = prox_distances_matrix.sum(0)[:len(indexes_p)]
		#	return number_of_enhancer_target_promoters_within_a_cluster_proximal_within_distance_interval
		
		#lower_bound = proximal_grid[0]
		#does_proximal_enhancers_targetting_the_promoter_exist = np.zeros((len(indexes_p), len(proximal_grid) - 1))
		#for ind, dist_upper_bound in enumerate(proximal_grid[1:]): 
		#	number_of_times_a_promoter_is_targetted_by_a_proximal_enhancer_within_distances = dist_filter_(abs(distances_matrix), lower_bound, dist_upper_bound); 
		#	lower_bound = dist_upper_bound; 
		#	does_proximal_enhancers_targetting_the_promoter_exist[number_of_times_a_promoter_is_targetted_by_a_proximal_enhancer_within_distances.astype(bool), ind] = dist_upper_bound;

		#---------------------- distances

		def dist_filter_new(distances_matrix, low_lim, up_lim):
			prox_distances_matrix = (distances_matrix <= up_lim) * (distances_matrix > low_lim) * interacting_mask
			maximum_value = np.max(distances_matrix)
			distances_matrix_constr = np.zeros_like(prox_distances_matrix).astype(float) + maximum_value
			distances_matrix_constr[prox_distances_matrix] = distances_matrix[prox_distances_matrix]
			prox_distances_matrix_minimal_values = distances_matrix_constr.min(0)[:len(indexes_p)]
			prox_distances_matrix_minimal_values[prox_distances_matrix_minimal_values == maximum_value] = 0.

			return prox_distances_matrix_minimal_values

		#---------------------- proximal

		does_proximal_enhancers_targetting_the_promoter_exist = dist_filter_new(abs(distances_matrix), 0., 40000.)

		number_of_enhancer_target_promoters_within_a_cluster = number_of_enhancer_target_promoters_within_a_cluster_distal + number_of_enhancer_target_promoters_within_a_cluster_proximal

		non_zero_gene_count_mask = (number_of_enhancer_target_promoters_within_a_cluster <> 0) + does_proximal_enhancers_targetting_the_promoter_exist.astype(bool)#.sum(1).astype(bool)

		non_zero_counts = number_of_enhancer_target_promoters_within_a_cluster[non_zero_gene_count_mask]

		non_zero_counts_probabilities = probability_of_enhancer_target_promoters_within_a_cluster[non_zero_gene_count_mask]

		gene_names_non_zero_count = gene_names[np.where(non_zero_gene_count_mask)[0] + total_p]

		return  gene_names_non_zero_count, number_of_enhancer_target_promoters_within_a_cluster_distal[non_zero_gene_count_mask], number_of_enhancer_target_promoters_within_a_cluster_proximal[[non_zero_gene_count_mask]], non_zero_counts_probabilities, does_proximal_enhancers_targetting_the_promoter_exist[non_zero_gene_count_mask]

	def get_genes_function(labels, option_, thresholds_test_est_FDR):

		posterior_correl_dist_true, posterior_correl_dist_false = classifiers_clean.posterior_producer([0], option_, total_posterior = False)

		for ind_,chrom in enumerate(chroms_to_infer):

			for label in labels:
	
				gene_names_non_zero_count, number_of_enhancer_target_promoters_within_a_cluster_distal, number_of_enhancer_target_promoters_within_a_cluster_proximal, non_zero_counts_probabilities, does_proximal_enhancers_targetting_the_promoter_exist = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, 1., thresholds_test_est_FDR, chrom, label, domain = False)

				#q
				to_save_chr_label = np.column_stack((gene_names_non_zero_count, number_of_enhancer_target_promoters_within_a_cluster_distal, number_of_enhancer_target_promoters_within_a_cluster_proximal, number_of_enhancer_target_promoters_within_a_cluster_distal + number_of_enhancer_target_promoters_within_a_cluster_proximal, non_zero_counts_probabilities, does_proximal_enhancers_targetting_the_promoter_exist))

				if cluster_specific: 
					to_save_chr_label = np.column_stack((to_save_chr_label, [label]*len(to_save_chr_label))) 

			if ind_ == 0: to_save = np.array([]).reshape(0, to_save_chr_label.shape[1])

			to_save	= np.r_[to_save, to_save_chr_label]

		to_save_sorted = to_save[np.lexsort((to_save[:, 4].astype(float), to_save[:, 3].astype(float)))][::-1]

		return to_save_sorted

	input_head = ("#" + "\t".join(["Data_set", "FDR", "Precision", "True_links_above_FDR", "Threshold"]) + "\n")

	option_ = selected_combinations[selection_option]
	
	comb = ",".join([dict_option[el] for el in option_])


	name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes = results_folder + "file_with_FDRs_{0}_{1}_smo_{2}_{3}".format("chr1", "chr2", config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_TSS"
	else:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_GENE"

	FDR_thresholds = np.loadtxt(name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes, dtype = str, delimiter = "\t")

	#proximal_grid = np.arange(0, 20001, 200)

	for FDR in config_variables.FDR:
	
		thresholds_test_est_FDR_dist_data, thresholds_test_est_FDR_dist, thresholds_test_est_FDR_data = FDR_thresholds[(FDR_thresholds[:,0] == comb) * (FDR_thresholds[:,2].astype(float) == FDR), -1]

		if cluster_specific:

			labels = config_variables.labels
			to_save = get_genes_function(labels, option_, float(thresholds_test_est_FDR_dist_data))


			np.savetxt("clusters_genes_vs_counts_prob_distant_cluster_specific.txt", to_save, delimiter = "\t", fmt ="%s")	
	
		else: 
			proximal_header = ["min proximal value"] #["{0}kB".format(el) for el in proximal_grid]

			labels = [0]
			to_save = get_genes_function(labels, option_, float(thresholds_test_est_FDR_dist_data))

			name_of_output_file = results_folder + "clusters_genes_vs_counts_prob_distant_all_{0}_{1}_smo_{2}_proximal_version_PR_met".format(FDR, ",".join([comb]), config_variables.use_smooth_prior_for_estimation)

			name_of_output_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

			if config_variables.disentagled_features_validation:
				name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_TSS"
			else:
				name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_GENE"

			np.savetxt(name_of_output_file, to_save, delimiter = "\t", fmt ="%s", header = "\t".join(["gene_name", "distal_enhancer_targets", "proximal_enhancer_targets", "sum_of_two", "probabilities"] + proximal_header))


