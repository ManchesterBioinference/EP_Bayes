def executor(selection_option):

	import config_variables
	import itertools
	import numpy as np
	from  prepare_interactions_clean import un_string

	results_folder = config_variables.results_folder

	name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes = results_folder + "file_with_FDRs_{0}_{1}_smo_{2}_{3}".format("chr1", "chr2", config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_TSS"
	else:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_GENE"

	FDR_thresholds = np.loadtxt(name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes, dtype = str, delimiter = "\t")

	classifiers_clean = config_variables.classifiers_clean

	#cluster_specific = False
	normalised = False
	cluster_specific = False

	chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
	chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh

	dict_chrom_enh_survived = config_variables.dict_chrom_enh_survived
	dict_chrom_distant = config_variables.dict_chrom_distant

	negative_interactions = config_variables.negative_interactions

	mode = config_variables.mode
	chroms_to_infer = config_variables.chroms_to_infer

	#gene_names = np.loadtxt("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered", dtype = str, usecols = (3,))	

	stuff = [1, 2, 3, 4]
	combinations = []

	for L in range(0, len(stuff)+1):
		for subset in itertools.combinations(stuff, L):
			if len(subset): combinations += [list(subset)]

	selected_combinations = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()

	dict_option = dict(zip(range(len(config_variables.datasets_names)), config_variables.datasets_names))

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

	def interactions_above_threshold(posterior_t, posterior_f, threshold_up, threshold_low, chrom, label = 0, domain = False):
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

		#probability_matrix += probability_matrix.T

		#--------------------------------------------------------------------------the part to look at for Paolo !
			
		if normalised:
			norm_factors = np.sum(probability_matrix, axis = 0)
			probability_matrix = probability_matrix/norm_factors
			probability_matrix[np.isnan(probability_matrix)] = 0.

		#--------------------------------------------------
		#how many true links for threshold estimated from test data at FDR
		
		interacting_mask = np.zeros_like(probability_matrix).astype(bool)
		interacting_mask[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = True

		number_of_true_interactions_at_prob_threshold_estimated_from_test_data_at_FDR_chrom = clustering_matrix_mask * (probability_matrix >= threshold_low) * (probability_matrix <= threshold_up) #* interacting_mask

		return number_of_true_interactions_at_prob_threshold_estimated_from_test_data_at_FDR_chrom.sum()

	def get_table_of_FDR_predictions(option_, thresholds_test_est_FDR):
		#add loop which would take FDR and mapped probability_thresholds_estimated_form_test_data_to_produce_files_specific_to_it

		posterior_correl_dist_true, posterior_correl_dist_false = classifiers_clean.posterior_producer([0], option_, total_posterior = False)

		number_of_true_interactions_at_prob_threshold_estimated_from_test_data_at_FDR = 0

		for chrom in chroms_to_infer:
			number_of_true_interactions_at_prob_threshold_estimated_from_test_data_at_FDR += interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, 1., thresholds_test_est_FDR, chrom, label = 0, domain = False)

		return number_of_true_interactions_at_prob_threshold_estimated_from_test_data_at_FDR

	input_head = ("#" + "\t".join(["Data_set", "FDR", "Precision", "True_links_above_FDR", "Threshold"]) + "\n")

	option_data_dist = selected_combinations[selection_option]
	comb = ",".join([dict_option[el] for el in option_data_dist])
	number_of_positives_FDR_data_dist, number_of_positives_FDR_dist = [], []	

	for FDR in config_variables.FDR[::-1]:

		#number_of_positives_FDR_data_dist, number_of_positives_FDR_dist = [], []

		thresholds_test_est_FDR_dist_data, thresholds_test_est_FDR_dist, thresholds_test_est_FDR_data = FDR_thresholds[(FDR_thresholds[:,0] == comb) * (FDR_thresholds[:,2].astype(float) == FDR), -1]

		#number_of_positives_FDR_data_dist += [get_table_of_FDR_predictions(option_data_dist, float(thresholds_test_est_FDR_dist_data))]


		if thresholds_test_est_FDR_dist_data <> 'N\A':
 			number_of_positives_FDR_data_dist += [get_table_of_FDR_predictions(option_data_dist, float(thresholds_test_est_FDR_dist_data))]
		else:
			number_of_positives_FDR_data_dist += [0]


		if thresholds_test_est_FDR_dist <> 'N\A': 
			number_of_positives_FDR_dist += [get_table_of_FDR_predictions([], float(thresholds_test_est_FDR_dist))]
		else: 
			number_of_positives_FDR_dist += [0] 

	to_save = np.column_stack((config_variables.FDR[::-1], number_of_positives_FDR_data_dist, number_of_positives_FDR_dist, np.array(number_of_positives_FDR_data_dist)/np.array(number_of_positives_FDR_dist).astype(float)))

	name_of_output_file_distance_vs_data_and_distance = results_folder + "table_distance_vs_{0}_distance_vs_FDR_smo_{1}_{2}".format(comb, config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_file_distance_vs_data_and_distance += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_TSS"
	else:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_GENE"

	np.savetxt(name_of_output_file_distance_vs_data_and_distance, to_save, delimiter = "\t", header = "\t".join(["FDR", "data/distance", "distance", "ratio"]))

