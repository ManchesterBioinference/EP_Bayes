def executor(selection_option, chrom_to_plot, FDR_thresholds_to_plot):

	import config_variables
	import itertools
	import numpy as np
	from  prepare_interactions_clean import un_string
	classifiers_clean = config_variables.classifiers_clean
	

	normalised = False
	filter_value = config_variables.filter_value
	#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
	chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh

	dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived
	dict_chrom_enh_survived = config_variables.dict_chrom_enh_survived
	dict_chrom_distant = config_variables.dict_chrom_distant

	dataset_time_series_dict = config_variables.dataset_time_series_dict
	negative_interactions = config_variables.negative_interactions

	link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
	name_of_time_series_promoter_file_for_TSS_start = config_variables.name_of_time_series_promoter_file_for_TSS_start

	mode = config_variables.mode
	chroms_to_infer = config_variables.chroms_to_infer


	def calculate_single_ROC_best_True_sensitivity(probabilities_true, probabilities_false, threshold_1, threshold_2, threshold_3):

		_True_positives_of_threshold = []
		_False_positives_of_threshold = []

		sorted_prob_true = np.sort(probabilities_true)
		sorted_prob_false = np.sort(probabilities_false)

		sorted_thresholds = np.sort(np.unique(np.r_[probabilities_true, probabilities_false]))
		sorted_thresholds = np.unique(np.r_[np.min(sorted_thresholds)*0.999, sorted_thresholds, np.max(sorted_thresholds)*1.01]) # min should probably vanish from here take a look at PR_top_MAP_dots_alternative_domain it should be np.unique(np.r_[sorted_thresholds, np.max(sorted_thresholds)*1.01])

		len_prob_true = len(probabilities_true)
		len_prob_false = len(probabilities_false)

		print 'len prob: ', len_prob_true, len_prob_false

		_True_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_true, sorted_thresholds)[0][::-1])

		_False_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_false, sorted_thresholds)[0][::-1])

		Precision = np.array(_True_positives_of_threshold, dtype = float)/(np.array(_True_positives_of_threshold, dtype = float) + np.array(_False_positives_of_threshold, dtype = float))
	
		True_positive_Rate = np.array(_True_positives_of_threshold)/float(len_prob_true)
		
		False_positive_Rate = np.array(_False_positives_of_threshold)/float(len_prob_false)

		print 'number of thresholds', len(True_positive_Rate), len(False_positive_Rate)

		num_of_probabilities_above_threshold_1 = np.where(True_positive_Rate >= threshold_1)[0][0]
		num_of_probabilities_above_threshold_2 = np.where(True_positive_Rate >= threshold_2)[0][0]
		num_of_probabilities_above_threshold_3 = np.where(True_positive_Rate >= threshold_3)[0][0]

		sorted_thresholds_reverse = sorted_thresholds[::-1]

		return sorted_thresholds_reverse[num_of_probabilities_above_threshold_1 + 1], sorted_thresholds_reverse[num_of_probabilities_above_threshold_2 + 1], sorted_thresholds_reverse[num_of_probabilities_above_threshold_3 + 1]

	def extract_TSS_coordinates(upstream):

		data = np.loadtxt(name_of_time_series_promoter_file_for_TSS_start, dtype = str,  delimiter = '\t')	
		plus_strand = data[:, 4] == '+'
		TSS_coordinates = np.zeros(len(plus_strand), int)
		TSS_coordinates[plus_strand] = data[plus_strand, 1].astype(int) + upstream
		TSS_coordinates[np.invert(plus_strand)] = data[np.invert(plus_strand), 2].astype(int) + upstream

		return TSS_coordinates



	TSS_coordinates = extract_TSS_coordinates(config_variables.upstream)

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
	def interactions_above_threshold(posterior_t, posterior_f, threshold_up, threshold_low, chrom, domain = False):
		enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)
		
		length_chr = len(indexes_p) + len(indexes_e)
		interaction_matrix = np.zeros((length_chr, length_chr))
		posterior_t, posterior_f = posterior_t[chrom], posterior_f[chrom]

		i_s_f, j_s_f = positive_negative_interactions_for_MAP(chrom)

		if domain:
			import interacting_domain_clean as interacting_domain
			if config_variables.TSS_or_intra_genic_for_domain_filter == "Intra_genic": coords_pro_domain = pro_coordinates[indexes_p]
			elif config_variables.TSS_or_intra_genic_for_domain_filter == "TSS_only": coords_pro_domain = np.column_stack((TSS_coordinates[indexes_p]-1, TSS_coordinates[indexes_p]+1))
			domain_matrix = interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'left', True)
			domain_matrix = domain_matrix + interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'right', True)
		else:
			domain_matrix = True

		if mode == "promoter_enhancer_interactions":

			if config_variables.disentagled_features_validation: 
				chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
			else:
				chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

			#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
			#true_inter_pro = un_string(chr_interactions_dict_pro_enh[chrom][:, :2]).astype(int)
			true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

			i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

			#chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
			#true_inter_pro = un_string(chr_interactions_dict_pro_enh[chrom][:, :2]).astype(int)
			#i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

			interaction_matrix[:,:] = np.min([np.min(posterior_t), np.min(posterior_f)])*0.999
			interaction_matrix[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = posterior_t
			interaction_matrix[i_s_f - total_p, j_s_f + len(indexes_p) - total_e] = posterior_f

			#--------------------------------------------------------------------------the part to look at for Paolo !
			#np.save("interaction_matrix_float", interaction_matrix) #I'm saving the interaction matrix for you
	
			if normalised:
				norm_factors = np.sum(interaction_matrix, axis = 0)
				interaction_matrix = interaction_matrix/norm_factors
				interaction_matrix[np.isnan(interaction_matrix)] = 0.

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = True

			true_pro_enh_inter_filtered = np.where(interacting_mask * (interaction_matrix >= threshold_low) * (interaction_matrix < threshold_up))
			i_s_t_filt, j_s_t_filt = true_pro_enh_inter_filtered[0] + total_p, true_pro_enh_inter_filtered[1] - len(indexes_p) + total_e # that line tells you which of the positive (green) ChIA-PET-confirmed interactions lay within threshold_low, threshold_up interval

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_f - total_p, j_s_f + len(indexes_p) - total_e] = True

			false_pro_enh_inter_filtered = np.where(interacting_mask * (interaction_matrix >= threshold_low) * (interaction_matrix < threshold_up))
			i_s_f_filt, j_s_f_filt = false_pro_enh_inter_filtered[0] + total_p, false_pro_enh_inter_filtered[1] - len(indexes_p) + total_e # that line tells you which of the negative (gray) interactions lay within threshold_low, threshold_up interval

			#--------------------------------------------------------------------------the part to look at for Paolo !

			return	i_s_t_filt, j_s_t_filt, i_s_f_filt, j_s_f_filt # the function takes threshold_up, threshold_low as an argument and returns predicted interactions with probabilities within  [threshold_low, threshold_up) interval
			# i_s_t_filt, j_s_t_filt, i_s_f_filt, j_s_f_filt legend: i_s_t are promoters of ChIA-PET confirmed interactions, j_s_t are enhancers of the interactions, i_s_f, j_s_f are (promoters, enhancers) interactions which aren't ChIA-PET confirmed
			#for paolo

		if mode == "enhancer_enhancer_interactions":

			chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh
			true_inter_enh = un_string(chr_interactions_dict_enh_enh[chrom][:, :2]).astype(int)
			i_s_t, j_s_t = true_inter_enh[:,0], true_inter_enh[:,1]

			interaction_matrix[:,:] = np.min([np.min(posterior_t), np.min(posterior_f)])*0.999
			interaction_matrix[i_s_t + len(indexes_p) - total_e, j_s_t + len(indexes_p) - total_e] = posterior_t
			interaction_matrix[i_s_f + len(indexes_p) - total_e, j_s_f + len(indexes_p) - total_e] = posterior_f
			interaction_matrix[j_s_t + len(indexes_p) - total_e, i_s_t + len(indexes_p) - total_e] = posterior_t # transpose to create a full matrix
			interaction_matrix[j_s_f + len(indexes_p) - total_e, i_s_f + len(indexes_p) - total_e] = posterior_f # transpose to create a full matrix


			if normalised: 
				norm_factors = np.sum(interaction_matrix, axis = 0)
				interaction_matrix = interaction_matrix/norm_factors
				interaction_matrix[np.isnan(interaction_matrix)] = 0.

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_t + len(indexes_p) - total_e, j_s_t + len(indexes_p) - total_e] = True

			true_enh_enh_inter_filtered = np.where(interacting_mask * (interaction_matrix >= threshold_low) * (interaction_matrix < threshold_up))
			i_s_t_filt, j_s_t_filt = true_enh_enh_inter_filtered[0] - len(indexes_p) + total_e, true_enh_enh_inter_filtered[1] - len(indexes_p) + total_e

			interacting_mask = np.zeros_like(interaction_matrix).astype(bool)
			interacting_mask[i_s_f + len(indexes_p) - total_e, j_s_f + len(indexes_p) - total_e] = True

			false_enh_enh_inter_filtered = np.where(interacting_mask * (interaction_matrix >= threshold_low) * (interaction_matrix < threshold_up))
			i_s_f_filt, j_s_f_filt = false_enh_enh_inter_filtered[0] - len(indexes_p) + total_e, false_enh_enh_inter_filtered[1] - len(indexes_p) + total_e

			return i_s_t_filt, j_s_t_filt, i_s_f_filt, j_s_f_filt

	stuff = [1, 2, 3, 4]
	combinations = []

	for L in range(0, len(stuff)+1):
		for subset in itertools.combinations(stuff, L):
			if len(subset): combinations += [list(subset)]

	selected_combinations = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()

	option_ = selected_combinations[selection_option]
	
	#comb = ",".join([dict_option_[el] for el in option_])

	chrom = chrom_to_plot

	posterior_correl_dist_true, posterior_correl_dist_false = classifiers_clean.posterior_producer([0], option_, total_posterior = False)

	posterior_correl_dist_true_total = list(itertools.chain.from_iterable([posterior_correl_dist_true[chrom_] for chrom_ in chroms_to_infer]))
	posterior_correl_dist_false_total = list(itertools.chain.from_iterable([posterior_correl_dist_false[chrom_] for chrom_ in chroms_to_infer]))

	#threshold_1, threshold_2, threshold_3 = calculate_single_ROC_best_True_sensitivity(posterior_correl_dist_true_total, posterior_correl_dist_false_total, 0.1, 0.2, 0.3)
	#print "threshold_TPR_10%: ", threshold_1, "threshold_TPR_20%: ", threshold_2, "threshold_TPR_30%: ", threshold_3

	comb = ",".join([config_variables.dict_option[el] for el in option_])

	name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes = "file_with_FDRs_{0}_{1}_smo_{2}_{3}".format("chr1", "chr2", config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_TSS"
	else:
		name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes += "_GENE"

	FDR_thresholds = np.loadtxt(name_of_output_file_with_thresholds_estimated_on_odd_even_chromosomes, dtype = str, delimiter = "\t")

	thresholds_test_est_FDR_dist_data, thresholds_test_est_FDR_dist, thresholds_test_est_FDR_data = np.zeros_like(config_variables.FDR), np.zeros_like(config_variables.FDR), np.zeros_like(config_variables.FDR)

	for ind, FDR in enumerate(FDR_thresholds_to_plot):

		thresholds_test_est_FDR_dist_data[ind] = FDR_thresholds[(FDR_thresholds[:,0] == comb) * (FDR_thresholds[:,2].astype(float) == FDR), -1][0]
		
	thresholds_test_est_FDR_dist_data = thresholds_test_est_FDR_dist_data.astype(float)

	i_s_t_filt_1, j_s_t_filt_1, i_s_f_filt_1, j_s_f_filt_1 = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, 1., thresholds_test_est_FDR_dist_data[0], chrom)
	i_s_t_filt_2, j_s_t_filt_2, i_s_f_filt_2, j_s_f_filt_2 = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, thresholds_test_est_FDR_dist_data[0], thresholds_test_est_FDR_dist_data[1], chrom)
	i_s_t_filt_3, j_s_t_filt_3, i_s_f_filt_3, j_s_f_filt_3 = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, thresholds_test_est_FDR_dist_data[1], thresholds_test_est_FDR_dist_data[2], chrom)

	#i_s_t_filt_1_dom, j_s_t_filt_1_dom, i_s_f_filt_1_dom, j_s_f_filt_1_dom = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, 1., threshold_1, chrom, domain = True)
	#i_s_t_filt_2_dom, j_s_t_filt_2_dom, i_s_f_filt_2_dom, j_s_f_filt_2_dom = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, threshold_1, threshold_2, chrom, domain = True)
	#i_s_t_filt_3_dom, j_s_t_filt_3_dom, i_s_f_filt_3_dom, j_s_f_filt_3_dom = interactions_above_threshold(posterior_correl_dist_true, posterior_correl_dist_false, threshold_2, threshold_3, chrom, domain = True)




	import matplotlib.pyplot as plt
	import matplotlib.patches as patch_figures
	import numpy as np
	from matplotlib.collections import PatchCollection
	import matplotlib

	name_of_enhancer_file_for_overlap = config_variables.name_of_enhancer_file_for_overlap
	name_of_time_series_promoter_file_for_TSS_start = config_variables.name_of_time_series_promoter_file_for_TSS_start
	dataset_time_series_dict = config_variables.dataset_time_series_dict
	link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
	dict_chrom_distant = config_variables.dict_chrom_distant
	dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived
	dict_chrom_enh_survived = config_variables.dict_chrom_enh_survived

	#------------------------------------------------------------------------------------- labels part

	enhancer_chroms, enhancer_coordinates, ts_1 = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"]["ER"]]
	promoter_chroms, promoter_coordinates, ts_2 = dataset_time_series_dict[link_data_set_name_to_file_name["promoters"]["ER"]]

	gene_names = np.loadtxt(name_of_time_series_promoter_file_for_TSS_start, dtype = str)

	

	chrom_distant_enh = dict_chrom_distant[chrom]
	chrom_enh_survived = dict_chrom_enh_survived[chrom]
	chrom_pro_survived = dict_chrom_pro_survived[chrom]


	chr_genes = gene_names[gene_names[:, 0] == chrom]
	chr_enhancers = enhancer_coordinates[enhancer_chroms == chrom]

	level_of_enhancer = 5
	higth_of_enhancer = .5
	level_of_promoter = 0
	higth_of_promoter = 1.

	def simpleaxis(ax):
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.get_xaxis().tick_bottom()
		ax.tick_params(
		axis='y',          # changes apply to the x-axis
		which='both',      # both major and minor ticks are affected
		left='off',      # ticks along the bottom edge are off
		right='off',         # ticks along the top edge are off
		labelleft='off') # labels along the bottom edge are off

		

	fig = plt.figure()
	ax = fig.add_subplot(111)
	simpleaxis(ax)
	ax.locator_params(nbins=4)
	x_limits = [4.175*10**7, 4.45*10**7]
	plt.xlim(x_limits)
	
	plt.figtext(0.02, 0.63+0.02,'enhancer', fontsize=20)
	plt.figtext(0.02, 0.21+0.02,'extended gene', fontsize=20)
	plt.figtext(0.125, 0.00, chrom, fontsize=20)
	plt.figtext(0.9215-0.025, 0.24, "+ strand", fontsize=20)
	plt.figtext(0.9215-0.0275, 0.19, " -  strand", fontsize=20)#0.925-0.025

	x = [min([int(chr_genes[:,1][0]), int(chr_enhancers[:,1][0])]) - 1000, max(int(chr_genes[:,2][-1]), int(chr_enhancers[:,1][-1]))  + 1000]
	line, = ax.plot(x, np.zeros_like(x), lw = 1.0, color='black')
	#line, = ax.plot(x, np.zeros_like(x) + 40., lw = 1.5, color='blue')
	line, = ax.plot(x, np.zeros_like(x) + level_of_enhancer, lw = 1.0, color='black')

	filtered_ehnacers_mask_survived_count = np.zeros_like(enhancer_chroms).astype(bool)
	filtered_ehnacers_mask_survived_dist = np.zeros_like(enhancer_chroms).astype(bool)
	filtered_ehnacers_mask_survived_count[chrom_enh_survived] = True
	filtered_ehnacers_mask_survived_dist[chrom_distant_enh] = True
	filtered_ehnacers_mask_survived = filtered_ehnacers_mask_survived_count*filtered_ehnacers_mask_survived_dist

	filtered_promoters_mask_survived = np.zeros_like(promoter_chroms).astype(bool)
	filtered_promoters_mask_survived[dict_chrom_pro_survived[chrom]] = True
	for gene in gene_names[filtered_promoters_mask_survived*(promoter_chroms == chrom)]:
		coords, gene_name, strand = gene[1:3].astype(float), gene[3], gene[4]

		if strand == '+': 
			el = patch_figures.Rectangle((coords[0], level_of_promoter), np.diff(coords)[0], higth_of_promoter, facecolor='b', alpha=0.8)
			#ax.text(coords[0] + np.diff(coords)[0]/2., level_of_promoter + higth_of_enhancer + .5, gene_name, fontsize=5)
		else: 
			el = patch_figures.Rectangle((coords[0], level_of_promoter-higth_of_promoter), np.diff(coords)[0], higth_of_promoter, facecolor='b', alpha=0.5)
			#ax.text(coords[0] + np.diff(coords)[0]/2., level_of_promoter - higth_of_enhancer - .5, gene_name, fontsize=5)
		ax.add_artist(el)
		el.set_clip_box(ax.bbox)
		


	for enhancer in enhancer_coordinates[filtered_ehnacers_mask_survived*(enhancer_chroms == chrom)]:
		coords = enhancer.astype(float)
		rect = patch_figures.Ellipse((coords[0] + np.diff(coords)[0]/2., level_of_enhancer), np.diff(coords)[0], higth_of_enhancer, facecolor='orange', alpha=0.5)
		ax.add_artist(rect)
		rect.set_clip_box(ax.bbox)


	interactions_of_threshold = [[i_s_t_filt_1, j_s_t_filt_1, i_s_f_filt_1, j_s_f_filt_1],
	[i_s_t_filt_2, j_s_t_filt_2, i_s_f_filt_2, j_s_f_filt_2],
	[i_s_t_filt_3, j_s_t_filt_3, i_s_f_filt_3, j_s_f_filt_3]]

	for alpha_, [i_s_t_filt_c_d, j_s_t_filt_c_d, i_s_f_filt_c_d, j_s_f_filt_c_d] in zip([1., 0.66, 0.33], interactions_of_threshold):

		for inter in zip(gene_names[i_s_t_filt_c_d], enhancer_coordinates[j_s_t_filt_c_d]):
			promoter_enhancer_interacts_with, enhancer = inter[0], inter[1]

			coords_promoter_enhancer_interacts_with, strand = promoter_enhancer_interacts_with[[1,2]].astype(float), promoter_enhancer_interacts_with[4]
			enhancer_centre = enhancer[0] + np.diff(enhancer)/2.

			if strand == "-": 
				coords_promoter_enhancer_interacts_with = coords_promoter_enhancer_interacts_with[[1,0]]
				level_of_promoter_strand_specific = level_of_promoter - higth_of_promoter
				#colour = "Silver"
			else:
				level_of_promoter_strand_specific = level_of_promoter + higth_of_promoter
				#colour = "Silver"

			line, = ax.plot([enhancer_centre, coords_promoter_enhancer_interacts_with[0]], [level_of_enhancer - higth_of_enhancer/2., level_of_promoter_strand_specific], lw = 3.0, color = 'green', alpha = alpha_)

		for inter in zip(gene_names[i_s_f_filt_c_d], enhancer_coordinates[j_s_f_filt_c_d]):
			promoter_enhancer_interacts_with, enhancer = inter[0], inter[1]

			coords_promoter_enhancer_interacts_with, strand = promoter_enhancer_interacts_with[[1,2]].astype(float), promoter_enhancer_interacts_with[4]
			enhancer_centre = enhancer[0] + np.diff(enhancer)/2.

			if strand == "-": 
				coords_promoter_enhancer_interacts_with = coords_promoter_enhancer_interacts_with[[1,0]]
				level_of_promoter_strand_specific = level_of_promoter - higth_of_promoter
				colour = "gray"
			else:
				level_of_promoter_strand_specific = level_of_promoter + higth_of_promoter
				colour = "gray"

			line, = ax.plot([enhancer_centre, coords_promoter_enhancer_interacts_with[0]], [level_of_enhancer - higth_of_enhancer/2., level_of_promoter_strand_specific], lw = 3.0, color = colour, alpha = alpha_)

	#plt.xlim([int(chr_genes[:,1][80]), int(chr_enhancers[:,2][133])])
	ax.set_ylim(higth_of_enhancer -2, higth_of_promoter + 7)

	from matplotlib.backends.backend_pdf import PdfPages

	name_of_output_cartoon_file = 'cartoon_{0}_{1}_{2}_{3}_average_PolII'.format(",".join([comb]), config_variables.one_sided_or_two_sided, config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_cartoon_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation: 

		name_of_output_cartoon_file += "_TSS" 
	else:
		name_of_output_cartoon_file += "_GENE"

	pdf = PdfPages(name_of_output_cartoon_file + ".pdf")

	pdf.savefig()
	pdf.close(); plt.show()


