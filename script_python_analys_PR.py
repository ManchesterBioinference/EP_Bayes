def executor(selection_option, FDR_level):
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	import itertools
	import config_variables

	data_folder = "./data/"
	results_folder = config_variables.results_folder

	def calculate_single_ROC_best_True_sensitivity(probabilities_true, probabilities_false, decreasing = True):

		_True_positives_of_threshold = []
		_False_positives_of_threshold = []

		sorted_prob_true = np.sort(probabilities_true)
		sorted_prob_false = np.sort(probabilities_false)

		sorted_thresholds = np.sort(np.unique(np.r_[probabilities_true, probabilities_false]))
		sorted_thresholds = np.unique(np.r_[sorted_thresholds, np.max(sorted_thresholds)*1.01])

		len_prob_true = len(probabilities_true)
		len_prob_false = len(probabilities_false)

		print 'len prob: ', len_prob_true, len_prob_false

		if decreasing:
			_True_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_true, sorted_thresholds)[0][::-1])
			_False_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_false, sorted_thresholds)[0][::-1])
		else: 
			_True_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_true, sorted_thresholds)[0])
			_False_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_false, sorted_thresholds)[0])

		Precision = np.array(_True_positives_of_threshold, dtype = float)/(np.array(_True_positives_of_threshold, dtype = float) + np.array(_False_positives_of_threshold, dtype = float))
	
		True_positive_Rate = np.array(_True_positives_of_threshold)/float(len_prob_true)
		
		False_positive_Rate = np.array(_False_positives_of_threshold)/float(len_prob_false)

		return True_positive_Rate, Precision

	stuff = [1, 2, 3, 4]
	combinations = []

	for L in range(0, len(stuff)+1):
		for subset in itertools.combinations(stuff, L):
			if len(subset): combinations += [list(subset)]

	selected_combinations = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()

	dict_option = dict(zip(range(len(config_variables.datasets_names)), config_variables.datasets_names))

	option_ = selected_combinations[selection_option]
	
	comb = ",".join([dict_option[el] for el in option_])

	name_of_output_file = results_folder + "clusters_genes_vs_counts_prob_distant_all_{0}_{1}_smo_{2}_proximal_version_PR_met".format(FDR_level, ",".join([comb]), config_variables.use_smooth_prior_for_estimation)

	name_of_output_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)


	enhancer_targets = np.loadtxt(name_of_output_file, delimiter = "\t", dtype = str)
	#enhancer_targets = np.loadtxt("clusters_genes_vs_counts_prob_distant_all_0.25_PolII,ER_smo_True_proximal_version_PR_met.txt", delimiter = "\t", dtype = str)

	enhancer_targets_genes = enhancer_targets[:, 0]
	probabilities = enhancer_targets[:, 4]

	probabilities_non_zero = probabilities.astype(float) <> 0.

	distal_enhancer_targets_genes = enhancer_targets_genes[probabilities_non_zero]
	distal_probabilities = probabilities[probabilities_non_zero]

	proximity = enhancer_targets[:,5].astype(float)
	proximal_enhancer_targets_genes = enhancer_targets_genes[(proximity <> 0.)]
	proximal_values = proximity[(proximity <> 0.)]
	#proximal_values[proximal_values == 0.] = proximal_values.max() + 1.
	#proximal_values = proximal_values.min(1)

	def simpleaxis(ax):
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		#ax.spines['bottom'].set_visible(False)
		ax.tick_params(top="off")
		#ax.tick_params(bottom="off")
		ax.tick_params(right="off")
		#ax.tick_params(left="off")
		

	fig = plt.figure()
	ax = fig.add_subplot(111)
	simpleaxis(ax)

	colours = iter(plt.rcParams['axes.color_cycle'])

	for data_set in np.array(["GRO", "RNA"])[[0]]:
		for FDR in ("0,001", "0,05", "0,01"):
			SEQ_genes = np.loadtxt(data_folder + "{0}seq_DE_{1}.csv.gz".format(data_set, FDR), dtype = str)
			#------------distance
			c=colours.next()
			positive_negatives_SEQ = np.in1d(distal_enhancer_targets_genes, SEQ_genes)
			True_positive_Rate_GRO, Precision_GRO = calculate_single_ROC_best_True_sensitivity(distal_probabilities[positive_negatives_SEQ].astype(float), distal_probabilities[np.invert(positive_negatives_SEQ)].astype(float))
			ax.plot(True_positive_Rate_GRO, Precision_GRO, label="{0}seq genes, FDR: {1}, model predictions".format(data_set, FDR), linewidth=2, color = c)
			#------------distance
			#------------proximity
			positive_negatives_SEQ_prox = np.in1d(proximal_enhancer_targets_genes, SEQ_genes)
			True_positive_Rate_GRO_prox, Precision_GRO_prox = calculate_single_ROC_best_True_sensitivity(proximal_values[positive_negatives_SEQ_prox], proximal_values[np.invert(positive_negatives_SEQ_prox)], decreasing = False)
			plt.plot(True_positive_Rate_GRO_prox, Precision_GRO_prox, label="{0}seq genes, FDR: {1}, proximity predictions".format(data_set, FDR), linewidth=3, color = c, linestyle='--')#markersize=8, marker="*")
			#------------proximity

			#gene_names = np.loadtxt("Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered", dtype = str, usecols = (3,))		
			#Percent_of_DE_genes = np.in1d(gene_names, SEQ_genes).sum()/float(len(gene_names))

			#print Percent_of_DE_genes

	plt.rcParams['xtick.labelsize'] = 18
	plt.rc('ytick', labelsize=18)	
	#plt.rc('xtick', labelsize=20)	
	plt.rcParams['figure.figsize'] = 15, 10 #figsize=(20,10)

	plt.xlabel('Recall(TPR)', fontsize=20)
	plt.ylabel('Precision', fontsize=20)
	#plt.title('1-prod(1-p_i) - distal vs proximal')

	plt.legend()

	name_of_output_FDR_file = results_folder + 'GRO_seq_active_GENES_{0}_{1}_{2}_{3}_{4}_average_PolII'.format(FDR_level, ",".join([comb]), config_variables.one_sided_or_two_sided, config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

	name_of_output_FDR_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation: 

		name_of_output_FDR_file += "_TSS" 
	else:
		name_of_output_FDR_file += "_GENE"

	pdf = PdfPages(name_of_output_FDR_file + ".pdf")

	pdf.savefig()
	plt.close("all")
	pdf.close(); #plt.show()


