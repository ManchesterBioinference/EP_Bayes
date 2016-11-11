def execute(prior_elements, plot_atr, plot_atr_kernel):

	import config_variables
	import matplotlib.pyplot as plt
	import numpy as np
	from matplotlib.backends.backend_pdf import PdfPages
	results_folder = config_variables.results_folder
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided
	mode = config_variables.mode
	dataset_names_option = config_variables.dataset_names_option

	if config_variables.mode_of_code == "FULL":
		name_of_histogram_file = results_folder + 'multipage_priors_average_{0}_{1}'.format(one_sided_or_two_sided, config_variables.mode_of_code)

	elif config_variables.Sample_MoG_classificator:
		name_of_histogram_file = results_folder + 'multipage_priors_average_{0}_{1}'.format(one_sided_or_two_sided, config_variables.mode_of_code)

	else:
		name_of_histogram_file = results_folder + 'multipage_priors_average_{0}_{1}'.format(one_sided_or_two_sided, "ODD")
		#name_of_histogram_file = 'multipage_priors_average{0}_{1}_{2}'.format(one_sided_or_two_sided, "_".join(chroms_in_prior), "_".join(chroms_to_infer))

	if plot_atr_kernel:
		name_of_histogram_file = name_of_histogram_file + "_kern"
	elif plot_atr:
		name_of_histogram_file = name_of_histogram_file + "_hist"
	else:
		pass

	name_of_histogram_file = name_of_histogram_file + "{0}_{1}_{2}_{3}".format(config_variables.number_of_bins, config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

	if config_variables.disentagled_features_validation: 
		name_of_histogram_file += "_TSS" 
	else:
		name_of_histogram_file += "_GENE"

	pdf = PdfPages(name_of_histogram_file + ".pdf")

	def plot_double_sided_dist_smooth_histograms(x, freq, color, label):

		plt.plot(x, freq, color + '-', lw=3, label = label)
		plt.fill_between(x,freq, where=None, color=color, alpha = 0.2)

	for label, classification_of_interactions, colour in [["positive interactions", "positive_interactions", "g"], ["negative interactions", "negative_interactions", "y"]]:
		if config_variables.Sample_MoG_classificator and classification_of_interactions == "negative_interactions":	continue	
		freq__ = prior_elements[mode][classification_of_interactions]["distance"]["prior_frequencies"]
		bins__ = prior_elements[mode][classification_of_interactions]["distance"]["prior_bins"]

		plot_double_sided_dist_smooth_histograms(bins__[:-1] + np.diff(bins__)/2, freq__, colour, label)

		plt.figure(1, figsize=(8, 6), dpi=200)
		plt.title("Distance prior", fontsize=28)
		plt.ylabel('density', fontsize=28)
		plt.xlabel("distance [B]", fontsize=28)
		tick_labels = [8, 5, 0 , 5, 8]
		string_labels = [r"$10^{%2d}$" % (i) for i in tick_labels]
		plt.rc('ytick', labelsize = 20)
		plt.rc('xtick', labelsize = 20)
		plt.xticks([-8., -5., 0., 5., 8.],  string_labels, fontsize=23)#["a", "b", "c", "d", "e"])#
		plt.xlim([-8.5, 8.5])
		plt.subplots_adjust(left=0.12, right=0.95, top=0.925, bottom=0.135)


	
	x1,x2,y1,y2 = plt.axis(); plt.axis((x1,x2,0.,y2*1.2)); plt.legend(); pdf.savefig();

	for i, data_set_name in enumerate(dataset_names_option):	
		print data_set_name
		for label, classification_of_interactions, colour in [["positive interactions", "positive_interactions", "g"], ["negative interactions", "negative_interactions", "y"]]:
			if config_variables.Sample_MoG_classificator and classification_of_interactions == "negative_interactions":	continue
			freq__ = prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_frequencies"]
			bins__ = prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_bins"]

			print len(freq__), len(bins__)
			plt.figure(i+2, figsize=(8, 6), dpi=200)			
			if data_set_name == "ER": plt.title(u'ER-\u03B1', fontsize=28)
			else: plt.title(data_set_name, fontsize=28)
			plt.ylabel('density', fontsize=28)
			plt.xlabel('correlation', fontsize=28)
			plot_double_sided_dist_smooth_histograms(bins__[:-1] + np.diff(bins__)/2, freq__, colour, label)
			ax = plt.gca(); cur_ylim = ax.get_ylim(); ax.set_ylim([0, cur_ylim[1]])
			#plt.subplots_adjust(left=0.1, right=0.95, top=0.925, bottom=0.125)
			plt.subplots_adjust(left=0.12, right=0.95, top=0.925, bottom=0.135)
		x1,x2,y1,y2 = plt.axis(); plt.axis((x1,x2,0.,y2*1.2)); plt.legend(); pdf.savefig();

	pdf.close(); plt.close('all')
