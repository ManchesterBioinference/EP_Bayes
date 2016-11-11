def convergence_checker(number_of_samples_arr, burn_in_start):

	import iter_loadtxt
	import numpy as np
	import matplotlib as mpl
	mpl.use('Agg')
	from matplotlib.backends.backend_pdf import PdfPages
	import config_variables
	import os
	results_folder = config_variables.results_folder
	chrom_names = config_variables.chrom_names
	dataset_time_series_dict = config_variables.dataset_time_series_dict
	link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
	import matplotlib.pyplot as plt
	from pylab import rcParams




	#blue_line = mlines.Line2D([], [], color='blue', marker='^', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data+prior')

	#name_std = 'MOG_distance_emprirical_mu_trace_of_std_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
	#name_mean = 'MOG_distance_emprirical_mu_trace_of_mean_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
	#output_mean = open(save_to_folder + name_mean, 'w')
	#output_std = open(save_to_folder + name_std, 'w')
	#mode_of_sampler = ["distance_prior", "distance_MOG", "dirichlet_MOG", "distance_MOG_empir_mu"][3]

	import selected_combinations as sel
	combinations, selected_combinations = sel.selected_combinations("SELECTIVE")

	#number_of_samples = 61000
	chains = [1,2,3]
	#number_of_samples_arr = [[62000]*3, [62000]*3, [240000]*3, [240000]*3, [1200000]*3]
	#number_of_samples_arr = [[80000]*3, [80000]*3, [80000]*3, [80000]*3, [160000]*3]
	min_number_of_samples = np.min(np.array(number_of_samples_arr))
	#burn_in_start = (np.array(number_of_samples_arr)/2.).astype(int)

	#burn_in_start[-1] = [300000]*3

	def opener(option_correl, chrom, number_of_samples_ar):

		comb = "_".join([config_variables.dict_option[el] for el in option_correl])
		kappa_0, mu_0, alpha_0, Beta_0 = config_variables.kappa_0, config_variables.mu_0, config_variables.alpha_0, config_variables.Beta_0

		save_to_folder = os.getcwd() + "/MOG_results_/"

		parameter_multiple_chains = {}

		for chain_number, number_of_samples in zip(chains, number_of_samples_ar):

			name = 'MOG_distance_emprirical_mu_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)

			name = save_to_folder + name
			print name
			parameter_multiple_chains[chain_number] = iter_loadtxt.iter_loadtxt(name, ",", dtype = int, skiprows=min(burn_in_start[index_opt]))

		n_ = parameter_multiple_chains[1].shape[0]
		j_ = parameter_multiple_chains[1].shape[1]
		m_ = len(chains)

		tensor_parameter_multiple_chains = np.zeros((m_, n_, j_))

		for chain_number in chains: tensor_parameter_multiple_chains[chain_number-1] = parameter_multiple_chains[chain_number]

		return tensor_parameter_multiple_chains


	def Gelman_Rubin_R_hat(tensor_parameter_multiple_chains):

		m, n, j = tensor_parameter_multiple_chains.shape
	
		means_j = tensor_parameter_multiple_chains.mean(1)

		def within():

			S_j_squared = np.var(tensor_parameter_multiple_chains, ddof=1, axis=1) #((tensor_parameter_multiple_chains-means_j.reshape(m,1,j))**2).sum(1)/float(n-1) #
			W = S_j_squared.mean(0)
			return W

		def between():

			mean_of_means = means_j.mean(0)
			B = np.var(means_j, ddof = 1, axis = 0)*n#n/float(m-1)*((means_j-mean_of_means)**2).sum(0) #
			return B

		W = within()
		B = between()

		Var_par = (1-n**(-1))*W + n**(-1)*B

		R_hat = (Var_par/W)**0.5

		return R_hat


	def inter_enhancer(chrom):

		pro_chroms, pro_coords, pro_time_series = dataset_time_series_dict[link_data_set_name_to_file_name["promoters"]["ER"]]
		enh_chroms, enh_coords, enh_time_series = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"]["ER"]]

		filtered_enhancers = config_variables.filtered_enhancers
		proximal_enhancers_mask = config_variables.proximal_enhancers_mask

		chrom_enh_survived = np.where((enh_chroms == chrom)*np.invert(proximal_enhancers_mask)*filtered_enhancers)[0]

		negative_interactions = config_variables.negative_interactions
		from  prepare_interactions_clean import un_string

		indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

		if config_variables.disentagled_features_validation: 
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
		else:
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

		true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

		i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

		interacting_enhancers_ = np.unique(j_s_t)-total_e

		enhancer_array_survived = np.zeros(len(indexes_e), bool)
		enhancer_array_interacting = np.zeros(len(indexes_e), bool)

		enhancer_array_survived[chrom_enh_survived-total_e] = True
		enhancer_array_interacting[interacting_enhancers_] = True

		mask_interacting_c = np.in1d(np.where(enhancer_array_survived)[0], np.where(enhancer_array_interacting)[0])

		return np.where(mask_interacting_c)[0]


	#option_correl_select = selected_combinations[1]
	#chrom_ = chrom_names[1]

	R_hat_dict = {}

	name_of_r_hat = results_folder + "r_hat_for_{0}_{1}".format("_".join(burn_in_start[:, 0].astype(str)), min_number_of_samples)
	pdf = PdfPages(name_of_r_hat + ".pdf")

	plt.rcParams['xtick.labelsize'] = 26
	plt.rc('ytick', labelsize = 26)
	f, ax = plt.subplots(1, len(selected_combinations), sharex=True, sharey=True, figsize=(20,10))
	f.subplots_adjust(left=0.085, bottom=0.15, right=0.965, top=0.925, hspace=0.1, wspace=0.05)
	marker_size_legend_box = 19
	legend_box_names_font_size = 22
	size_of_combination_name = 35
	size_of_y_label = 35
	ax[0].set_ylabel('Density', fontsize = size_of_y_label)

	ax[0].locator_params(axis = 'x', nbins = 6)
	f.text(0.525, 0.04, r'$ \hat R $', ha='center', fontsize=35)

	#number_of_samples_arr = [120000]*3

	for index_opt, option_correl_select, number_of_samples_ar in zip(range(len(selected_combinations)), selected_combinations, number_of_samples_arr):
		#if index_opt == 4: continue

		comb_plot = ",".join([config_variables.dict_option[el] for el in option_correl_select])
		if option_correl_select == combinations[-1]: comb_plot = "All"
		ax[index_opt].set_title(comb_plot, fontsize = size_of_combination_name)
	
		R_hat_total = []
		for chrom_ in chrom_names[:-1]:

			tensor_parameter_multiple_chains = opener(option_correl_select, chrom_, number_of_samples_ar)
			#ll
			if config_variables.interacting_enhancers_only_MOG: 
				interacting_enhancers = inter_enhancer(chrom_)
				tensor_parameter_multiple_chains = tensor_parameter_multiple_chains[:, :, interacting_enhancers]
		
			R_hat = Gelman_Rubin_R_hat(tensor_parameter_multiple_chains[:,:,:])
	
			print np.where(np.isnan(R_hat))[0]
			R_hat[np.isnan(R_hat)] = 1.
			R_hat_total += R_hat.tolist()
	
		R_hat_dict[index_opt] = R_hat_total

		ax[index_opt].hist(R_hat_dict[index_opt], bins = np.arange(0.95, 1.5, 0.01), facecolor='green', alpha=0.5, normed=1, histtype='bar')
		ax[index_opt].set_xlim([0.9, 1.6])
		ax[index_opt].set_ylim([0., 90.])

	pdf.savefig()
	pdf.close()
	plt.close("all")

#def opener_likelihood(option_correl, chrom, number_of_samples_ar):
#
#	comb = "_".join([config_variables.dict_option[el] for el in option_correl])
#	kappa_0, mu_0, alpha_0, Beta_0 = config_variables.kappa_0, config_variables.mu_0, config_variables.alpha_0, config_variables.Beta_0
#
#	save_to_folder = os.getcwd() + "/MOG_results_/"
#
#	parameter_multiple_chains_likelihood = {}
#
#	for chain_number, number_of_samples in zip(chains, number_of_samples_ar):
#
#		name_likelihood = 'MOG_distance_emprirical_mu_trace_of_likelihood_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
#
#		name = save_to_folder + name_likelihood
#		print name
#		parameter_multiple_chains_likelihood[chain_number] =  np.ravel(iter_loadtxt.iter_loadtxt(name, ",", dtype = float))[burn_in_start:min_number_of_samples]

#	n_ = parameter_multiple_chains_likelihood[1].shape[0]
#	m_ = len(chains)
#	tensor_parameter_multiple_chains_likelihood = np.zeros((m_, n_))
#
#	for chain_number in chains: tensor_parameter_multiple_chains_likelihood[chain_number-1] = parameter_multiple_chains_likelihood[chain_number]
#
#	return tensor_parameter_multiple_chains_likelihood
#
#
#tensor_parameter_multiple_chains_likelihood_dict = {}
#number_of_samples_arr = [[61000, 70000, 60000], [61000, 70000, 60000], [120000]*3, [120000]*3, [150000]*3]

#name_of_likelihood = results_folder + "likelihood_for_{0}_{1}".format(burn_in_start, min_number_of_samples)
#pdf_likelihood = PdfPages(name_of_likelihood + ".pdf")
#plt.rcParams['xtick.labelsize'] = 26
#plt.rc('ytick', labelsize = 26)
#f, ax = plt.subplots(1, len(selected_combinations), sharex=True, sharey=True, figsize=(20,10))
#f.subplots_adjust(left=0.085, bottom=0.15, right=0.965, top=0.925, hspace=0.1, wspace=0.05)
#marker_size_legend_box = 19
#legend_box_names_font_size = 22
#size_of_combination_name = 35
#size_of_y_label = 35
#ax[0].set_ylabel('Density', fontsize = size_of_y_label)
#ax[0].locator_params(axis = 'x', nbins = 6)
#f.text(0.525, 0.04, "Likelihood", ha='center', fontsize=35)


#def movingaverage (values, window):
#    weights = np.repeat(1.0, window)/window
#    sma = np.convolve(values, weights, 'valid')
#    return sma
#
#pdf_likelihood = PdfPages(name_of_likelihood + "_3.pdf")
#for index_opt, option_correl_select, number_of_samples_ar in zip(range(len(selected_combinations)), selected_combinations, number_of_samples_arr)[2:3]:
#
#	#f, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20,10))
#	index_opt = 0
#
#	tensor_parameter_multiple_chains_likelihood_dict[index_opt] = {}
#
#	comb_plot = ",".join([config_variables.dict_option[el] for el in option_correl_select])
#	if option_correl_select == combinations[-1]: comb_plot = "All"
#	#ax[index_opt].set_title(comb_plot, fontsize = size_of_combination_name)
#
#	small_burn_in = 1000	
#
#	for chrom_ in chrom_names[:-1]:
#
#		f, ax = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(20,10))
#		ax.set_title(comb_plot, fontsize = size_of_combination_name)
#
#		tensor_parameter_multiple_chains_likelihood_dict[index_opt][chrom_] = opener_likelihood(option_correl_select, chrom_, number_of_samples_ar)
#
#		for chain_number in chains:	
#
#			mins = []
#			maxs = []
#
#			N = np.arange(min_number_of_samples-burn_in_start) + 1
#
#			every_nth_sample = tensor_parameter_multiple_chains_likelihood_dict[index_opt][chrom_][chain_number - 1]
#
#			cum_mean_of_every_nth_sample = every_nth_sample[small_burn_in:].cumsum()/N[:-small_burn_in]
#
#			frame_ = 500
#			moving_average_n = movingaverage(every_nth_sample, frame_)
#
#			ax.plot(N[frame_/2:-frame_/2], moving_average_n[:-1], alpha = 0.5)
#			ax.plot(N[small_burn_in/2:-small_burn_in/2], cum_mean_of_every_nth_sample, alpha = 0.8)
#
#			mins += [np.min(cum_mean_of_every_nth_sample)]
#			maxs += [np.max(cum_mean_of_every_nth_sample)]
#
#		print min(mins), max(maxs)
#
#		#ax.set_xlim([burn_in_start, min_number_of_samples])
#		#ax.set_ylim([min(mins), max(maxs)])
#
#		pdf_likelihood.savefig()
#		plt.close()
#
#pdf_likelihood.close()
#plt.close("all")

