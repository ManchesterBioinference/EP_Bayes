def prior_bins_prob_and_plotter(prior_elements, low_dist, up_dist, use_smooth_prior_for_estimation, plot_atr, plot_atr_kernel, Sample_MoG_classificator =  False):
	import matplotlib as mpl
	#mpl.use('Agg')
	import matplotlib.pyplot as plt
	from matplotlib.backends.backend_pdf import PdfPages
	import config_variables
	domain = config_variables.domain
	mode = config_variables.mode
	chroms_in_prior = config_variables.chroms_in_prior
	np = config_variables.np
	dataset_names_option = config_variables.dataset_names_option
	low_cor, up_cor = -1., 1.
	one_sided_or_two_sided = config_variables.one_sided_or_two_sided	
	log_distances = config_variables.log_distances
	likelihood_cross_validation = config_variables.likelihood_cross_validation
	distant_enh_only = config_variables.distant_enh_only
	interacting_enhancers_only = config_variables.interacting_enhancers_only
	chroms_to_infer=config_variables.chroms_to_infer
	results_folder = config_variables.results_folder

	def bins_prep_adaptive(array, l_limit, u_limit, how_many_in_bin):

		array_sorted=sorted(np.concatenate((array, [l_limit], [u_limit])))
		bins=[]
		bins = [array_sorted[n] for n in range(len(array_sorted)) if (n % how_many_in_bin)==0]
	
		if bins[-1] <> array_sorted[-1]: bins.append(array_sorted[-1])	

		return bins


	def profile_histogram_adaptive_domains(array, l_limit, u_limit, how_many_in_bin, possible_distances_counts):#, f, colour):
		"""
		gives the empirical probabilities of the interactions based on the mean interactions in "IHH015M_ipet.tsv"
		"""

		bins = bins_prep_adaptive(array, l_limit, u_limit, how_many_in_bin)

		distance_allocations = np.digitize(array, bins=bins)
		distance_allocations = distance_allocations - 1
		probabilities = np.zeros(len(bins)-1)
		norm = np.sum(1/possible_distances_counts.astype(float))

		for index, el in enumerate(np.unique(distance_allocations)): probabilities[index] = np.sum(1/possible_distances_counts[el == distance_allocations].astype(float))

		differences = np.diff(bins)

		prob = probabilities/norm/differences

		#if one_sided_or_two_sided == "single_sided": 
		#	prob, bins = two_sided_to_one_sided_domain(prob, bins, number_of_samples = 800000, number_of_new_bins = 600)

		return prob, bins		

	def profile_histogram_adaptive(array, l_limit, u_limit, how_many_in_bin):
		"""
		gives the empirical probabilities of the interactions based on the mean interactions in "IHH015M_ipet.tsv"
		"""
		bins = bins_prep_adaptive(array, l_limit, u_limit, how_many_in_bin)

		n, bins = np.histogram(array, bins=bins, density = True)
		return n, bins



	import itertools

	def create_priors_domains():

		if not(domain): 
			if mode == 'promoter_enhancer_interactions': 
				if interacting_enhancers_only: initiatie_number_of_bins = iter([25, 20, 2000, 2000])
				else: initiatie_number_of_bins = iter([25, 20, 2000, 2000])

			else: 
				if interacting_enhancers_only: initiatie_number_of_bins = iter([80, 70, 20000, 20000])
				else: initiatie_number_of_bins = iter([80, 70, 200000, 200000])

		else: 
			if mode == 'promoter_enhancer_interactions': initiatie_number_of_bins = iter([25, 20, 40, 50])
			else: initiatie_number_of_bins = iter([70, 70, 200, 200])
			

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			for attribute_of_interaction in ["distance", "correlation"]:

				number_in_bin = initiatie_number_of_bins.next()
				prior_elements[mode][classification_of_interactions][attribute_of_interaction]["number_in_bin_of_histogram"] = number_in_bin

				if attribute_of_interaction == "distance":
					total_array = [prior_elements[mode][classification_of_interactions][attribute_of_interaction]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
					total_array = np.array(list(itertools.chain.from_iterable(total_array)))
					if not(domain):
				
						if one_sided_or_two_sided == "double_sided":

							for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):

								prior_elements[mode][classification_of_interactions][attribute_of_interaction][positive_or_negative_side] = {}

								[prior_elements[mode][classification_of_interactions][attribute_of_interaction][positive_or_negative_side]["prior_frequencies"],
								prior_elements[mode][classification_of_interactions][attribute_of_interaction][positive_or_negative_side]["prior_bins"]] = profile_histogram_adaptive(total_array[sign*total_array > 0], l_limit = low_dist[positive_or_negative_side], u_limit = up_dist[positive_or_negative_side], how_many_in_bin = int(number_in_bin/2.))

	
						else:
							[prior_elements[mode][classification_of_interactions][attribute_of_interaction]["prior_frequencies"], 
							prior_elements[mode][classification_of_interactions][attribute_of_interaction]["prior_bins"]] = profile_histogram_adaptive(total_array, low_dist, up_dist, number_in_bin)
					
					else:
						possible_distances_counts = prior_elements[mode][classification_of_interactions][attribute_of_interaction]["possible_distances_counts"]
						possible_distances_counts = np.array(possible_distances_counts)[possible_distances_counts <> 0]
						total_array = np.array(total_array)[possible_distances_counts <> 0]
						if one_sided_or_two_sided == "double_sided":
							for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):
								prior_elements[mode][classification_of_interactions][attribute_of_interaction][positive_or_negative_side] = {}

								[prior_elements[mode][classification_of_interactions][attribute_of_interaction][positive_or_negative_side]["prior_frequencies"], 
								prior_elements[mode][classification_of_interactions][attribute_of_interaction][positive_or_negative_side]["prior_bins"]] = profile_histogram_adaptive_domains(total_array[sign*total_array > 0], low_dist[positive_or_negative_side], up_dist[positive_or_negative_side], int(number_in_bin/2.), possible_distances_counts[sign*total_array > 0])
						elif one_sided_or_two_sided == "single_sided":
								[prior_elements[mode][classification_of_interactions][attribute_of_interaction]["prior_frequencies"], 
								prior_elements[mode][classification_of_interactions][attribute_of_interaction]["prior_bins"]] = profile_histogram_adaptive_domains(total_array, low_dist, up_dist, number_in_bin, possible_distances_counts)


				if attribute_of_interaction == "correlation":
					for data_set_name in dataset_names_option:
						total_array = [prior_elements[mode][classification_of_interactions][attribute_of_interaction][data_set_name]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
						total_array = list(itertools.chain.from_iterable(total_array))
						[prior_elements[mode][classification_of_interactions][attribute_of_interaction][data_set_name]["prior_frequencies"], 
						prior_elements[mode][classification_of_interactions][attribute_of_interaction][data_set_name]["prior_bins"]] = profile_histogram_adaptive(total_array, low_cor, up_cor, number_in_bin)


	create_priors_domains()


	def join_two_sides_of_prior_together(classification_of_interactions):
	
		prob = {}
		bins = {}	
		for positive_or_negative_side in ["positive_side", "negative_side"]:
			prob[positive_or_negative_side], bins[positive_or_negative_side] = [prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_frequencies"],
																				prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_bins"]]

		bins_ = np.r_[bins["negative_side"], bins["positive_side"]]

		if log_distances: 
			prob_ = np.r_[prob["negative_side"], [0], prob["positive_side"]]

		else:
			prob_ = np.r_[prob["negative_side"], [(prob["negative_side"][-1] + prob["positive_side"][0])/2.], prob["positive_side"]]

		prob_ /= sum(prob_*np.diff(bins_))

		return prob_, bins_

	
	def calculate_kern(attribute_of_interaction_, sample_, l_limit, up_limit, number_of_bins, colour = ("r", "b"), weights_ = None, bandwidth_pos = None, bandwidth_neg = None):
	
		prob_ = {}
		bins_ = {}
		import kern_density_est
		kern_density_est.plot_atr = plot_atr_kernel

		xgrid = [[],[]]
		xgrid[0] = np.linspace(l_limit, up_limit, number_of_bins[0])
		xgrid[1] = np.linspace(l_limit, up_limit, number_of_bins[1])

		if domain:
	
			if attribute_of_interaction_ == "distance":

				prob_["positive_interactions"], bins_["positive_interactions"] = kern_density_est.kern_scipy_gaus_weighted(sample_["positive_interactions"], colour[0], xgrid[0], weights = weights_["positive_interactions"], bandwidth = bandwidth_pos)
				prob_["negative_interactions"], bins_["negative_interactions"] = kern_density_est.kern_scipy_gaus_weighted(sample_["negative_interactions"], colour[1], xgrid[1], weights = weights_["negative_interactions"], bandwidth = bandwidth_neg)

			else:

				prob_["positive_interactions"], bins_["positive_interactions"] = kern_density_est.kern_scipy_gaus(sample_["positive_interactions"], colour[0], xgrid[0], bandwidth = bandwidth_pos)
				prob_["negative_interactions"], bins_["negative_interactions"] = kern_density_est.kern_scipy_gaus(sample_["negative_interactions"], colour[1], xgrid[1], bandwidth = bandwidth_neg)	
		else:


			prob_["positive_interactions"], bins_["positive_interactions"] = kern_density_est.kern_scipy_gaus(sample_["positive_interactions"], colour[0], xgrid[0], bandwidth=bandwidth_pos)
			prob_["negative_interactions"], bins_["negative_interactions"] = [], []
			if not(Sample_MoG_classificator): prob_["negative_interactions"], bins_["negative_interactions"] = kern_density_est.kern_scipy_gaus(sample_["negative_interactions"], colour[1], xgrid[1], bandwidth=bandwidth_neg)

		if use_smooth_prior_for_estimation:	return  prob_, bins_
		else: return  [[], []], [[], []]


	def calculate_bandwidth_dist(sample_):
		import kern_density_est
		optimum = {}
		if domain: print "positive_interactions"

		if likelihood_cross_validation:
			if domain:
				print "negative_interactions"
				optimum["positive_interactions"] = kern_density_est.chrom_cross_validation_distance(prior_elements, thresholds = np.linspace(0.01, .4, 200), classification_of_interactions = "positive_interactions", weights = weights)
				optimum["negative_interactions"] = kern_density_est.chrom_cross_validation_distance(prior_elements, thresholds = np.linspace(0.01, .4, 200), classification_of_interactions = "negative_interactions", weights = weights)
			else: 
				optimum["positive_interactions"] = kern_density_est.chrom_cross_validation_distance(prior_elements, thresholds = np.linspace(0.01, .4, 200), classification_of_interactions = "positive_interactions")
				optimum["negative_interactions"] = {}
				for positive_or_negative_side in ["positive_side", "negative_side"]: 
					optimum["negative_interactions"][positive_or_negative_side] = "scott"

		else: 
			optimum["positive_interactions"] = {}
			optimum["negative_interactions"] = {}
			for positive_or_negative_side in ["positive_side", "negative_side"]: 
				optimum["positive_interactions"][positive_or_negative_side] = kern_density_est.cross_validation(sample_[positive_or_negative_side]["positive_interactions"])
				optimum["negative_interactions"][positive_or_negative_side] = "scott"

		return optimum

	def calculate_bandwidth_correl(sample_):
		import kern_density_est
		optimum = {}
		
		if likelihood_cross_validation:
			if domain:
				print "negative_interactions"
				optimum["positive_interactions"] = kern_density_est.chrom_cross_validation_correlation(prior_elements, data_set_name, thresholds = np.linspace(0.01, .4, 200), classification_of_interactions = "positive_interactions")
				optimum["negative_interactions"] = kern_density_est.chrom_cross_validation_correlation(prior_elements, data_set_name, thresholds = np.linspace(0.01, .4, 200), classification_of_interactions = "negative_interactions")
			else: 
				optimum["positive_interactions"] = kern_density_est.chrom_cross_validation_correlation(prior_elements, data_set_name, thresholds = np.linspace(0.01, .4, 200), classification_of_interactions = "positive_interactions")
				optimum["negative_interactions"] = "scott"
		else: 
				
			optimum["positive_interactions"] = kern_density_est.cross_validation(sample_["positive_interactions"])
			optimum["negative_interactions"] = "scott"

		return optimum

	prob={}
	bins={}
	optimum = {}
	prob_smooth={}
	bins_smooth={}

	attribute_ = {}
	weights = {}

	number_of_bins = config_variables.number_of_bins

	if not one_sided_or_two_sided == "double_sided":
		print "not implemented - requires recalculation of optimal bandwidths and changes in cross-validation functions, since the cross-validation functions are now written to calculate first one side and than another - it's easy to fix that, the normal valdation (not all but one) is almost implemented "
		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			prob[classification_of_interactions] = prior_elements[mode][classification_of_interactions]["distance"]["prior_frequencies"]
			bins[classification_of_interactions] = prior_elements[mode][classification_of_interactions]["distance"]["prior_bins"]

			total_array = [prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
			total_array = np.array(list(itertools.chain.from_iterable(total_array)))

			
		if use_smooth_prior_for_estimation:
			
			prob_smooth, bins_smooth = calculate_kern(total_array, low_dist, up_dist, number_of_bins, colour = ("g", "y"))

			for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
				prob_smooth_ = prob_smooth[classification_of_interactions][:-1] + np.diff(prob_smooth[classification_of_interactions])/2.
				prob_smooth_ /= sum(prob_smooth_*np.diff(bins_smooth[classification_of_interactions]))

				[prior_elements[mode][classification_of_interactions]["distance"]["prior_frequencies"], 
				prior_elements[mode][classification_of_interactions]["distance"]["prior_bins"]] = prob_smooth_, bins_smooth[classification_of_interactions]
	

	elif one_sided_or_two_sided == "double_sided":		

		if use_smooth_prior_for_estimation:

			for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):

				prob[positive_or_negative_side] = {}
				bins[positive_or_negative_side] = {}

				attribute_[positive_or_negative_side] = {}
				weights[positive_or_negative_side] = {}

				prob_smooth[positive_or_negative_side] = {}
				bins_smooth[positive_or_negative_side] = {}

				for classification_of_interactions in ["positive_interactions", "negative_interactions"]:

					prob[positive_or_negative_side][classification_of_interactions] = prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_frequencies"]
					bins[positive_or_negative_side][classification_of_interactions] = prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_bins"]
			
					total_array = [prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
					total_array = np.array(list(itertools.chain.from_iterable(total_array)))

					attribute_[positive_or_negative_side][classification_of_interactions] = total_array[sign*total_array > 0]

					if domain: 
						possible_distances_counts = prior_elements[mode][classification_of_interactions]["distance"]["possible_distances_counts"]
						possible_distances_counts = np.array(possible_distances_counts)[possible_distances_counts <> 0]
						total_array = np.array(total_array)[possible_distances_counts <> 0]
						weights[positive_or_negative_side][classification_of_interactions] = (1./possible_distances_counts[sign*total_array > 0])/np.sum(1./possible_distances_counts[sign*total_array > 0])
						attribute_[positive_or_negative_side][classification_of_interactions] = total_array[sign*total_array > 0]


			optimum = calculate_bandwidth_dist(attribute_)

			for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):
				print positive_or_negative_side

				prob_smooth[positive_or_negative_side], bins_smooth[positive_or_negative_side] = calculate_kern("distance", attribute_[positive_or_negative_side], low_dist[positive_or_negative_side], up_dist[positive_or_negative_side], number_of_bins, colour = ("g", "y"), weights_ = weights[positive_or_negative_side], bandwidth_pos = optimum["positive_interactions"][positive_or_negative_side], bandwidth_neg = optimum["negative_interactions"][positive_or_negative_side])
				

			for positive_or_negative_side in ["positive_side", "negative_side"]:
				for classification_of_interactions in ["positive_interactions", "negative_interactions"]:				
					if Sample_MoG_classificator and classification_of_interactions == "negative_interactions": continue

					prob_smooth_ = prob_smooth[positive_or_negative_side][classification_of_interactions][:-1] + np.diff(prob_smooth[positive_or_negative_side][classification_of_interactions])/2.
					prob_smooth_ /= sum(prob_smooth_*np.diff(bins_smooth[positive_or_negative_side][classification_of_interactions]))

					[prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_frequencies"], 
					prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_bins"]] = prob_smooth_, bins_smooth[positive_or_negative_side][classification_of_interactions] # make the histograms smooth

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			[prior_elements[mode][classification_of_interactions]["distance"]["prior_frequencies"], 
			prior_elements[mode][classification_of_interactions]["distance"]["prior_bins"]] = join_two_sides_of_prior_together(classification_of_interactions)


	if use_smooth_prior_for_estimation:

		prob={}
		bins={}
		attribute_={}

		number_of_bins = config_variables.number_of_bins

		

		for i, data_set_name in enumerate(dataset_names_option):
			print data_set_name
			for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
				prob[classification_of_interactions] = prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_frequencies"]
				bins[classification_of_interactions] = prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_bins"]

				total_array = [prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
				total_array = np.array(list(itertools.chain.from_iterable(total_array)))

				attribute_[classification_of_interactions] = total_array

			optimum = calculate_bandwidth_correl(attribute_)

			prob_smooth, bins_smooth = calculate_kern("correlation", attribute_, -1., 1., number_of_bins, colour = ("g", "y"), bandwidth_pos = optimum["positive_interactions"], bandwidth_neg = optimum["negative_interactions"])

	
			for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
				if Sample_MoG_classificator and classification_of_interactions == "negative_interactions": continue
				prob_smooth_ = prob_smooth[classification_of_interactions][:-1] + np.diff(prob_smooth[classification_of_interactions])/2.
				prob_smooth_ /= sum(prob_smooth_*np.diff(bins_smooth[classification_of_interactions]))

				[prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_frequencies"], 
				prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["prior_bins"]] = prob_smooth_, bins_smooth[classification_of_interactions]

		
	return prior_elements

