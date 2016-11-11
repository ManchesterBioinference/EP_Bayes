def executor(prior_elements, low_dist, up_dist):

	import config_variables
	import itertools
	mode = config_variables.mode
	chroms_in_prior = config_variables.chroms_in_prior
	np = config_variables.np
	log_distances = config_variables.log_distances

	prob={}
	bins={}
	optimum = {}
	prob_smooth={}
	bins_smooth={}

	attribute_ = {}

	for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):

		prob[positive_or_negative_side] = {}
		bins[positive_or_negative_side] = {}

		attribute_[positive_or_negative_side] = {}

		prob_smooth[positive_or_negative_side] = {}
		bins_smooth[positive_or_negative_side] = {}

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:

			prob[positive_or_negative_side][classification_of_interactions] = prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_frequencies"]
			bins[positive_or_negative_side][classification_of_interactions] = prior_elements[mode][classification_of_interactions]["distance"][positive_or_negative_side]["prior_bins"]
			
			total_array = [prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"][chrom_] for chrom_ in chroms_in_prior]
			total_array = np.array(list(itertools.chain.from_iterable(total_array)))

			attribute_[positive_or_negative_side][classification_of_interactions] = total_array[sign*total_array > 0]

	optimum = "scott"
	colour = ("g", "y")
	number_of_bins = config_variables.number_of_bins
	mode = config_variables.mode

	def calculate_kern(sample_, l_limit, up_limit, number_of_bins):
		import kern_density_est

		xgrid = [[],[]]
		xgrid[0] = np.linspace(l_limit, up_limit, number_of_bins[0])

		prob_, bins_ = kern_density_est.kern_scipy_gaus(np.r_[sample_["positive_interactions"], sample_["negative_interactions"]], colour[0], xgrid[0], bandwidth=optimum)

		return prob_, bins_
	
	for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):
		print positive_or_negative_side

		prob_smooth[positive_or_negative_side], bins_smooth[positive_or_negative_side] = calculate_kern(attribute_[positive_or_negative_side], low_dist[positive_or_negative_side], up_dist[positive_or_negative_side], number_of_bins)
				

	prior_elements[mode]["MOG_distance"] = {}
	for positive_or_negative_side in ["positive_side", "negative_side"]:

		prior_elements[mode]["MOG_distance"][positive_or_negative_side] = {}	

		prob_smooth_ = prob_smooth[positive_or_negative_side][:-1] + np.diff(prob_smooth[positive_or_negative_side])/2.
		prob_smooth_ /= sum(prob_smooth_*np.diff(bins_smooth[positive_or_negative_side]))

		[prior_elements[mode]["MOG_distance"][positive_or_negative_side]["prior_frequencies"], 
		prior_elements[mode]["MOG_distance"][positive_or_negative_side]["prior_bins"]] = prob_smooth_, bins_smooth[positive_or_negative_side] # make the histograms smooth


	def join_two_sides_of_prior_together():
	
		prob = {}
		bins = {}	
		for positive_or_negative_side in ["positive_side", "negative_side"]:
			prob[positive_or_negative_side], bins[positive_or_negative_side] = [prior_elements[mode]["MOG_distance"][positive_or_negative_side]["prior_frequencies"],
																					prior_elements[mode]["MOG_distance"][positive_or_negative_side]["prior_bins"]]

		bins_ = np.r_[bins["negative_side"], bins["positive_side"]]

		if log_distances: 
			prob_ = np.r_[prob["negative_side"], [0], prob["positive_side"]]

		else:
			prob_ = np.r_[prob["negative_side"], [(prob["negative_side"][-1] + prob["positive_side"][0])/2.], prob["positive_side"]]

		prob_ /= sum(prob_*np.diff(bins_))

		return prob_, bins_

	
	prior_elements[mode]["MOG_distance"]["prior_frequencies"], prior_elements[mode]["MOG_distance"]["prior_bins"] = join_two_sides_of_prior_together()

	return prior_elements[mode]["MOG_distance"]["prior_frequencies"], prior_elements[mode]["MOG_distance"]["prior_bins"]


