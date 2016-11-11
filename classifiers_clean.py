import config_variables
import itertools
from prepare_interactions_clean import un_string
import chrom_specific_negative_interactions
chroms_to_infer = config_variables.chroms_to_infer
dataset_names_option = config_variables.dataset_names_option
np = config_variables.np
dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived
classificator_elements = config_variables.classificator_elements
mode = config_variables.mode
filter_value = config_variables.filter_value
alternative_classificator = config_variables.alternative_classificator
import time


def classifier(array_of_merged_probabilities_true, array_of_merged_probabilities_false, K):

	prob = (1 + (K-1)*np.exp((np.log(array_of_merged_probabilities_false)-np.log(array_of_merged_probabilities_true)).sum(1)))**-1		

	return prob


classifiers_elements = {}


#-----------------------------------------------------------------
#initialises global dictionary
for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
	classifiers_elements[classification_of_interactions] = {}
	for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
		classifiers_elements[classification_of_interactions][probability_of_being_positive_or_negative] = {}
		classifiers_elements[classification_of_interactions][probability_of_being_positive_or_negative]["column_stacked_probabilities"] = {}
		classifiers_elements[classification_of_interactions][probability_of_being_positive_or_negative]["posterior_values"] = {}
		for chrom_ in chroms_to_infer:
			classifiers_elements[classification_of_interactions][probability_of_being_positive_or_negative]["column_stacked_probabilities"][chrom_] = []
			classifiers_elements[classification_of_interactions][probability_of_being_positive_or_negative]["posterior_values"][chrom_] = []

#fill in the global dictionary:
#-----------------------------------------------------------------
for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
	for probability_of_being_positive_or_negative in ["probabilities_of_being_positive_interactions", "probabilities_of_being_negative_interactions"]:
		for chrom_ in chroms_to_infer:
			column_stacked_probabilities = classificator_elements[filter_value][mode][classification_of_interactions]["distance"][probability_of_being_positive_or_negative]["posterior_component_values"][chrom_]
			for data_set_name in dataset_names_option: 
				column_stacked_probabilities = np.column_stack((column_stacked_probabilities, classificator_elements[filter_value][mode][classification_of_interactions]["correlation"][probability_of_being_positive_or_negative][data_set_name]["posterior_component_values"][chrom_]))

			classifiers_elements[classification_of_interactions][probability_of_being_positive_or_negative]["column_stacked_probabilities"][chrom_] = column_stacked_probabilities

#calculates the posterior with different options
#-----------------------------------------------------------------

def interactions_extractor(chrom):

	if config_variables.disentagled_features_validation: 
		chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
	else:
		chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

	true_pro_enh_indexes = un_string(chr_interactions_pro_enh[:, :2])

	#true_pro_enh_indexes = un_string(config_variables.chr_interactions_dict_pro_enh[chrom])

		
	prom_enh_false_interactions = chrom_specific_negative_interactions.chrom_specific_negative_interactions(chrom, mode)

	enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e = chrom_specific_negative_interactions.initialise_variables(chrom)
	true_pro_enh_indexes[:,0] = true_pro_enh_indexes[:,0] - total_p
	true_pro_enh_indexes[:,1] = true_pro_enh_indexes[:,1] - total_e

	return true_pro_enh_indexes, prom_enh_false_interactions


def MOG_classifier(mode_of_sampler, number_of_samples, burn_in, comb = "ER", kappa_0=1.0, mu_0=1.0, alpha_0=1.0, Beta_0=1.0, total_posterior = False, pairwise_number_in_pack = 150, chain = 1):
	import os

	def name_creator_and_pair_wise_exist_check(chrom, mode_of_sampler, comb, number_of_samples, burn_in, chain):

		if mode_of_sampler == "distance_prior": name = 'prior_distance_trace_of_c_{0}_{1}'.format(chrom, number_of_samples)
		elif mode_of_sampler == "distance_MOG": name = 'MOG_distance_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples)
		elif mode_of_sampler == "dirichlet_MOG": name = 'MOG_dirichlet_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples)
		elif mode_of_sampler == "distance_MOG_empir_mu": name = 'MOG_distance_emprirical_mu_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples)

		if chain: name = name + "_{0}".format(chain); print chain

		output_folder = "./MOG_results_/"

		temp_output = output_folder + "Pairwise_prob/"

		if not os.path.exists(temp_output): os.makedirs(temp_output)

		def does_pair_wise_exists():

			pair_wise_prob = temp_output + name + "_{0}_pairwise_prob.npy".format(burn_in)

			pair_wise_prob_exist = os.path.exists(pair_wise_prob)
			
			return pair_wise_prob_exist, pair_wise_prob

		name_of_MOG_chain_file = output_folder + name

		print name

		pair_wise_prob_exist, name_of_MOG_chain_file_pair_wise_prob = does_pair_wise_exists()

		print name, pair_wise_prob_exist

		return name_of_MOG_chain_file, name_of_MOG_chain_file_pair_wise_prob, pair_wise_prob_exist 

	def loads_MoG_results(chrom, name):

		import iter_loadtxt

		_c_trace_raw = iter_loadtxt.iter_loadtxt(name, ",", dtype = int) # saves memory

		if mode_of_sampler == "dirichlet_MOG":

			_c_trace = _c_trace_raw

		else:

			num_of_promoters = len(config_variables.dict_chrom_pro_survived[chrom])	
			promoters_fixed_labels = np.zeros((len(_c_trace_raw), num_of_promoters),dtype = int)
			promoters_fixed_labels[:] = np.arange(num_of_promoters, dtype = int)
			_c_trace = np.c_[promoters_fixed_labels, _c_trace_raw]

		return _c_trace

	def cluster_estimator_similarity(_c_trace):

		from multiprocessing import Pool
		pool = Pool(processes = 6)

		#pack = number_of_samples
		#dim = _c_trace.shape
		#incr = int(dim[0]/pack)

		#number_of_samples = 

		bins = range(0, number_of_samples, pairwise_number_in_pack)
	
		if bins[-1] <> number_of_samples: bins.append(number_of_samples)	

		a = [_c_trace[bins[i]:bins[i+1]] for i in range(len(bins[:-1]))]

		import pararell_methods

		start = time.time()

		total_matrix = sum(pool.imap_unordered(pararell_methods.pararell_calc_ne, a))
		pool.close()
		pool.join()
		end = time.time()
		print end-start
		

		return total_matrix

	def standard_size_converter(total_matrix, chrom):

		indexes_p, indexes_e, total_p, total_e = chrom_specific_negative_interactions.initialise_variables(chrom)[2:]
		length_chr = len(indexes_p) + len(indexes_e)

		interaction_matrix = np.zeros((length_chr, length_chr), bool)

		chrom_pro_not_survived = config_variables.dict_chrom_pro_not_survived[chrom]
		chrom_enh_not_survived = config_variables.dict_chrom_enh_not_survived[chrom]
		dict_chrom_proximal = config_variables.dict_chrom_proximal

		if len(chrom_pro_not_survived): interaction_matrix[chrom_pro_not_survived - total_p, :] = True
		if len(chrom_enh_not_survived): interaction_matrix[:, len(indexes_p) + chrom_enh_not_survived - total_e] = True # gets rid of filtered out enhancers which could be causing nans due to their correlations 
		if config_variables.distant_enh_only and len(dict_chrom_proximal[chrom]): interaction_matrix[:, len(indexes_p) + dict_chrom_proximal[chrom] - total_e] = True

		interaction_matrix = np.invert(interaction_matrix + interaction_matrix.T)
		temp_expanded_total_matrix = np.zeros(length_chr * length_chr, int)

		temp_expanded_total_matrix[np.ravel(interaction_matrix)] = np.ravel(total_matrix)
		expanded_total_matrix = temp_expanded_total_matrix.reshape(length_chr, length_chr)

		return expanded_total_matrix

	posterior_of_option = {}
	chrom_posterior = {}

	for classification_of_interactions in ["positive_interactions", "negative_interactions"]: posterior_of_option[classification_of_interactions] = {}

	for chrom_ in chroms_to_infer:

		name_of_MOG_chain_file, name_of_MOG_chain_file_pair_wise_prob, pair_wise_prob_exist = name_creator_and_pair_wise_exist_check(chrom_, mode_of_sampler, comb, number_of_samples, burn_in, chain)

		if not pair_wise_prob_exist:

			_c_trace_distance = loads_MoG_results(chrom_, name_of_MOG_chain_file)

			_c_trace_distance = _c_trace_distance[burn_in:]

			total_matrix = cluster_estimator_similarity(_c_trace_distance)

			np.save(name_of_MOG_chain_file_pair_wise_prob, total_matrix)

		else:

			total_matrix = np.load(name_of_MOG_chain_file_pair_wise_prob)


		total_matrix = standard_size_converter(total_matrix, chrom_)

		total_matrix = total_matrix / float(number_of_samples - burn_in)

		print total_matrix

		true_pro_enh_indexes, prom_enh_false_interactions = interactions_extractor(chrom_)

		indexes_p, indexes_e, total_p, total_e = chrom_specific_negative_interactions.initialise_variables(chrom_)[2:]
		chrom_posterior["positive_interactions"] = total_matrix[true_pro_enh_indexes[:,0], true_pro_enh_indexes[:,1] + len(indexes_p)]
		chrom_posterior["negative_interactions"] = total_matrix[prom_enh_false_interactions[:,0], prom_enh_false_interactions[:,1] + len(indexes_p)]

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]: posterior_of_option[classification_of_interactions][chrom_] = chrom_posterior[classification_of_interactions]

	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		if total_posterior: posterior_of_option[classification_of_interactions] = list(itertools.chain.from_iterable([posterior_of_option[classification_of_interactions][chrom__] for chrom__ in chroms_to_infer]))

	return posterior_of_option["positive_interactions"], posterior_of_option["negative_interactions"]

def posterior_producer_non_domain(option_dist, option_correl, total_posterior = False):
	
	option_2 = option_dist + ( np.array(option_correl) + 1).tolist()
	print option_2
	posterior_of_option = {}

	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		posterior_of_option[classification_of_interactions] = {}

		for chrom_ in chroms_to_infer:
			probs_true = classifiers_elements[classification_of_interactions]["probabilities_of_being_positive_interactions"]["column_stacked_probabilities"][chrom_][:, option_2]
			probs_false = classifiers_elements[classification_of_interactions]["probabilities_of_being_negative_interactions"]["column_stacked_probabilities"][chrom_][:, option_2]
			K = len(dict_chrom_pro_survived[chrom_])
			chrom_posterior = classifier(probs_true, probs_false, K)	
			posterior_of_option[classification_of_interactions][chrom_] = chrom_posterior

		if total_posterior: posterior_of_option[classification_of_interactions] = list(itertools.chain.from_iterable([posterior_of_option[classification_of_interactions][chrom_] for chrom_ in chroms_to_infer]))

	return	posterior_of_option["positive_interactions"], posterior_of_option["negative_interactions"]


def posterior_producer_domain(option_dist, option_correl, total_posterior = False):
	import domain_allocator_clean

	option_2 = option_dist + ( np.array(option_correl) + 1).tolist()
	print option_2
	posterior_of_option = {}

	for classification_of_interactions in ["positive_interactions", "negative_interactions"]: posterior_of_option[classification_of_interactions] = {}

	one_totals = 0
	total_of_non_empty = 0
	for chrom_ in chroms_to_infer:
		allocations = {}
		allocations["positive_interactions"], allocations["negative_interactions"], filtered_promoters_in_sub_domains, size_domains = domain_allocator_clean.allocator(chrom_)
		
		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:

			probs_true = classifiers_elements[classification_of_interactions]["probabilities_of_being_positive_interactions"]["column_stacked_probabilities"][chrom_][:, option_2]
			probs_false = classifiers_elements[classification_of_interactions]["probabilities_of_being_negative_interactions"]["column_stacked_probabilities"][chrom_][:, option_2]

			chrom_posterior = np.zeros(len(probs_true)) # will store classifiers result for each positive or negative interaction respecively

			for sub_domain in np.unique(allocations[classification_of_interactions]):
				
				K = filtered_promoters_in_sub_domains[sub_domain]
				if K == 1: one_totals +=1
				if K: total_of_non_empty += 1
				interactions_of_subdomain = allocations[classification_of_interactions]  == sub_domain

				if sum(interactions_of_subdomain) == 0: continue
				print K
				'make sure that the oprobabilities of the interactions are in the right corresponding entries of the vector so that interaction <--> correct probability' 
				chrom_posterior[interactions_of_subdomain] = classifier(probs_true[interactions_of_subdomain], probs_false[interactions_of_subdomain], K)
			posterior_of_option[classification_of_interactions][chrom_] = chrom_posterior

	print "number of single promoter domains:", one_totals, total_of_non_empty, one_totals/float(total_of_non_empty)
	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		if total_posterior: posterior_of_option[classification_of_interactions] = list(itertools.chain.from_iterable([posterior_of_option[classification_of_interactions][chrom_] for chrom_ in chroms_to_infer]))

	return	posterior_of_option["positive_interactions"], posterior_of_option["negative_interactions"]


def classifier_alternative(probs_true, probs_false, mask_of_existing_interactions_of_enhancer):

	single_k_model_ratios = {}

	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		mask = mask_of_existing_interactions_of_enhancer[classification_of_interactions]
		single_k_model_ratios[classification_of_interactions] = (probs_true[classification_of_interactions][mask] / probs_false[classification_of_interactions][mask]).prod(1)

	denominator = np.sum(np.r_[single_k_model_ratios["positive_interactions"], single_k_model_ratios["negative_interactions"]])	

	return single_k_model_ratios["positive_interactions"]/denominator, single_k_model_ratios["negative_interactions"]/denominator


def posterior_producer_one_promoter_model(option_dist, option_correl, total_posterior = False):


	option_2 = option_dist + ( np.array(option_correl) + 1 ).tolist()
	print option_2
	posterior_of_option = {}


	for classification_of_interactions in ["positive_interactions", "negative_interactions"]: posterior_of_option[classification_of_interactions] = {}

	for chrom_ in chroms_to_infer:

		interactions = {}
		interactions["positive_interactions"], interactions["negative_interactions"] = interactions_extractor(chrom_)

		probs_true = {}
		probs_false = {}
		chrom_posterior = {}

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:

			probs_true[classification_of_interactions] = classifiers_elements[classification_of_interactions]["probabilities_of_being_positive_interactions"]["column_stacked_probabilities"][chrom_][:, option_2]
			probs_false[classification_of_interactions] = classifiers_elements[classification_of_interactions]["probabilities_of_being_negative_interactions"]["column_stacked_probabilities"][chrom_][:, option_2]
			chrom_posterior[classification_of_interactions] = np.zeros(len(interactions[classification_of_interactions]))

		mask_of_existing_interactions_of_enhancer = {}

		enhancers_to_infer = np.unique(np.r_[interactions["positive_interactions"][:,1], interactions["negative_interactions"][:,1]])

		for enhancer in enhancers_to_infer:
			for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
				mask_of_existing_interactions_of_enhancer[classification_of_interactions] = interactions[classification_of_interactions][:,1] == enhancer

			chrom_posterior["positive_interactions"][mask_of_existing_interactions_of_enhancer["positive_interactions"]], chrom_posterior["negative_interactions"][mask_of_existing_interactions_of_enhancer["negative_interactions"]] = classifier_alternative(probs_true, probs_false, mask_of_existing_interactions_of_enhancer)

		for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
			posterior_of_option[classification_of_interactions][chrom_] = chrom_posterior[classification_of_interactions]
	for classification_of_interactions in ["positive_interactions", "negative_interactions"]:
		if total_posterior: posterior_of_option[classification_of_interactions] = list(itertools.chain.from_iterable([posterior_of_option[classification_of_interactions][chrom_] for chrom_ in chroms_to_infer]))

	return	posterior_of_option["positive_interactions"], posterior_of_option["negative_interactions"]


if alternative_classificator: posterior_producer = posterior_producer_one_promoter_model
elif config_variables.domain: posterior_producer = posterior_producer_domain
else: posterior_producer = posterior_producer_non_domain




