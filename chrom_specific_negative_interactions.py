import config_variables
from prepare_interactions_clean import un_string

link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
dataset_time_series_dict = config_variables.dataset_time_series_dict
np = config_variables.np
dict_total_pro = config_variables.dict_total_pro
dict_total_enh = config_variables.dict_total_enh
dict_chrom_pro_not_survived = config_variables.dict_chrom_pro_not_survived
dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived
dict_chrom_enh_not_survived = config_variables.dict_chrom_enh_not_survived
dict_chrom_enh_survived = config_variables.dict_chrom_enh_survived
upstream = config_variables.upstream
domain = config_variables.domain
chr_interactions_dict_pro_enh = config_variables.chr_interactions_dict_pro_enh
chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh
interacting_negatives = config_variables.interacting_negatives
TSS_or_intra_genic_for_domain_filter = config_variables.TSS_or_intra_genic_for_domain_filter
one_sided_or_two_sided = config_variables.one_sided_or_two_sided
interacting_enhancers_only = config_variables.interacting_enhancers_only
dict_chrom_proximal = config_variables.dict_chrom_proximal
distant_enh_only = config_variables.distant_enh_only
name_of_time_series_promoter_file_for_TSS_start = config_variables.name_of_time_series_promoter_file_for_TSS_start

if domain: import interacting_domain

def calculate_correlations(pro_ts, enh_ts):

	time_series = np.concatenate((pro_ts, enh_ts), axis = 0)
	correlations = np.corrcoef(time_series)
	return correlations

def calculate_distances(domain, points_p, points_e):

	promoter_enhancers_point_coordinates = np.r_[points_p, points_e]

	distances_matrix = promoter_enhancers_point_coordinates[:, None] - promoter_enhancers_point_coordinates
	return distances_matrix

def extract_TSS_coordinates(upstream):

	data = np.loadtxt(config_variables.name_of_time_series_promoter_file_for_TSS_start, dtype = str,  delimiter = '\t')	
	plus_strand = data[:, 4] == '+'
	TSS_coordinates = np.zeros(len(plus_strand), int)
	TSS_coordinates[plus_strand] = data[plus_strand, 1].astype(int) - upstream
	TSS_coordinates[np.invert(plus_strand)] = data[np.invert(plus_strand), 2].astype(int) + upstream

	return TSS_coordinates



TSS_coordinates = extract_TSS_coordinates(upstream)

def initialise_variables(chrom):

	name_of_pro_t_s = link_data_set_name_to_file_name["promoters"]["ER"]
	name_of_enh_t_s = link_data_set_name_to_file_name["enhancers"]["ER"]

	pro_chroms, pro_coordinates, ts_ = dataset_time_series_dict[name_of_pro_t_s]
	enh_chroms, enh_coordinates, ts = dataset_time_series_dict[name_of_enh_t_s]

	indexes_p = np.where(pro_chroms==chrom)[0]	# gives the number of promoters for a chromosome	
	indexes_e = np.where(enh_chroms==chrom)[0]	# gives the number of enhancers for a chromosome

	total_p = dict_total_pro[chrom]
	total_e = dict_total_enh[chrom]

	return enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e 

def chrom_specific_negative_interactions(chrom, mode, prior_mode = False):

	enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e = initialise_variables(chrom)

	length_chr = len(indexes_p) + len(indexes_e)

	chrom_pro_not_survived = dict_chrom_pro_not_survived[chrom]
	chrom_pro_survived = dict_chrom_pro_survived[chrom]
	chrom_enh_not_survived = dict_chrom_enh_not_survived[chrom]
	chrom_enh_survived = dict_chrom_enh_survived[chrom]	

	
	
	# interaction_domains_adjustments beginning -----------------------
	if domain:
		if TSS_or_intra_genic_for_domain_filter == "Intra_genic": coords_pro_domain = pro_coordinates[indexes_p]
		elif TSS_or_intra_genic_for_domain_filter == "TSS_only": coords_pro_domain = np.column_stack((TSS_coordinates[indexes_p]-1, TSS_coordinates[indexes_p]+1))
		domain_matrix = interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'left', True)
		domain_matrix = domain_matrix + interacting_domain.interacting_domains(coords_pro_domain, enh_coordinates[indexes_e], chrom, 'right', True)
	else:
		domain_matrix = True
	# interaction_domains_adjustments_ending -----------------------

	def promoter_enhancer_interactions_generator():

	
		chr_interactions_pro_enh = chr_interactions_dict_pro_enh[chrom]

		if config_variables.alternative_classificator_outside_enhancers: chrom_interacting_enhancers_pro = config_variables.chrom_interacting_enhancers_pro[chrom]
		else: chrom_interacting_enhancers_pro = np.unique(un_string(chr_interactions_dict_pro_enh[chrom])[:,1])
		chrom_interacting_promoters_pro = np.unique(un_string(chr_interactions_dict_pro_enh[chrom])[:,0])

		interaction_matrix = np.zeros((length_chr, length_chr), bool)
	
		interaction_matrix[range(length_chr), range(length_chr)] = True	# gets rid of diagonal
		interaction_matrix[np.tril_indices(length_chr)] = True # gets rid of symmetric interactions
	
		interaction_matrix[0:len(indexes_p), 0:len(indexes_p)] = True # gets rid of promoter_promoter_interactions

		features = np.array(['p{0}'.format(ind) for ind in indexes_p] + ['e{0}'.format(ind) for ind in indexes_e]) # creates a frame with chromosome specific interactions

		true_pro_enh_indexes = un_string(chr_interactions_pro_enh)

		print 'number of pro_enh true interactions: ', len(chr_interactions_pro_enh)

		if len(chrom_pro_not_survived): interaction_matrix[chrom_pro_not_survived - total_p, :] = True # gets rid of negative interactions which could be generated by filtered promoters

		if interacting_negatives:
			mask_interacting_promoters = np.zeros(length_chr).astype(bool)# we don't have to filter out enhancers which didn't pass the filter thresold. Since we consider only the interacting enhancers that's a subset of survived enhnacers.
			mask_interacting_promoters[chrom_interacting_promoters_pro - total_p] = True
			mask_non_interacting_promoters = np.invert(mask_interacting_promoters)	
			interaction_matrix[mask_non_interacting_promoters, len(indexes_p):] = True # it's equivalent to interacting_enhancers_mask_invert

		#if config_variables.disentagled_features_validation: 
			#true_pro_enh_indexes = un_string(config_variables.chr_interactions_dict_pro_enh_TSS[chrom])
			#chrom_interacting_enhancers_pro = np.unique(true_pro_enh_indexes[:, 1])

		mask_interacting_enhancers = np.zeros(length_chr).astype(bool)# we don't have to filter out enhancers which didn't pass the filter thresold. Since we consider only the interacting enhancers that's a subset of survived enhnacers.
		mask_interacting_enhancers[chrom_interacting_enhancers_pro - total_e + len(indexes_p)] = True
		mask_non_interacting_enhancers = np.invert(mask_interacting_enhancers)	
		#interaction_matrix[:len(indexes_p), mask_non_interacting_enhancers] = True # it's equivalent to interacting_enhancers_mask_invert

		if interacting_enhancers_only or prior_mode: interaction_matrix[:len(indexes_p), mask_non_interacting_enhancers] = True # it's equivalent to interacting_enhancers_mask_invert
		elif len(chrom_enh_not_survived): interaction_matrix[:len(indexes_p), len(indexes_p) + chrom_enh_not_survived - total_e] = True # gets rid of filtered out enhancers which could be causing nans due to their correlations 
		if distant_enh_only and len(dict_chrom_proximal[chrom]): interaction_matrix[:len(indexes_p), len(indexes_p) + dict_chrom_proximal[chrom] - total_e] = True
	
		interaction_matrix[true_pro_enh_indexes[:, 0] - total_p, true_pro_enh_indexes[:, 1] - total_e + len(indexes_p)] = True

		interaction_matrix[len(indexes_p): len(indexes_p) + len(indexes_e), len(indexes_p): len(indexes_p) + len(indexes_e)] = True # gets rid of enhancers-enhancer block

		
		indexes_of_zero_interactions = np.where(True == np.invert(interaction_matrix)*domain_matrix)
		column_1st = indexes_of_zero_interactions[0]
		column_2nd = indexes_of_zero_interactions[1] - len(indexes_p)

		prom_enh_false_interactions = np.concatenate((column_1st[:,None], column_2nd[:,None]), axis=1)
		#pro-enh interactions end-----------------------------------------------------------------
		return prom_enh_false_interactions


	if "promoter_enhancer_interactions" in mode: 
		negative_of_type_of_interactions = promoter_enhancer_interactions_generator()
		print 'number of pro_enh false interactions: ', len(negative_of_type_of_interactions)

	def enhancer_enhancer_interactions_generator():

		chr_interactions_enh_enh = chr_interactions_dict_enh_enh[chrom]
		chrom_interacting_enhancers_enh = np.unique(un_string(chr_interactions_dict_enh_enh[chrom])[:,0])
		chrom_interacting_enhancers_enh = np.unique(np.r_[chrom_interacting_enhancers_enh, np.unique(un_string(chr_interactions_dict_enh_enh[chrom])[:,1])])

		interaction_matrix = np.zeros((length_chr, length_chr), bool)
		interaction_matrix[range(length_chr), range(length_chr)] = True	# gets rid of diagonal
		interaction_matrix[0:len(indexes_p), 0:len(indexes_p)] = True # gets rid of promoter_promoter_interactions
		interaction_matrix[:len(indexes_p), len(indexes_p): len(indexes_p) + len(indexes_e)] = True # gets rid of promoter-enhancer block

		print 'number of enh_enh true interactions: ', len(chr_interactions_enh_enh)
		#enh-enh interactions start-----------------------------------------------------------------

	
		if len(chrom_enh_not_survived): interaction_matrix[len(indexes_p) + chrom_enh_not_survived - total_e, :] = True # sorts out raws
		if distant_enh_only and len(dict_chrom_proximal[chrom]): interaction_matrix[len(indexes_p) + dict_chrom_proximal[chrom] - total_e, :] = True

		if interacting_negatives:
		
			mask_interacting_enhancers = np.zeros(length_chr).astype(bool)
			mask_interacting_enhancers[chrom_interacting_enhancers_enh - total_e + len(indexes_p)] = True
			mask_non_interacting_enhancers = np.invert(mask_interacting_enhancers)	
			interaction_matrix[mask_non_interacting_enhancers, len(indexes_p):] = True


		#sort out columns--------------------------------------
		mask_interacting_enhancers = np.zeros(length_chr).astype(bool)
		mask_interacting_enhancers[chrom_interacting_enhancers_enh - total_e + len(indexes_p)] = True
		mask_non_interacting_enhancers = np.invert(mask_interacting_enhancers)	

		if interacting_enhancers_only or prior_mode: interaction_matrix[len(indexes_p):, mask_non_interacting_enhancers] = True # it's equivalent to interacting_enhancers_mask_invert
		elif len(chrom_enh_not_survived): interaction_matrix[len(indexes_p):, len(indexes_p) + chrom_enh_not_survived - total_e] = True # gets rid of filtered out enhancers which could be causing nans due to their correlations 
		if distant_enh_only and len(dict_chrom_proximal[chrom]): interaction_matrix[len(indexes_p):, len(indexes_p) + dict_chrom_proximal[chrom] - total_e] = True

			
		#sort out columns--------------------------------------end

		true_enh_enh_indexes = un_string(chr_interactions_enh_enh)
		interaction_matrix[true_enh_enh_indexes[:,0] - total_e + len(indexes_p), true_enh_enh_indexes[:,1] - total_e + len(indexes_p)] = True
		interaction_matrix[true_enh_enh_indexes[:,1] - total_e + len(indexes_p), true_enh_enh_indexes[:,0] - total_e + len(indexes_p)] = True
		interaction_matrix[np.tril_indices(length_chr)] = True # gets rid of symmetric interactions
	
		indexes_of_zero_interactions = np.where(True == np.invert(interaction_matrix)*domain_matrix)
		column_1st = indexes_of_zero_interactions[0] - len(indexes_p)
		column_2nd = indexes_of_zero_interactions[1] - len(indexes_p)

		enh_enh_false_interactions = np.concatenate((column_1st[:, None], column_2nd[:, None]), axis=1)
				

		return enh_enh_false_interactions

	if "enhancer_enhancer_interactions" in mode: 
		negative_of_type_of_interactions = enhancer_enhancer_interactions_generator()
		print 'number of enh_enh false interactions: ', len(negative_of_type_of_interactions)


	return negative_of_type_of_interactions
	#------------------------------------------------------------------------------------------------------------------


def distance_and_correl_of_interactions(chrom, negative_of_type_of_interactions, mode, data_set_name = "default", attribute_distance = False, attribute_correlation = False):

	print chrom, data_set_name

	enh_coordinates, pro_coordinates, indexes_p, indexes_e, total_p, total_e = initialise_variables(chrom)
	#--------------------------------
	#initialises variables

	#distance & correlation
	def uniqueness_of_enh_enh_atributes():
		mask_feature_uniqueness_enh_enh = np.zeros((len(indexes_e), len(indexes_e)), bool)
		mask_feature_uniqueness_enh_enh[negative_of_type_of_interactions[:,0], negative_of_type_of_interactions[:,1]] = True
		mask_feature_uniqueness_enh_enh = mask_feature_uniqueness_enh_enh + mask_feature_uniqueness_enh_enh.T
		mask_feature_uniqueness_enh_enh[np.tril_indices(len(mask_feature_uniqueness_enh_enh))] = False

		return mask_feature_uniqueness_enh_enh

	if ("enhancer_enhancer_interactions" in mode) * (attribute_distance + attribute_correlation): mask_feature_uniqueness_enh_enh = uniqueness_of_enh_enh_atributes()

	def distance_calculator():

		point_coordinates_promoter, point_coordinates_enhancer = TSS_coordinates[indexes_p], np.mean(enh_coordinates[indexes_e], axis = 1)
		distances_matrix = calculate_distances(domain, point_coordinates_promoter, point_coordinates_enhancer)


		if "promoter_enhancer_interactions" in mode:

			if config_variables.disentagled_features_validation: 
				chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
			else:
				chr_interactions_pro_enh = chr_interactions_dict_pro_enh[chrom]

			pro_enh_indexes = un_string(chr_interactions_pro_enh[:, :2])
			dist_from_enhancer_to_promoters = distances_matrix[len(indexes_p):, 0: len(indexes_p)]
			
			#if domain:
			#	random_negative_interaction_negative_distances = negative_of_type_of_interactions[np.random.choice(len(negative_of_type_of_interactions), int(len(negative_of_type_of_interactions)/2.), replace = False)]
			#	dist_from_enhancer_to_promoters[random_negative_interaction_negative_distances[:,1], random_negative_interaction_negative_distances[:,0]] *= -1

			dist_of_true_pro_enh_inter = dist_from_enhancer_to_promoters[pro_enh_indexes[:, 1]  - total_e, pro_enh_indexes[:, 0]  - total_p]
			dist_of_false_pro_enh_inter = dist_from_enhancer_to_promoters[negative_of_type_of_interactions[:, 1], negative_of_type_of_interactions[:, 0]]
			
			return dist_from_enhancer_to_promoters, dist_of_true_pro_enh_inter, dist_of_false_pro_enh_inter


		if "enhancer_enhancer_interactions" in mode:
			chr_interactions_enh_enh = chr_interactions_dict_enh_enh[chrom]
			enh_enh_indexes = un_string(chr_interactions_enh_enh[:, :2])

			dist_from_enhancer_to_enhancers = distances_matrix[len(indexes_p):, len(indexes_p):]

			if domain:
				random_interaction_negative_distances = np.ones(len(indexes_e), bool)
				random_interaction_negative_distances = np.c_[random_interaction_negative_distances[0][:, None], random_interaction_negative_distances[1][:, None]]
				random_interaction_negative_distances = random_interaction_negative_distances[np.random.choice(len(random_interaction_negative_distances), int(len(random_interaction_negative_distances)/2.), replace = False)]
				dist_from_enhancer_to_enhancers[random_interaction_negative_distances[:, 0], random_interaction_negative_distances[:, 1]] *= -1
				dist_from_enhancer_to_enhancers[random_interaction_negative_distances[:, 1], random_interaction_negative_distances[:, 0]] *= -1

			dist_of_false_enh_enh_inter = dist_from_enhancer_to_enhancers[mask_feature_uniqueness_enh_enh]	
			dist_of_true_enh_enh_inter = dist_from_enhancer_to_enhancers[enh_enh_indexes[:, 0] - total_e, enh_enh_indexes[:, 1] - total_e]

			return dist_from_enhancer_to_enhancers, dist_of_true_enh_enh_inter, dist_of_false_enh_enh_inter

	if attribute_distance: dist_from_type_of_interaction, dist_of_true_type_of_interaction, dist_of_false_type_of_interaction = distance_calculator()

	def correlation_calculator(data_set_name):

		name_of_pro_t_s = link_data_set_name_to_file_name["promoters"][data_set_name]
		name_of_enh_t_s = link_data_set_name_to_file_name["enhancers"][data_set_name]

		pro_time_series = dataset_time_series_dict[name_of_pro_t_s][2]
		enh_time_series = dataset_time_series_dict[name_of_enh_t_s][2]


		print 'correlator: start'
		pro_ts = pro_time_series[indexes_p]
		enh_ts = enh_time_series[indexes_e]
		chrom_correlations_matrix = calculate_correlations(pro_ts, enh_ts)
		print 'correlator: stop'

		if "promoter_enhancer_interactions" in mode:

			if config_variables.disentagled_features_validation: 
				chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
			else:
				chr_interactions_pro_enh = chr_interactions_dict_pro_enh[chrom]

			pro_enh_indexes = un_string(chr_interactions_pro_enh[:, :2])

			correl_from_enhancer_to_promoters = chrom_correlations_matrix[len(indexes_p):len(indexes_p) + len(indexes_e), 0: len(indexes_p)]
			correl_of_true_pro_enh_inter = correl_from_enhancer_to_promoters[pro_enh_indexes[:, 1] - total_e, pro_enh_indexes[:, 0] - total_p]
			correl_of_false_pro_enh_inter = correl_from_enhancer_to_promoters[negative_of_type_of_interactions[:, 1], negative_of_type_of_interactions[:, 0]]	

			print 'nan_correl', (np.isnan(correl_of_false_pro_enh_inter)).any(), (np.isnan(correl_of_true_pro_enh_inter)).any()

			return correl_from_enhancer_to_promoters, correl_of_true_pro_enh_inter, correl_of_false_pro_enh_inter

		if "enhancer_enhancer_interactions" in mode:
			chr_interactions_enh_enh = chr_interactions_dict_enh_enh[chrom]
			enh_enh_indexes = un_string(chr_interactions_enh_enh[:, :2])
	
			correl_from_enhancer_to_enhancers = chrom_correlations_matrix[len(indexes_p): len(indexes_p) + len(indexes_e), len(indexes_p): len(indexes_p) + len(indexes_e)]
			correl_of_true_enh_enh_inter = correl_from_enhancer_to_enhancers[enh_enh_indexes[:, 0] - total_e, enh_enh_indexes[:, 1] - total_e]
			correl_of_false_enh_enh_inter = correl_from_enhancer_to_enhancers[mask_feature_uniqueness_enh_enh]
		
			print 'nan_correl', (np.isnan(correl_of_false_enh_enh_inter)).any(), (np.isnan(correl_of_true_enh_enh_inter)).any()
	
			return correl_from_enhancer_to_enhancers, correl_of_true_enh_enh_inter, correl_of_false_enh_enh_inter

	if attribute_correlation: correl_from_type_of_interaction, correl_of_true_type_of_interaction, correl_of_false_type_of_interaction = correlation_calculator(data_set_name)

	
	if attribute_correlation: return correl_of_true_type_of_interaction, correl_of_false_type_of_interaction
	elif attribute_distance: return dist_of_true_type_of_interaction, dist_of_false_type_of_interaction
		
	




