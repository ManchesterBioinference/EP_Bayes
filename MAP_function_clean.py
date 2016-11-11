import config_variables 
from prepare_interactions_clean import un_string

np = config_variables.np
mode = config_variables.mode

negative_interactions = config_variables.negative_interactions
mode = config_variables.mode
domain = config_variables.domain
alternative_classificator = config_variables.alternative_classificator

def positive_negative_interactions_for_MAP(chrom):

	indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

	if mode == "promoter_enhancer_interactions":

		false_inter_pro = negative_interactions.chrom_specific_negative_interactions(chrom, mode)

		i_s_f, j_s_f = false_inter_pro[:,0] + total_p, false_inter_pro[:,1] + total_e

	if mode == "enhancer_enhancer_interactions":

		false_inter_enh = negative_interactions.chrom_specific_negative_interactions(chrom, mode)
	
		i_s_f, j_s_f = false_inter_enh[:,0] + total_e, false_inter_enh[:,1] + total_e

	return i_s_f, j_s_f

def MAP(posterior_t, posterior_f, chrom, i_s_f, j_s_f):

	indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]
		
	length_chr = len(indexes_p) + len(indexes_e)
	interaction_matrix = np.zeros((length_chr, length_chr))
	posterior_t, posterior_f = posterior_t[chrom], posterior_f[chrom]


	if mode == "promoter_enhancer_interactions":

		if config_variables.disentagled_features_validation: 
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
		else:
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

		true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

		i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]


		interaction_matrix[:,:] = np.min([np.min(posterior_t), np.min(posterior_f)])*0.99
		interaction_matrix[i_s_t - total_p, j_s_t + len(indexes_p) - total_e] = posterior_t
		interaction_matrix[i_s_f - total_p, j_s_f + len(indexes_p) - total_e] = posterior_f

		MAP_indexes = np.argmax(interaction_matrix, axis = 0)

		if config_variables.alternative_classificator_outside_enhancers: 
			chrom_interacting_enhancers_pro = config_variables.chrom_interacting_enhancers_pro[chrom]
			max_poster_is_pro = MAP_indexes[len(indexes_p) + chrom_interacting_enhancers_pro - total_e] + total_p	
			MAP_predicted_intereactions_pro = np.column_stack((max_poster_is_pro, chrom_interacting_enhancers_pro))

		else: 
			max_poster_is_pro = MAP_indexes[len(indexes_p) + np.unique(j_s_t) - total_e] + total_p # gives a maximum aposteriori promoter to each infered enhancer
			MAP_predicted_intereactions_pro = np.column_stack((max_poster_is_pro, np.unique(j_s_t)))

		infered_promoters_pro_enh = MAP_indexes[len(indexes_p):] + total_p
		MAP_probabilites_pro_enh = interaction_matrix[MAP_indexes, range(len(indexes_p) + len(indexes_e))][len(indexes_p):]

		probabilities_for_promoters_of_interacting_enhancers = interaction_matrix[:len(indexes_p), len(indexes_p) + np.unique(j_s_t) - total_e]

		link_exists_pro = [ind for ind, el in enumerate(MAP_predicted_intereactions_pro.tolist()) if el in true_inter_pro.tolist()]
		mask_link_exists_pro = np.zeros(len(MAP_predicted_intereactions_pro), bool)
		mask_link_exists_pro[link_exists_pro] = True
	
		return mask_link_exists_pro, infered_promoters_pro_enh, MAP_probabilites_pro_enh, probabilities_for_promoters_of_interacting_enhancers

	if mode == "enhancer_enhancer_interactions":

		chr_interactions_dict_enh_enh = config_variables.chr_interactions_dict_enh_enh
		true_inter_enh = un_string(chr_interactions_dict_enh_enh[chrom][:, :2]).astype(int)
		i_s_t, j_s_t = true_inter_enh[:,0], true_inter_enh[:,1]

		interaction_matrix[:,:] = np.min([np.min(posterior_t), np.min(posterior_f)])*0.99
		interaction_matrix[i_s_t + len(indexes_p) - total_e, j_s_t + len(indexes_p) - total_e] = posterior_t
		interaction_matrix[i_s_f + len(indexes_p) - total_e, j_s_f + len(indexes_p) - total_e] = posterior_f
		interaction_matrix[j_s_t + len(indexes_p) - total_e, i_s_t + len(indexes_p) - total_e] = posterior_t # transpose to create a full matrix
		interaction_matrix[j_s_f + len(indexes_p) - total_e, i_s_f + len(indexes_p) - total_e] = posterior_f # transpose to create a full matrix

		MAP_indexes = np.argmax(interaction_matrix, axis = 0)
		max_poster_is_enh = MAP_indexes[len(indexes_p) + np.unique(j_s_t) - total_e] - len(indexes_p) + total_e

		infered_enhancers_enh_enh = MAP_indexes[len(indexes_p):] + total_e 

		MAP_probabilites_enh_enh = interaction_matrix[MAP_indexes, range(len(indexes_p) + len(indexes_e))][len(indexes_p):]

		MAP_predicted_intereactions_enh = np.column_stack((max_poster_is_enh, np.unique(j_s_t)))

		link_exists_enh = [ind for ind, el in enumerate(MAP_predicted_intereactions_enh.tolist()) if el in true_inter_enh.tolist()]
		mask_link_exists_enh = np.zeros(len(MAP_predicted_intereactions_enh), bool)
		mask_link_exists_enh[link_exists_enh] = True

		return mask_link_exists_enh, infered_enhancers_enh_enh, MAP_probabilites_enh_enh
	

