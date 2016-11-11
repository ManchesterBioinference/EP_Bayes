
from prepare_interactions_clean import un_string

def inter_enhancer(chrom):
	negative_interactions = config_variables.negative_interactions
	indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

	if config_variables.disentagled_features_validation: 
		chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
	else:
		chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

	true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

	i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]
	interacting_enhancers = np.unique(j_s_t)-total_e
	return interacting_enhancers

suma = 0
for chrom in chrom_names:
	suma += len(inter_enhancer(chrom))

