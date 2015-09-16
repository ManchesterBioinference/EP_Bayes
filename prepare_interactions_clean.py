import config_variables
np = config_variables.np
re = config_variables.re

def un_string(array_to_clean):  return np.column_stack((map(lambda x: int(re.findall('\d+', x)[0]), array_to_clean[:,0]), map(lambda x: int(re.findall('\d+', x)[0]), array_to_clean[:,1])))

def filter_true_interactions_of_promoters_and_enhancers_which_didnt_survive_filtering(f_name):
	dataset_time_series_dict = config_variables.dataset_time_series_dict
	link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
	re = config_variables.re
	filtered_promoters = config_variables.filtered_promoters
	filtered_enhancers = config_variables.filtered_enhancers
	path = config_variables.path



	#--------------------------------------------------------------
	def interactions_chrom_splitter_dict(interactions): # performs chromosome split and filters out paths longer than threshold.
		chr_interactions_dict={}		
		for i in np.r_[np.array(range(1, 23), dtype='S2'), ['X'], ['Y']]: chr_interactions_dict['chr{0}'.format(i)] = interactions[(interactions[:,0] == 'chr{0}'.format(i))*(interactions[:,3].astype(int) <= path)][:, 1:]
		return chr_interactions_dict

	true_interactions_all = np.loadtxt(f_name, usecols = range(0,4), dtype = str)

	true_interactions_letters = np.array(map(lambda x: re.findall('\D+', x)[0], true_interactions_all[:,1]))

	global_pro_enh = true_interactions_all[true_interactions_letters == 'p']
	global_pro_enh_ind = un_string(global_pro_enh[:, 1:])
	
	global_enh_enh = true_interactions_all[true_interactions_letters == 'e']
	global_enh_enh_ind = un_string(global_enh_enh[:, 1:])	

	survived_extra_filtering_pro = np.in1d(global_pro_enh_ind[:, 0], np.where(filtered_promoters)[0])*np.in1d(global_pro_enh_ind[:, 1], np.where(filtered_enhancers)[0])
	survived_extra_filtering_enh = np.in1d(global_enh_enh_ind[:, 0], np.where(filtered_enhancers)[0])*np.in1d(global_enh_enh_ind[:, 1], np.where(filtered_enhancers)[0])

	chr_interactions_dict_enh_enh = interactions_chrom_splitter_dict(true_interactions_all[true_interactions_letters == 'e'][survived_extra_filtering_enh]) # applies data_specific filters to interactions
	chr_interactions_dict_pro_enh = interactions_chrom_splitter_dict(true_interactions_all[true_interactions_letters == 'p'][survived_extra_filtering_pro]) # applies data_specific filters to interactions

	
	def total_for_chromosomes():
	
		pro_chroms = dataset_time_series_dict[link_data_set_name_to_file_name["promoters"]['ER']][0] # comment: I should keep it "ER" independent in case someone deleted "ER"
		enh_chroms = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"]['ER']][0]
		total_e, total_p, dict_total_enh, dict_total_pro = 0, 0, {}, {} 

		for i in np.concatenate((np.array(range(1, 23), dtype='S2'), ['X'], ['Y'])): 
			dict_total_enh['chr{0}'.format(i)], dict_total_pro['chr{0}'.format(i)] = total_e, total_p 
			total_e += np.sum(enh_chroms=='chr{0}'.format(i)); total_p += np.sum(pro_chroms=='chr{0}'.format(i))

		return dict_total_enh, dict_total_pro

	dict_total_enh, dict_total_pro = total_for_chromosomes()

	return chr_interactions_dict_pro_enh, chr_interactions_dict_enh_enh, dict_total_enh, dict_total_pro
