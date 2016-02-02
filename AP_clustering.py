def concatenator(cluster_mode, merge_time_series_option, count_filter_each_data_set, pol2_rep_correl_filt = False, distant_enh_only = False):

	import config_variables
	dataset_time_series_dict = config_variables.dataset_time_series_dict
	dataset_time_series_dict_mean_std = config_variables.dataset_time_series_dict_mean_std
	time_points = config_variables.time_points
	name_of_time_series_file = config_variables.name_of_time_series_file
	np = config_variables.np
	link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name

	name = name_of_time_series_file[cluster_mode]
	name_of_overlap_file_dict = config_variables.name_of_overlap_file_dict

	merge_time_series_option_ = [name + '_' + el for el in merge_time_series_option]
	total_mask = np.ones(len(dataset_time_series_dict[link_data_set_name_to_file_name[cluster_mode]['ER']][0]), bool)

	if type(pol2_rep_correl_filt) == float:
		rep_1_ts, rep_2_ts = dataset_time_series_dict[link_data_set_name_to_file_name[cluster_mode]['PolII_2012-03']][2][:, :time_points], dataset_time_series_dict[link_data_set_name_to_file_name[cluster_mode]['PolII']][2][:, :time_points]
		correlations_1 = np.zeros(len(rep_1_ts))

		for ind, el in enumerate(np.c_[rep_1_ts, rep_2_ts]): correlations_1[ind] = np.corrcoef(el[:time_points], el[time_points:])[0,1]

		total_mask = total_mask*(correlations_1 > pol2_rep_correl_filt)	

	for opt in merge_time_series_option_: total_mask = total_mask*(dataset_time_series_dict[opt][2].sum(1) >= count_filter_each_data_set)

	
	if distant_enh_only:
		overlapping_peaks = np.loadtxt(name_of_overlap_file_dict[cluster_mode], dtype = str)[:,-1].astype(int)
		total_mask[overlapping_peaks] = False

	concat_set = dataset_time_series_dict_mean_std[link_data_set_name_to_file_name[cluster_mode][merge_time_series_option[0]]][2][total_mask]

	for opt in merge_time_series_option[1:]: concat_set = np.c_[concat_set, dataset_time_series_dict_mean_std[link_data_set_name_to_file_name[cluster_mode][opt]][2][total_mask]]

	name_of_merged_time_series_to_cluster = name + '_concat'
	for el in merge_time_series_option: name_of_merged_time_series_to_cluster = name_of_merged_time_series_to_cluster + '_' + el
	name_of_merged_time_series_to_cluster = "{0}_{1}".format(name_of_merged_time_series_to_cluster, count_filter_each_data_set)

	if type(pol2_rep_correl_filt) == float: name_of_merged_time_series_to_cluster = name_of_merged_time_series_to_cluster + "_cor_{0}".format(pol2_rep_correl_filt)
	if distant_enh_only: name_of_merged_time_series_to_cluster = name_of_merged_time_series_to_cluster + "_distant_only"

	survived_indexes = np.where(total_mask)[0]
	import os as os
	path_to_R = os.getcwd() + "/R_scripts/"
	np.savetxt(path_to_R + name_of_merged_time_series_to_cluster + "_survived_indexes", survived_indexes)
	np.savetxt(path_to_R + name_of_merged_time_series_to_cluster, concat_set, delimiter=',')

	return name_of_merged_time_series_to_cluster


def AP_clustering(name_of_file_clustered, number_of_clusters):
	
	import os as os  #or popen
	cwd = os.getcwd()
	path_to_R = os.getcwd() + "/R_scripts/"
	os.chdir(path_to_R)
	os.system("Rscript AP_clustering_executor.R {0} {1}".format(name_of_file_clustered, number_of_clusters))
	os.chdir(cwd)


