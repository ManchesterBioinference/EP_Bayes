#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)

import os as os
cwd = os.getcwd()
path_to_R = cwd + "/R_scripts/"
os.chdir(path_to_R)
os.system("tar xvzf data_temp_output_for_cluster_figures.tar.gz")
os.chdir(cwd)

figure_index = 0
for figure_index in range(5):
	cluster_figure_selection = ["cluster_ER_enhancer", "cluster_Pol2s_enhancer", "cluster_Pol2s_promoter", "cluster_ER_promoter", "cluster_Pol2s_ER_enhancer"][figure_index]
	ENHANCER_GENE_INTERACTOR.executor(do_clustering = True, redo_raw_CHIA_PET_interactions = False, cluster_figure_selection = cluster_figure_selection)

