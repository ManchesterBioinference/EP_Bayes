#------------------------------------------------------------------------------------------------

#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)


PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = True, upstream = 300)

#The code will produce figures 2, 3b, 3d, 3f 

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][0]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = False, upstream = 300)


#The code will produce figure 3a, 4a, 12s

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][2]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = False, upstream = 300)


#the code will produce figure 1a,b,c,d and 1s 

#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)

#figure_index = 0
for figure_index in range(5):
	cluster_figure_selection = ["cluster_ER_enhancer", "cluster_Pol2s_enhancer", "cluster_Pol2s_promoter", "cluster_ER_promoter", "cluster_Pol2s_ER_enhancer"][figure_index]
	ENHANCER_GENE_INTERACTOR.executor(do_clustering = True, redo_raw_CHIA_PET_interactions = False, cluster_figure_selection = cluster_figure_selection)


# redo clustering
#ENHANCER_GENE_INTERACTOR.executor(do_clustering = True, re_do_clustering = True, redo_raw_CHIA_PET_interactions = False, cluster_figure_selection = "cluster_Pol2s_ER_enhancer_test")


#------------------------------------------------------------------------------------------------------------------------------

#FULL figures:

#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)


PR_CURVES = ["SELECTIVE", "ALL"][1]
mode_of_code = ["ODD", "EVEN", "FULL"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = True, upstream = 300)

#The code will produce figures 2, 3b, 3d, 3f 

PR_CURVES = ["SELECTIVE", "ALL"][1]
mode_of_code = ["ODD", "EVEN", "FULL"][0]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = False, upstream = 300)


#The code will produce figure 3a, 4a, 12s

PR_CURVES = ["SELECTIVE", "ALL"][1]
mode_of_code = ["ODD", "EVEN", "FULL"][2]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = False, upstream = 300)

