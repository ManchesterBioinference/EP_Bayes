#------------------------------------------------------------------------------------------------

#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)


PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = True, upstream = 1500, downstream = 0, upstream_t_s = 300)

#The code will produce figures 2, 3b, 3d, 3f 

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][0]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 0, upstream_t_s = 300)

#The code will produce figure 3a, 4a, 12s

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][2]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 0, upstream_t_s = 300)



#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)


PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][1]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = True, upstream = 1500, downstream = 1500, upstream_t_s = 300)

#The code will produce figures 2, 3b, 3d, 3f 

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][0]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][1]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 1500, upstream_t_s = 300)

#The code will produce figure 3a, 4a, 12s

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][2]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][1]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 1500, upstream_t_s = 300)



#ENRICHMENTS:

ENHANCER_GENE_INTERACTOR.executor(rising_of_falling_POL2_or_ER = ["ER", "Pol2"][1], plot_TF_enrichments_in_cluster = True, DB_version = False, redo_raw_CHIA_PET_interactions = False)

ENHANCER_GENE_INTERACTOR.executor(rising_of_falling_POL2_or_ER = ["ER", "Pol2"][1], plot_TF_enrichments_in_cluster = True, DB_version = True, redo_raw_CHIA_PET_interactions = False)


#------------------------------------------------------------------------------------------------

#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)

#"ALL_option"

PR_CURVES = ["SELECTIVE", "ALL"][1]
mode_of_code = ["ODD", "EVEN", "FULL"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = True, upstream = 1500, downstream = 0, upstream_t_s = 300)

#The code will produce figures 2, 3b, 3d, 3f 

PR_CURVES = ["SELECTIVE", "ALL"][1]
mode_of_code = ["ODD", "EVEN", "FULL"][0]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 0, upstream_t_s = 300)

#The code will produce figure 3a, 4a, 12s

PR_CURVES = ["SELECTIVE", "ALL"][1]
mode_of_code = ["ODD", "EVEN", "FULL"][2]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][0]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 0, upstream_t_s = 300)

#"ALL_option"

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][1]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = True, upstream = 1500, downstream = 1500, upstream_t_s = 300)

#The code will produce figures 2, 3b, 3d, 3f 

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][0]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][1]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 1500, upstream_t_s = 300)

#The code will produce figure 3a, 4a, 12s

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL"][2]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][1]
GENE_OR_PROMOTER_MODE = ["GENE_MODE", "TSS_MODE"][1]

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, GENE_OR_PROMOTER_MODE = GENE_OR_PROMOTER_MODE, redo_raw_CHIA_PET_interactions = False, upstream = 1500, downstream = 1500, upstream_t_s = 300)





