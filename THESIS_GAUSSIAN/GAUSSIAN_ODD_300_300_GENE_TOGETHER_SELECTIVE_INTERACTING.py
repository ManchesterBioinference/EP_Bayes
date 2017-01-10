#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)

PR_CURVES = ["SELECTIVE", "ALL"][0]
mode_of_code = ["ODD", "EVEN", "FULL", "GAUSSIAN_SAMPLE", "MK_PAIRWISE", "convergence_checker"][0]
mode_of_code_2 = ["WITHOUT", "ADD_GAUSSIAN_VALIDATION"][1]
mode_of_features_and_interactions = ["FEATURES_AND_INTERACTIONS_TOGETHER", "FEATURES_AND_INTERACTIONS_SEPERATE"][0]


number_of_samples_correl, burn_in_correl = [80000, 80000, 80000, 80000, 160000], [40000, 40000, 40000, 40000, 80000] #[70000]*5, [10000]*5
number_of_samples_dist, burn_in_dist = 70000, 10000

kappa_0, mu_0, alpha_0, Beta_0 = 1.0, 0.0, 2.0, 2.0

mode_of_sampler = ["distance_prior", "distance_MOG", "dirichlet_MOG", "distance_MOG_empir_mu"][3]

interacting_enhancers_only_MOG = True

ENHANCER_GENE_INTERACTOR.executor(PR_CURVES, mode_of_code, mode_of_features_and_interactions, redo_raw_CHIA_PET_interactions = False, upstream = 300, option_correl_select = option_correl_select, kappa_0 = kappa_0, mu_0 = mu_0, alpha_0 = alpha_0, Beta_0 = Beta_0, mode_of_sampler = mode_of_sampler, csf_mode = True, interacting_enhancers_only_MOG = interacting_enhancers_only_MOG, number_of_samples_correl = number_of_samples_correl, burn_in_correl = burn_in_correl)


