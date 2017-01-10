#The code will produce figures 2, 3b, 3d, 3f 
import graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
ENHANCER_GENE_INTERACTOR = graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned
#graph_plotter_all_PR_PCA_true_only_MAP_per_million_multi_chroms_multi_filters_non_domain_inter_background_pro_dist_alt_mixer_cleaned.executor(redo_raw_CHIA_PET_interactions = False, upstream = 300)


mode_of_code = ["ODD","EVEN", "FULL", "GAUSSIAN_SAMPLE", "MK_PAIRWISE", "CLUSTERING"][-1]

re_do_clustering = False #choose whether you want to redo clustering or just plot figures based on the labels from previous clustering.
plot_TF_enrichments_in_cluster = False

cluster_figure_selections = ["cluster_ER_enhancer", "cluster_Pol2s_enhancer", "cluster_Pol2s_promoter", "cluster_ER_promoter", "cluster_Pol2s_ER_enhancer", "cluster_Pol2s_ER_enhancer_test"]

for cluster_figure_selection in cluster_figure_selections:

	ENHANCER_GENE_INTERACTOR.executor(mode_of_code = mode_of_code, cluster_figure_selection = cluster_figure_selection, re_do_clustering = re_do_clustering, plot_TF_enrichments_in_cluster = plot_TF_enrichments_in_cluster)


#plots enrichments
cluster_figure_selection = cluster_figure_selections[0] #any really, it doesn't matter for the task
plot_TF_enrichments_in_cluster = True
ENHANCER_GENE_INTERACTOR.executor(mode_of_code = mode_of_code, cluster_figure_selection = cluster_figure_selection, re_do_clustering = re_do_clustering, plot_TF_enrichments_in_cluster = plot_TF_enrichments_in_cluster)


