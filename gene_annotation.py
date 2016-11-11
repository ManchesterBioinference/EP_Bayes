
data_folder = config_variables.data_folder

#genes = np.loadtxt(data_folder + "Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed.gz", dtype = str)[:,3]

gene_names = np.loadtxt(config_variables.name_of_time_series_promoter_file_for_TSS_start, dtype = str, usecols = (0,1,2,3))

labels = np.loadtxt("./R_scripts/AP_clustering_output/Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_ER_100_labels", dtype = str)

labels = np.array([el.split(",")[1] for el in labels])[1:].astype(int)

survived_indexes = np.loadtxt("./R_scripts/AP_clustering_output/Homo_sapiens.GRCh37.75.gtf_filtered_gene_joint_2_cleaned_chr_sorted_sorted_ordered_0_indexed_300_unfiltered_count_concat_ER_100_survived_indexes").astype(int)



survived_genes = gene_names[survived_indexes]

for label in range(min(labels),max(labels)):

	genes_in_cluster = survived_genes[labels == label]

	np.savetxt("R_scripts/AP_clustering_output/genes_in_cluster/genes_in_cluster_{0}.txt".format(label), genes_in_cluster, fmt ="%s")	


Predicted_genes = np.loadtxt("./results/clusters_genes_vs_counts_prob_distant_all_0.2_PolII,ER_smo_True_proximal_version_PR_met_300_0_300", dtype = str)

Predicted_genes_coordinates = gene_names[np.in1d(gene_names[:,-1], Predicted_genes[Predicted_genes[:,4].astype(float) > 0][:,0])]

np.savetxt("./results/Predicted_genes_coordinates_gene_annotation", Predicted_genes_coordinates, fmt = "%s")

