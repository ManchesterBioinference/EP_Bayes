args <- commandArgs(TRUE)

library(matrixStats)
cwd <- getwd()

update_cwd = paste(cwd, "AP_clustering_output" ,sep = "/")
print(update_cwd)
setwd(update_cwd)

print(getwd())

print(args[1])
print(args[2])
require(methods)

name_of_file_to_cluster <- args[1]
number_of_clusters <- as.integer(args[2])

library("apcluster")
library('matrixStats')
data = read.table(name_of_file_to_cluster, sep=",")
print("calculates distance matrix")
S = negDistMat(data,r=2)
print("starts clustering")
res = apclusterK(S,K=number_of_clusters)
labels_ = labels(res, type="enum")
print("ends clustering")
#"common_region_peaks_extended_less_time_points_sorted_unfiltered_count_concat_ER_PolII_300_-1.0_labels"
clustered_data_set_and_its_labels = paste(name_of_file_to_cluster, "labels", sep="_")
print("saves labels")
write.csv(labels_, file = clustered_data_set_and_its_labels)

