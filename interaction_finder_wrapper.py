import os
import shlex, subprocess
import config_variables

data_folder = config_variables.data_folder
temp_output = config_variables.temp_output

name_of_promoter_file_for_overlap = config_variables.name_of_promoter_file_for_overlap
name_of_enhancer_file_for_overlap = config_variables.name_of_enhancer_file_for_overlap

upstream = config_variables.upstream
downstream = config_variables.downstream
promoter_overlaps_enhancer_file = config_variables.promoter_overlaps_enhancer_file


command_line = "windowBed -a {0} -b {1} -sw -l {2} -r {3} > {4}".format(name_of_promoter_file_for_overlap, name_of_enhancer_file_for_overlap, upstream, downstream, promoter_overlaps_enhancer_file)
os.system(command_line)

name_of_overlap_file_pro = config_variables.name_of_overlap_file_pro
name_of_overlap_file_enh = config_variables.name_of_overlap_file_enh

os.system("cut -f1-5 {0} > {1}".format(promoter_overlaps_enhancer_file, name_of_overlap_file_pro))
os.system("cut -f6-11 {0} > {1}".format(promoter_overlaps_enhancer_file, name_of_overlap_file_enh))

import non_ER_promoters_clean
import non_proximal_enhancer_clean

os.system("cat {0}Non_ER_controled_promoters_pindexed.txt {0}ER_controled_promoters_pindexed.txt {0}distal_ER_peaks_pindexed.txt > {0}all_features_without_distance_ones".format(temp_output))

os.system("sort -k 1,1d -k 2,2n -k 3,3n {0}all_features_without_distance_ones > {0}all_features_without_distance_ones_s".format(temp_output))

import orderer_clean
orderer_clean.executor("{0}all_features_without_distance_ones_s".format(temp_output), config_variables.chrom_names)

os.system("pairToBed -a {0}merged_interactions_CHiA_PET_ordered_0_ind.gz -b {1}all_features_without_distance_ones_s_ordered -type both > {1}all_features_s_ordered_concat_new_interactions_{2}_{3}".format(data_folder, temp_output, upstream, downstream))

import interaction_checker_changed_clean
interaction_checker_changed_clean.executor("{0}all_features_s_ordered_concat_new_interactions_{1}_{2}".format(temp_output, upstream, downstream), upstream, downstream, temp_output)



