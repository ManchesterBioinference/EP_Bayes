import numpy as np
import config_variables
temp_output = config_variables.temp_output
name_of_overlap_file_enh = config_variables.name_of_overlap_file_enh
name_of_enhancer_file_for_overlap = config_variables.name_of_enhancer_file_for_overlap

ER_peaks = np.loadtxt(name_of_enhancer_file_for_overlap, usecols = range(4), dtype = str)

ER_proximal_binding_sites_indexes = np.loadtxt(name_of_overlap_file_enh, dtype = int, usecols =(3,), delimiter = '\t')

unique_ER_proximal_binding_sites = ER_peaks[np.unique(ER_proximal_binding_sites_indexes)]

ER_distal_binding_sites_indexes = np.ones(len(ER_peaks), bool)
ER_distal_binding_sites_indexes[ER_proximal_binding_sites_indexes] = False
ER_distal_binding_sites = ER_peaks[ER_distal_binding_sites_indexes]


column = np.core.defchararray.add('enh', ER_distal_binding_sites[:, 3])

test = np.concatenate((ER_distal_binding_sites[:,[0,1,2]], column[:, None]), axis = 1)

np.savetxt(temp_output + 'distal_ER_peaks_pindexed.txt', test, fmt='%s', delimiter='\t')

column_2 = np.core.defchararray.add('ovenh', unique_ER_proximal_binding_sites[:, 3])

test_2 = np.concatenate((unique_ER_proximal_binding_sites[:, [0,1,2]], column_2[:, None]), axis = 1)

np.savetxt(temp_output + 'proximal_ER_peaks_pindexed.txt', test_2, fmt='%s', delimiter='\t')

