import numpy as np
import config_variables
name_of_promoter_file_for_overlap = config_variables.name_of_promoter_file_for_overlap
name_of_overlap_file_pro = config_variables.name_of_overlap_file_pro

Promoters = np.loadtxt(name_of_promoter_file_for_overlap, dtype = str)
ER_controled_promoters_indexes = np.loadtxt(name_of_overlap_file_pro, dtype = int, usecols =(4,), delimiter = '\t')

unique_ER_controled_promoters = Promoters[np.unique(ER_controled_promoters_indexes)]

Non_ER_controled_promoters_indexes = np.ones(len(Promoters), bool)
Non_ER_controled_promoters_indexes[ER_controled_promoters_indexes] = False

Non_ER_controled_promoters = Promoters[Non_ER_controled_promoters_indexes]

column = np.core.defchararray.add('ERpro', unique_ER_controled_promoters[:, 4])

test = np.concatenate((unique_ER_controled_promoters[:, [0,1,2]], column[:, None]), axis = 1)

np.savetxt('ER_controled_promoters_pindexed.txt', test, fmt='%s', delimiter='\t')

column_2 = np.core.defchararray.add('nERpro', Non_ER_controled_promoters[:, 4])

test_2 = np.concatenate((Non_ER_controled_promoters[:, [0,1,2]], column_2[:, None]), axis = 1)

np.savetxt('Non_ER_controled_promoters_pindexed.txt', test_2, fmt='%s', delimiter='\t')

