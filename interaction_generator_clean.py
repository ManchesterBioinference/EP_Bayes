def generator(pro_survived, enh_survived, domain, max_path):

	import copy
	import numpy as np
	import re
	import config_variables
	promoter_overlaps_enhancer_file = config_variables.promoter_overlaps_enhancer_file
	upstream = config_variables.upstream
	downstream = config_variables.downstream
	link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
	dataset_time_series_dict = config_variables.dataset_time_series_dict
	TSS_or_intra_genic_for_domain_filter = config_variables.TSS_or_intra_genic_for_domain_filter
	name_of_time_series_promoter_file_for_TSS_start = config_variables.name_of_time_series_promoter_file_for_TSS_start
	temp_output = config_variables.temp_output
	
	#parameters-------------------------
	#ovenh_ovenh_pro_pro_version = False
	#max_pro_enh_mode = True
	#-----------------------------------

	enhancer_enhancer_inter = np.loadtxt(temp_output + 'enhancer_enhancer_interactions_{0}_{1}'.format(upstream, downstream), usecols = (0,1,2), dtype = str, delimiter = '\t') 
	promoter_promoter_inter = np.loadtxt(temp_output + 'promoter_promoter_interactions_{0}_{1}'.format(upstream, downstream), usecols = (0,1,2), dtype = str, delimiter = '\t') 
	promoter_enhancer_inter = np.loadtxt(temp_output + 'promoter_enhancer_interactions_{0}_{1}'.format(upstream, downstream), usecols = (0,1,2), dtype = str, delimiter = '\t') 

	un_stringer = lambda x: int(re.findall('\d+', x)[0])
	
	def un_featurer(array): f = lambda x: re.findall('\D+', x)[0]; return np.array(map(f, array))


	def un_string(array_): return np.c_[np.array(map(un_stringer, array_[:, 0]))[:, None], np.array(map(un_stringer, array_[:, 1]))[:, None]]

	enh_enh_indexes_list = un_string(enhancer_enhancer_inter[:,1:])
	pro_enh_indexes_list = un_string(promoter_enhancer_inter[:,1:])
	pro_pro_indexes_list = un_string(promoter_promoter_inter[:,1:])


	def filter_(array_, filt_1, filt_2): return np.in1d(map(un_stringer, array_[:, 1]), filt_1) * np.in1d(map(un_stringer, array_[:, 2]), filt_2)

	#cleans interactions:
	#----------------------------------------------------------------------------------------------------------------------------------------
	promoter_promoter_inter = promoter_promoter_inter[filter_(promoter_promoter_inter, pro_survived, pro_survived)]
	enhancer_enhancer_inter = enhancer_enhancer_inter[filter_(enhancer_enhancer_inter, enh_survived, enh_survived)]
	promoter_enhancer_inter = promoter_enhancer_inter[filter_(promoter_enhancer_inter, pro_survived, enh_survived)]

	#----------------------------------------------------------------------------------------------------------------------------------------

	enh_enh_indexes_list = un_string(enhancer_enhancer_inter[:,1:])
	pro_enh_indexes_list = un_string(promoter_enhancer_inter[:,1:])
	pro_pro_indexes_list = un_string(promoter_promoter_inter[:,1:])



	def prepare_overlaps():
		
		overlaps = np.loadtxt(promoter_overlaps_enhancer_file, delimiter = '\t', usecols = (4, 8), dtype = int)
		overlaps_promoter_enhancer_inter = overlaps 
		diction_overlaps_ovenh = {}

		for overl in overlaps[:, 1]: promoters = list(overlaps[overl == overlaps[:, 1], 0]); diction_overlaps_ovenh[overl] = promoters

		return diction_overlaps_ovenh

	diction_overlaps_ovenh = prepare_overlaps()


	def promoter_promoter_adder(diction_overlaps_ovenh, index_1, index_2, chro, pro_pro_indexes_list, promoter_promoter_inter):
		
		print 'ovenh-ovenh'
		promoters_1 = diction_overlaps_ovenh[index_1]
		promoters_2 = diction_overlaps_ovenh[index_2]

		promoters_1 = promoters_1[np.in1d(promoters_1, pro_survived)] # converts ER signals which overlap filtered out promoters into distant peaks
		promoters_2 = promoters_2[np.in1d(promoters_2, pro_survived)]

		legend = np.r_[promoters_1, promoters_2]

		if len(legend) > 1:
			counts, bins = np.histogram(legend, np.arange(min(legend), max(legend) + 2))
			digitize = np.digitize(legend, np.arange(min(legend), max(legend) + 2)) - 1
			if len(np.where(counts > 1)[0]): print 'ambigous allocation - peak overlaps two promoters which interact through ChIA-PET index {0}, {1}'.format(index_1, index_2)
		
			matrix = np.ones((len(legend), len(legend)), dtype = bool)			
			matrix[:len(promoters_1), :len(promoters_1)] = False
			matrix[len(promoters_1): len(promoters_1) + len(promoters_2), len(promoters_1): len(promoters_1) + len(promoters_2)] = False
			matrix[np.tril_indices(len(legend))] = False # upper triangular elements
			pro_inter = np.c_[legend[np.where(matrix)[0]], legend[np.where(matrix)[1]]]
			pro_inter_symbolic = [[chro ,'ovpro{0}'.format(ind1), 'ovpro{0}'.format(ind2)] for ind1, ind2 in pro_inter]

			for el_1, el_2 in zip(pro_inter, pro_inter_symbolic):
				lista = list(el_1)
				if lista not in map(list, pro_pro_indexes_list):
					pro_pro_indexes_list = np.r_[pro_pro_indexes_list, el_1]
					promoter_promoter_inter = np.r_[promoter_promoter_inter, el_2]

		return pro_pro_indexes_list, promoter_promoter_inter
				

	def ER_enhancer_Non_enhancer_pro_adder(pro_enh_indexes_list, promoter_enhancer_inter, pro_pro_indexes_list, promoter_promoter_inter):

		enhancer_enhancer_inter_filtered = []
		pro_enh_indexes_list_added = []
		pro_enh_indexes_list_symbolic_added = []

		for el in enhancer_enhancer_inter:
			chro = el[0]
			index_1 = int(un_stringer(el[1]))
			index_2 = int(un_stringer(el[2]))
			feature_1 = re.findall('\D+', el[1])[0]
			feature_2 = re.findall('\D+', el[2])[0]

			if feature_1 == 'ovenh' and feature_2 == 'enh':	
				for pro_dict in diction_overlaps_ovenh[index_1]: # takes promoters corresponding to ER overlapping with them, previus step will clean out overnh which doesnt have signal
					pro_enh_int = [pro_dict, index_2]
					pro_enh_int_symbolic = [chro, 'ovpro{0}'.format(pro_dict), 'enh{0}'.format(index_2)]
					if pro_enh_int not in map(list, pro_enh_indexes_list) and pro_dict in pro_survived:
						pro_enh_indexes_list_added += [pro_enh_int]
						pro_enh_indexes_list_symbolic_added += [pro_enh_int_symbolic]

			elif feature_1 == 'enh' and feature_2 == 'ovenh':	# here you've got to check what index has got the promoter.
				for pro_dict in diction_overlaps_ovenh[index_2]:
					pro_enh_int = [pro_dict, index_1] 				
					pro_enh_int_symbolic = [chro, 'ovpro{0}'.format(pro_dict), 'enh{0}'.format(index_1)]
					if pro_enh_int not in map(list, pro_enh_indexes_list) and pro_dict in pro_survived:
						pro_enh_indexes_list_added += [pro_enh_int]
						pro_enh_indexes_list_symbolic_added += [pro_enh_int_symbolic]

			elif feature_1 == 'ovenh' and feature_2 == 'ovenh':

				if ovenh_ovenh_pro_pro_version:	pro_pro_indexes_list, promoter_promoter_inter = promoter_promoter_adder(diction_overlaps_ovenh, index_1, index_2, chro, pro_pro_indexes_list, promoter_promoter_inter) # converts ovenh-ovenh interactions to promoter-promoter interactions
				else: enhancer_enhancer_inter_filtered.append(el) # keeps ovenh-ovenh interactions as enh-enh interactions

			else: 	
				enhancer_enhancer_inter_filtered.append(el)

		if len(pro_enh_indexes_list_added): 
			pro_enh_indexes_list = np.r_[pro_enh_indexes_list, pro_enh_indexes_list_added]
			promoter_enhancer_inter = np.r_[promoter_enhancer_inter, pro_enh_indexes_list_symbolic_added]

		return pro_enh_indexes_list, promoter_enhancer_inter, pro_pro_indexes_list, promoter_promoter_inter

	#pro_enh_indexes_list, promoter_enhancer_inter, pro_pro_indexes_list, promoter_promoter_inter = ER_enhancer_Non_enhancer_pro_adder(pro_enh_indexes_list, promoter_enhancer_inter, pro_pro_indexes_list, promoter_promoter_inter) # it doesn't do a thing if you take only distal enhancers. Overlapping enhancers are then not among distal enhancers which you allow to form with with promoters. 


	def ER_ovenh_pro_pro_adder():

		filtered_promoter_enhancer_inter = []
		overlaps = np.loadtxt(promoter_overlaps_enhancer_file, delimiter = '\t', usecols = (4, 8), dtype = int)
		overlaps_promoter_enhancer_inter = overlaps 
		diction_overlaps_ovenh = {}

		for overl in overlaps[:, 1]:	

			promoters = list(overlaps[overl == overlaps[:, 1], 0])
			diction_overlaps_ovenh[overl] = promoters

		for el in promoter_enhancer_inter:
			
			index_1 = int(re.findall('\d+', el[1])[0])
			index_2 = int(re.findall('\d+', el[2])[0])
			feature_2 = re.findall('\D+', el[1])[0]

			if feature_2 == 'ovenh':	
				for pro_dict in diction_overlaps_ovenh[index_2]:
					#pro_pro_int = [index_1, dict_pro_survived[pro_dict]] 
					pro_pro_int = [index_1, pro_dict]				
					if pro_pro_int not in filtered_promoter_promoter_inter and pro_dict in pro_survived:
						filtered_promoter_promoter_inter.append(pro_pro_int)

			else:
				filtered_promoter_enhancer_inter.append(el)
		return np.array(filtered_promoter_enhancer_inter)

	#if not(max_pro_enh_mode): promoter_enhancer_inter = ER_ovenh_pro_pro_adder()

	#if the aim is to maximise number of promoter-enhancer interactions it could be best to set it to false

	def scan_through_pro_pro_inter_and_add_pro_enh_inter(): # it's kind of arbitrary what we are trying to get rid of.. so this part can actually be coded, so that when pro-pro then pro-enh and enh-enh... but it may be best to do it with interaction matrix  
		pass

	def into_string(array, s1, s2): 
		array_ = un_string(array[:,1:])
		array_ = np.c_[array[:,0][:,None], np.array([s1 + str(index) for index in array_[:,0]])[:, None], np.array([s2 + str(index) for index in array_[:,1]])[:, None]]	
		return array_

	promoter_enhancer_inter = into_string(promoter_enhancer_inter,'p', 'e')
	promoter_promoter_inter = into_string(promoter_promoter_inter,'p', 'p')
	enhancer_enhancer_inter = into_string(enhancer_enhancer_inter,'e', 'e')


	#--------------------------------------------------------------------------------------------------------------------------------------------------------


	

	def stringer(array, st): return np.array(map(lambda x: '{0}{1}'.format(st,x), array))


	def filter_domains(symb_1, symb_2, chr_interactions, length_chr, row_indexes_plus, column_indexes_plus):

		if len(chr_interactions) > 0:	

			interaction_matrix = np.zeros((length_chr, length_chr), bool)
			true_indexes = un_string(chr_interactions)
			
			interaction_matrix[true_indexes[:,0] + row_indexes_plus, true_indexes[:,1] + column_indexes_plus] = True
			interaction_matrix[true_indexes[:,1] + column_indexes_plus, true_indexes[:,0] + row_indexes_plus] = True
			interaction_matrix[np.tril_indices(length_chr)] = False # gets rid of symmetric interactions		
			print 'number of original: {0} {1} - {2}  true survived domain filtering: {3}'.format(len(true_indexes), symb_1, symb_2, np.sum(interaction_matrix*domain_matrix))
			domain_raws = (np.where(interaction_matrix*domain_matrix)[0] - row_indexes_plus)
			domain_columns = (np.where(interaction_matrix*domain_matrix)[1] - column_indexes_plus)
			true_indexes_dom_symb = np.c_[stringer(domain_raws, symb_1)[:, None], stringer(domain_columns, symb_2)[:, None]]
			s_1 = len(true_indexes_dom_symb)
			return np.c_[np.array(['chr{0}'.format(i)]*s_1)[:,None], true_indexes_dom_symb]
			
		else: return []	

	if domain:
		import interacting_domain
		def extract_TSS_coordinates(upstream):

			data = np.loadtxt(name_of_time_series_promoter_file_for_TSS_start, dtype = str,  delimiter = '\t')	
			plus_strand = data[:, 4] == '+'
			TSS_coordinates = np.zeros(len(plus_strand), int)
			TSS_coordinates[plus_strand] = data[plus_strand, 1].astype(int) - upstream
			TSS_coordinates[np.invert(plus_strand)] = data[np.invert(plus_strand), 2].astype(int) + upstream

			return TSS_coordinates

		TSS_coordinates = extract_TSS_coordinates(upstream)

		def initialise_variables(chrom):

			name_of_pro_t_s = link_data_set_name_to_file_name["promoters"]["ER"]
			name_of_enh_t_s = link_data_set_name_to_file_name["enhancers"]["ER"]

			pro_chroms, pro_coordinates, ts_p = dataset_time_series_dict[name_of_pro_t_s]
			enh_chroms, enh_coordinates, ts_e = dataset_time_series_dict[name_of_enh_t_s]

			indexes_p = np.where(pro_chroms==chrom)[0]	# gives the number of promoters for a chromosome	
			indexes_e = np.where(enh_chroms==chrom)[0]	# gives the number of enhancers for a chromosome

			return pro_chroms, enh_chroms, pro_coordinates, enh_coordinates, indexes_p, indexes_e

		enhancer_enhancer_inter_dom = {}
		promoter_enhancer_inter_dom = {}
		promoter_promoter_inter_dom = {}	
		total_p = total_e = 0
		for i in np.concatenate((np.array(range(1, 23), dtype='S2'), ['X'], ['Y'])):

			chr_interactions_enh_enh = enhancer_enhancer_inter[enhancer_enhancer_inter[:,0]=='chr{0}'.format(i)][:,1:]
			chr_interactions_pro_enh = promoter_enhancer_inter[promoter_enhancer_inter[:,0]=='chr{0}'.format(i)][:,1:]
			chr_interactions_pro_pro = promoter_promoter_inter[promoter_promoter_inter[:,0]=='chr{0}'.format(i)][:,1:]

			if not(len(chr_interactions_enh_enh)*len(chr_interactions_pro_enh)*len(chr_interactions_pro_pro)): continue

			pro_chroms, enh_chroms, pro_coord, enh_coord, indexes_p, indexes_e = initialise_variables('chr{0}'.format(i))

			if TSS_or_intra_genic_for_domain_filter == "TSS_only":
				pro_coord = np.column_stack((TSS_coordinates, TSS_coordinates + 2))
			
			length_chr = len(indexes_p) + len(indexes_e)
			chrom_pro_coord = pro_coord[indexes_p]
			chrom_ER_coord = enh_coord[indexes_e]
		
			domain_matrix = interacting_domain.interacting_domains(chrom_pro_coord, chrom_ER_coord, 'chr{0}'.format(i), 'left')
			domain_matrix = domain_matrix + interacting_domain.interacting_domains(chrom_pro_coord, chrom_ER_coord, 'chr{0}'.format(i), 'right')				

			promoter_enhancer_inter_dom['chr{0}'.format(i)] = filter_domains('p','e', chr_interactions_pro_enh, length_chr, - total_p, - total_e + len(indexes_p)) 
			enhancer_enhancer_inter_dom['chr{0}'.format(i)] = filter_domains('e','e', chr_interactions_enh_enh, length_chr, - total_e + len(indexes_p), - total_e + len(indexes_p)) 
			promoter_promoter_inter_dom['chr{0}'.format(i)] = filter_domains('p','p', chr_interactions_pro_pro, length_chr, - total_p, - total_p)
			# interaction_domains_adjustments end -----------------------

			total_p += len(indexes_p)
			total_e += len(indexes_e)
	
		promoter_enhancer_inter_dom['chrY'] = []
		enhancer_enhancer_inter_dom['chrY'] = []
		promoter_promoter_inter_dom['chrY'] = []
			
		enhancer_enhancer_inter = enhancer_enhancer_inter_dom['chr1']
		promoter_enhancer_inter = promoter_enhancer_inter_dom['chr1']
		promoter_promoter_inter = promoter_promoter_inter_dom['chr1']
		for i in np.r_[np.arange(2, 23).astype('S2'), ['X'], ['Y']]: 
			add_1 = promoter_enhancer_inter_dom['chr{0}'.format(i)]
			if len(add_1): promoter_enhancer_inter = np.r_[promoter_enhancer_inter, add_1]
			add_2 = enhancer_enhancer_inter_dom['chr{0}'.format(i)]
			if len(add_2): enhancer_enhancer_inter = np.r_[enhancer_enhancer_inter, add_2]
			add_3 = promoter_promoter_inter_dom['chr{0}'.format(i)]
			if len(add_3): promoter_promoter_inter = np.r_[promoter_promoter_inter, add_3]

		interactions_to_save = np.r_[promoter_enhancer_inter, enhancer_enhancer_inter]
		interactions_to_save = np.c_[interactions_to_save, np.ones(len(interactions_to_save))[:, None]]

	#if not(generate_intermediates): return interactions_of_path

	def prepares_reverse_map_and_uniqueness():

		promoters = np.r_[promoter_promoter_inter[:,[0,1]], promoter_promoter_inter[:,[0,2]], promoter_enhancer_inter[:,[0,1]]]
		promoters = np.array(list(set(map(tuple, promoters))))
		promoters_sort_indexes = np.argsort(map(un_stringer, promoters[:,1]))
		promoters = promoters[promoters_sort_indexes]

		enhancers = np.r_[enhancer_enhancer_inter[:,[0,1]], enhancer_enhancer_inter[:,[0,2]], promoter_enhancer_inter[:,[0,2]]]
		enhancers = np.array(list(set(map(tuple, enhancers))))
		enhancers_sort_indexes = np.argsort(map(un_stringer, enhancers[:,1]))
		enhancers = enhancers[enhancers_sort_indexes]

		pro_enh_unique = np.r_[promoters, enhancers]

		chroms_frame, pro_enh_unique_ordered = [], []
		for i in np.r_[np.arange(1, 23).astype('S2'), ['X'], ['Y']]: chr_mask = pro_enh_unique[:,0] == 'chr{0}'.format(i); chroms, chroms_features = pro_enh_unique[chr_mask,0], pro_enh_unique[chr_mask,1]; chroms_frame = np.r_[chroms_frame, chroms]; pro_enh_unique_ordered = np.r_[pro_enh_unique_ordered, chroms_features]

		dict_inter = {}
		for index, el in enumerate(pro_enh_unique_ordered): dict_inter[el] = index
		return pro_enh_unique_ordered, chroms_frame, dict_inter

	unique_features, chroms_frame, dict_inter = prepares_reverse_map_and_uniqueness()

	def difference_interactions_prod(difference_arr, path):
		indexes_of_lower_diagonal = np.tril_indices(len(unique_features))
		difference_arr[indexes_of_lower_diagonal] = False
		#difference_arr[range(len(unique_features)), range(len(unique_features))] = False
		indexes_of_non_zero_interactions = np.where(difference_arr)

		column_chr = chroms_frame[indexes_of_non_zero_interactions[0]].astype(str)
		column_1st = unique_features[indexes_of_non_zero_interactions[0]]
		column_2nd = unique_features[indexes_of_non_zero_interactions[1]]
		column_path = np.array([path]*len(indexes_of_non_zero_interactions[0]), str)
		

		diff_interactions = np.c_[column_chr[:,None], column_1st[:,None], column_2nd[:,None], column_path[:,None]]
		#print diff_interactions
		return diff_interactions

	a, b = [], []

	matrix_of_interactions = np.zeros((len(unique_features),len(unique_features)), bool)
	matrix_of_interactions[range(len(unique_features)), range(len(unique_features))] = True
	cumulative_old = np.zeros((len(unique_features),len(unique_features)), bool)
	concat = np.r_[promoter_enhancer_inter, enhancer_enhancer_inter, promoter_promoter_inter]

	indexes_1, indexes_2 = [], []
	for chr_, feature_1, feature_2 in concat: indexes_1.append(dict_inter[feature_1]), indexes_2.append(dict_inter[feature_2])
	matrix_of_interactions[[indexes_1, indexes_2], [indexes_2, indexes_1]] = True

	matrix_of_interactions_so_far = matrix_of_interactions
	cumulative = matrix_of_interactions
	collective = cumulative.astype(int)

	#ll = [el for el in map(list, concat) if el not in map(list, interactions_of_path[:,:3])]

	path = 1
	interactions_of_path = difference_interactions_prod(cumulative-cumulative_old, path) 

	def gets_rid_of_promoter_promoter_inter(array_): return np.prod(np.c_[un_featurer(array_[:,0])[:,None], un_featurer(array_[:,1])[:,None]] == ['p','p'], axis = 1) == False

	print 'path:' ,path, 'size =', (sum(sum(np.array(matrix_of_interactions_so_far))) - len(unique_features))/2

	while not(np.array_equal(cumulative, cumulative_old)):

		path += 1
		
		cumulative_old = cumulative
		matrix_of_interactions_so_far = np.dot(matrix_of_interactions, matrix_of_interactions_so_far)
		cumulative = matrix_of_interactions_so_far
		difference = cumulative - cumulative_old
		collective = collective + difference.astype(int)*path

		print 'path:' ,path, 'size =', (sum(sum(np.array(matrix_of_interactions_so_far))) - len(unique_features))/2

		if path == max_path: break
		
		interactions_of_path = np.r_[interactions_of_path, difference_interactions_prod(difference, path)]

	mask = gets_rid_of_promoter_promoter_inter(interactions_of_path[:,[1,2]])
	interactions_of_path = interactions_of_path[mask]	

	return interactions_of_path


