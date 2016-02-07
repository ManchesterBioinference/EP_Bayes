def executor(interactions_to_dissect, upstream, downstream):
	import numpy as np
	import re
	import itertools
	import config_variables
	temp_output = config_variables.temp_output

	int_flts=np.loadtxt(interactions_to_dissect, usecols=(1,2,4,5,8,9), dtype=int)
	int_strs=np.loadtxt(interactions_to_dissect, usecols=(0,3,6,7,10), dtype=str)
	lines = np.loadtxt(interactions_to_dissect, dtype=str)

	#full = open('ER_promoters_and_distal_ER_peaks_concat_ER_concat_new_interactions_column_changed.bed', 'w')
	#clean = open('ER_promoters_and_distal_ER_peaks_concat_ER_concat_new_interactions_column_changed_cleaned.bed', 'w')
	#inter_only = open('pure_true_interactions', 'w')
	#clean_inter_only = open('pure_true_interactions_cleaned', 'w')
	#clean_inter_only_single = open('pure_true_interactions_cleaned_single', 'w')
	removed = open(temp_output + 'removed_both_and_wrong_{0}_{1}'.format(upstream, downstream), 'w')


	promoter_promoter_interactions = open(temp_output + 'promoter_promoter_interactions_{0}_{1}'.format(upstream, downstream), 'w')
	promoter_enhancer_interactions = open(temp_output + 'promoter_enhancer_interactions_{0}_{1}'.format(upstream, downstream), 'w')
	enhancer_enhancer_interactions = open(temp_output + 'enhancer_enhancer_interactions_{0}_{1}'.format(upstream, downstream), 'w')



	index = -1
	last_index = np.shape(int_flts)[0]-1


	def check_interaction_type_and_save(chrom, index_1, index_2, save_2):

		if 'ERpro'==re.match('\D+', index_1).group():

			if 'enh'==re.match('\D+', index_2).group():
				#save_2 = '{0}\t{1}\t{2}\n'.format(chrom, index_1, index_2)
				promoter_enhancer_interactions.write(save_2 + '\n')
			
			elif 'ovenh'==re.match('\D+', index_2).group():
				#save_2 = '{0}\t{1}\t{2}\n'.format(chrom, index_1, index_2)
				promoter_enhancer_interactions.write(save_2 + '\n')

			elif 'nERpro'==re.match('\D+', index_2).group():
				promoter_promoter_interactions.write(save_2 + '\n')

			elif 'ERpro'==re.match('\D+', index_2).group():
				promoter_promoter_interactions.write(save_2 + '\n')
		
		elif 'nERpro'==re.match('\D+', index_1).group():

			if 'enh'==re.match('\D+', index_2).group():
				#save_2 = '{0}\t{1}\t{2}\n'.format(chrom, index_1, index_2)
				promoter_enhancer_interactions.write(save_2 + '\n')
			
			elif 'ovenh'==re.match('\D+', index_2).group():
				#save_2 = '{0}\t{1}\t{2}\n'.format(chrom, index_1, index_2)
				promoter_enhancer_interactions.write(save_2 + '\n')

			elif 'nERpro'==re.match('\D+', index_2).group():
				promoter_promoter_interactions.write(save_2 + '\n')

			elif 'ERpro'==re.match('\D+', index_2).group():
				promoter_promoter_interactions.write(save_2 + '\n')


		elif 'enh'==re.match('\D+', index_1).group():

			if 'enh'==re.match('\D+', index_2).group():
			
				enhancer_enhancer_interactions.write(save_2 + '\n')
			
			elif 'ovenh'==re.match('\D+', index_2).group():
			
				enhancer_enhancer_interactions.write(save_2 + '\n')	

			elif 'nERpro'==re.match('\D+', index_2).group():
				save_2 = '{0}\t{2}\t{1}'.format(chrom, index_1, index_2)			
				promoter_enhancer_interactions.write(save_2 + '\n')

			elif 'ERpro'==re.match('\D+', index_2).group():
				save_2 = '{0}\t{2}\t{1}'.format(chrom, index_1, index_2)
				promoter_enhancer_interactions.write(save_2 + '\n')

		elif 'ovenh'==re.match('\D+', index_1).group():

			if 'enh'==re.match('\D+', index_2).group():
			
				enhancer_enhancer_interactions.write(save_2 + '\n')
			
			elif 'ovenh'==re.match('\D+', index_2).group():
			
				enhancer_enhancer_interactions.write(save_2 + '\n')

			elif 'nERpro'==re.match('\D+', index_2).group():
				save_2 = '{0}\t{2}\t{1}'.format(chrom, index_1, index_2)			
				promoter_enhancer_interactions.write(save_2 + '\n')

			elif 'ERpro'==re.match('\D+', index_2).group():
				save_2 = '{0}\t{2}\t{1}'.format(chrom, index_1, index_2)
				promoter_enhancer_interactions.write(save_2 + '\n')


		else:
			print 'something wrong'

		if index_1 == index_2:
			print 'something_wrong', save_2


	string_checker=[]
	failed =[]
	while index < last_index-1: #look at ** comment below 
		lista=[]
	
		while index < last_index-1: # last index in the array -1.
		
			index+=1
		
			if int_strs[index, 2] == int_strs[index + 1, 2]:
				lista.append(index)
				if index==last_index-1:
					lista.append(last_index)
			else:
				inter = int_strs[index, 2]
				lista.append(index)
				break
	
		current_interaction = int_strs[index, 2]		
	
		lista_2=[]		
	
		for el_index in lista:

			set_left = set(range(int_flts[el_index,0], int_flts[el_index,1]+1))
			set_right = set(range(int_flts[el_index,2], int_flts[el_index,3]+1))
			set_peak = set(range(int_flts[el_index,4], int_flts[el_index,5]+1))

			cardinality_1 = len(set_left & set_peak)
			cardinality_2 = len(set_right & set_peak)		

			if cardinality_1 <> 0 and cardinality_2==0:
				state=(el_index,'left')

			elif cardinality_2 <> 0 and cardinality_1==0:	
				state=(el_index,'right')

			elif cardinality_2 <> 0 and cardinality_1<>0:	
				state=(el_index,'both')

			else:
				state=(el_index,'wrong')
				print state
			lista_2.append(state)

		for el in itertools.combinations(lista_2, 2):

			#save = '{0[0]}\t{1[0]}\t{1[1]}\t{0[1]}\t{1[2]}\t{1[3]}\t{0[2]}\t{0[3]}\t{1[4]}\t{1[5]}\t{0[4]}'.format(int_strs[el[0][0]], int_flts[el[0][0]]) # reconstruct the line in the database
		
						
			if (el[0][1]=='right' and el[1][1]=='left') or (el[0][1]=='left' and el[1][1]=='right'):
				
				#save = '{0}\t{1}\t{2}'.format(save, int_strs[el[0][0], 4], int_strs[el[1][0], 4])
						
				save_2 = '{0}\t{1}\t{2}'.format(int_strs[el[0][0], 3], int_strs[el[0][0], 4], int_strs[el[1][0], 4])

				anty_save_2 = '{0}\t{2}\t{1}'.format(int_strs[el[0][0], 3], int_strs[el[0][0], 4], int_strs[el[1][0], 4])

				#clean_inter_only.write(save_2 + '\n')	# saves only the lines which satisfied the logical statement.		

				#inter_only.write(save_2 + '\n') # saves the lines which satisfied the logical statement or the two below.

				if save_2 == anty_save_2:
					print 'something_wrong', save_2, el 
					failed.append(lista_2)


				if save_2 not in string_checker and anty_save_2 not in string_checker: 
	
					string_checker.append(save_2)	
					string_checker.append(anty_save_2)
					#clean_inter_only_single.write(save_2 + '\n') #saves only lines which haven't occured before. (anty interactions)
					check_interaction_type_and_save(int_strs[el[0][0], 3],int_strs[el[0][0], 4], int_strs[el[1][0], 4], save_2) # checks a type of interaction and saves to interaction type file
				
					#test_clean = True

								
			elif el[0][1]=='both' or el[1][1]=='both':
				
				#save = '{0}\t{1}\t{2}\t{3}\t{4}'.format(save, int_strs[el[0][0], 4], int_strs[el[1][0], 4], el[1][1], el[1])

				save_2 = '{0}\t{1}\t{2}'.format(int_strs[el[0][0], 3], int_strs[el[0][0], 4], int_strs[el[1][0], 4])
				anty_save_2 = '{0}\t{2}\t{1}'.format(int_strs[el[0][0], 3], int_strs[el[0][0], 4], int_strs[el[1][0], 4])
				#inter_only.write(save_2 + '\n')

				if save_2 not in string_checker and anty_save_2 not in string_checker: 
	
					string_checker.append(save_2)	
					string_checker.append(anty_save_2)
					#clean_inter_only_single.write(save_2 + '\n') #saves only lines which haven't occured before. (anty interactions)
					#check_interaction_type_and_save(int_strs[el[0][0], 3],int_strs[el[0][0], 4], int_strs[el[1][0], 4], save_2)


					removed.write(save_2 + '\n')

					#test_clean = False
				
			elif el[0][1]=='wrong' or el[1][1]=='wrong':
				
				save_2 = '{0}\t{1}\t{2}'.format(int_strs[el[0][0], 3], int_strs[el[0][0], 4], int_strs[el[1][0], 4])
				anty_save_2 = '{0}\t{2}\t{1}'.format(int_strs[el[0][0], 3], int_strs[el[0][0], 4], int_strs[el[1][0], 4])
				#inter_only.write(save_2 + '\n')

				if save_2 not in string_checker and anty_save_2 not in string_checker: 
	
					string_checker.append(save_2)	
					string_checker.append(anty_save_2)
					#clean_inter_only_single.write(save_2 + '\n') #saves only lines which haven't occured before. (anty interactions)
					#check_interaction_type_and_save(int_strs[el[0][0], 3],int_strs[el[0][0], 4], int_strs[el[1][0], 4], save_2)


					removed.write(save_2 + '\n')
				
			
		

			#if test_clean == True: clean.write(save + '\n') 
			
			#full.write(save + '\n')

		
		
	#full.close()
	#clean.close()
	#inter_only.close()
	#clean_inter_only.close()
	#clean_inter_only_single.close()
	promoter_promoter_interactions.close()
	enhancer_enhancer_interactions.close()
	promoter_enhancer_interactions.close()
	removed.close()	
