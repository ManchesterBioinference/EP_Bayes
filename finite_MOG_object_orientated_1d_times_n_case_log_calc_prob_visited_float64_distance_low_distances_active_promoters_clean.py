import numpy as np
#import pylab as pb
from scipy import stats
#from math import *
import copy
#from operator import itemgetter			
#from sys import argv
#import bisect as bis
import time
import os

#import threading
import config_variables
from numpy.random import random_sample
import numexpr as ne

#mode = config_variables.mode
dict_chrom_distant = config_variables.dict_chrom_distant
proximal_enhancers_mask = config_variables.proximal_enhancers_mask
dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived
dict_chrom_enh_survived = config_variables.dict_chrom_enh_survived
filtered_promoters = config_variables.filtered_promoters
filtered_enhancers = config_variables.filtered_enhancers
probabilities_of_a_bin = config_variables.probabilities_of_a_bin
adequate_histogram_bins = config_variables.adequate_histogram_bins
dataset_time_series_dict = config_variables.dataset_time_series_dict
link_data_set_name_to_file_name = config_variables.link_data_set_name_to_file_name
datasets_names = config_variables.datasets_names
log_distances = config_variables.log_distances
test_prior = config_variables.test_prior
interacting_enhancers_only_MOG = config_variables.interacting_enhancers_only_MOG



class finite_mix:
	def __init__(self, probabilities_of_a_bin, adequate_histogram_bins, active_promoters, active_promoters_time_series, peaks, X, chrom_enh_survived, enh_coords, chrom, kappa_0, mu_0, alpha_0, Beta_0, output, output_likelihood, mode_of_sampler): #output_mean, output_std, mode_of_sampler):
		
		self.mode_of_sampler = mode_of_sampler
		self.chrom = chrom
		#self.emp_distances = emp_distances
		self.output = output
		

		self.output_likelihood = output_likelihood
		#self.output_mean = output_mean
		#self.output_std = output_std

		self.chrom_enh_survived = chrom_enh_survived
		self.enh_coords = enh_coords	
	
		self.peaks = peaks
		
		self.X = X # same size as peaks

		self.D = np.shape(self.X)[1]

		#self.X = self.data_mean_std_normalisation(self.X)
		
		self.active_promoters = active_promoters# independent of peaks

		self.active_promoters_time_series = active_promoters_time_series			
	
		#self.active_promoters_time_series = self.data_mean_std_normalisation(self.active_promoters_time_series)
	
		#self.string = '{0}_{1}_{2}_{3}_{4}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom)

		#prior parameters
		#self._c_trace = [] # stores consequtive c outcomes.
		#self._mean_trace = []
		#self._standard_dev_trace=[]
		
		self.kappa_0 = float(kappa_0) 
		self.mu_0 = float(mu_0) 
		self.alpha_0 = float(alpha_0)
		self.Beta_0 = 1./float(Beta_0) # converts scale to rate, parameters should be given in scale. Derivations are in rate
		#self.mu_0 = np.zeros(self.D)

		if mode_of_sampler == "distance_MOG_empir_mu":
			self.mu_0 = self.active_promoters_time_series
			
		self.time_series = np.r_[self.active_promoters_time_series, self.X]
		
		self.K = len(self.active_promoters)
		self.J = len(self.peaks)
		self.N = self.K + self.J
				
		self.gaussian_const = -0.5*float(self.D)*np.log(2.0*np.pi)

		if test_prior: self.total_dist = []

		self.vec_gamma = np.vectorize(np.random.gamma)
		self.vec_normal = np.vectorize(np.random.normal)

		if mode_of_sampler == "dirichlet_MOG":
			self.alpha = 1.
			self.cl = self.K # ammend that if you want less clusters "dirichlet_MOG"
			self.labels = np.arange(self.cl)
			self.c = self.labels[np.random.randint(0, self.cl, self.N)]

			self.initialise_probabilities_and_lookup_arrays()
		


		else:

			self.dist = self.calculate_distances_between_promoters_and_non_promoters()
			self.empirical_prob_inter_profile = probabilities_of_a_bin #"distance_MOG"
			self.bins = adequate_histogram_bins #"distance_MOG"

			neg_dist = self.dist < 0

			if log_distances:
				negative_distances = -1*self.dist[neg_dist]
				positive_distances = self.dist[np.invert(neg_dist)]
				negative_distances = -1*np.log10(negative_distances)
				positive_distances = np.log10(positive_distances)
			else: 
				negative_distances = self.dist[neg_dist]
				positive_distances = self.dist[np.invert(neg_dist)]

			self.dist[neg_dist] = negative_distances
			self.dist[np.invert(neg_dist)] = positive_distances 

			distance_prob_component_neg = self.allocate_to_empirical_probabilities_in_histogram(negative_distances, adequate_histogram_bins, probabilities_of_a_bin)[1]
			distance_prob_component_pos = self.allocate_to_empirical_probabilities_in_histogram(positive_distances, adequate_histogram_bins, probabilities_of_a_bin)[1]

			self.distance_prob_component = np.zeros_like(self.dist)
			self.distance_prob_component[neg_dist] = distance_prob_component_neg
			self.distance_prob_component[np.invert(neg_dist)] = distance_prob_component_pos

			self.cl = self.K
			self.labels = np.arange(self.cl)
			self.c = self.labels[np.random.randint(0, self.K, self.J)] # generates randint(1,K) times self.N. Randint generates from as if from range(1,21) = list <1,20>
			self.c = np.r_[np.arange(self.K), self.c]

			if interacting_enhancers_only_MOG:

				interacting_enhancers = self.inter_enhancer(chrom)
				self.interacting_enhancers = interacting_enhancers

				self.c = np.zeros(self.J, int)-1
				self.c[interacting_enhancers] = self.labels[np.random.randint(0, self.K, len(interacting_enhancers))] # generates randint(1,K) times self.N. Randint generates from as if from range(1,21) = list <1,20>
				self.c = np.r_[np.arange(self.K), self.c]

				

		#self.matching_label_matrix = np.zeros((self.cl, self.N))
		if mode_of_sampler <> "distance_prior":
			self.standard_dev=np.ones((self.cl, self.D))
			self.mean=np.zeros((self.cl, self.D))
			self.variance = self.standard_dev**2.0
	

# Iterator function
	def sample(self,iters):
		"""
		does the actual job to produce final results
		"""
		t_new = time.time()
			
		for k in range(iters):
			t_iter = time.time()
			#self._c_trace.append(copy.deepcopy(self.c[self.K:]))
			#self._mean_trace.append(copy.deepcopy(self.mean))
			#self._standard_dev_trace.append(copy.deepcopy(self.standard_dev))
			#self._alpha_trace.append(self.alpha)
			self.update_c_new()
			#print self.c
			if test_prior: 
				distances_of_interacting_enhancers = self.test_priors_funct(self.dist[np.arange(self.J), self.c[self.K:]])
				self.total_dist += distances_of_interacting_enhancers.tolist()
			self.save_to_to_file_by_line()
			#t_new = time.time()
			
			if self.mode_of_sampler <> "distance_prior": 
				self.calculate_likelihood_of_X_max_ne()
				self.update_mean_and_standard_deviation_2()
				
				#if k%10==0:self.save_to_to_file_by_line_mean_and_std()
			#print "update_std =", time.time() - t_new
			#self.update_alpha() - we fix it for a time being since it has to be sampled from adaptive rejection sampling
			#self.update_beta()    
	
			if k%1==0: 
				t_old = t_new
				t_new = time.time()				
				print k, t_new - t_old, t_new - t_iter
			print k, t_new - t_old, t_new - t_iter

	def test_priors_funct(self, distances_per_sample):
		from prepare_interactions_clean import un_string
		chrom = self.chrom
		negative_interactions = config_variables.negative_interactions
		indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

		if config_variables.disentagled_features_validation: 
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
		else:
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

		true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

		i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]
		interacting_enhancers = np.unique(j_s_t)-total_e

		take_distances_array = np.zeros(len(indexes_e))
		take_distances_array[self.chrom_enh_survived-total_e] = distances_per_sample
		distaces_of_interacting_enhancers = take_distances_array[interacting_enhancers]
		return distaces_of_interacting_enhancers
		#interacting_enhancers_mask = np.zeros(indexes_e).astype(bool)
		#interacting_enhancers_mask[interacting_enhancers] = True


	def inter_enhancer(self, chrom):

		negative_interactions = config_variables.negative_interactions
		from  prepare_interactions_clean import un_string

		indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

		if config_variables.disentagled_features_validation: 
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
		else:
			chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

		true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

		i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

		interacting_enhancers_ = np.unique(j_s_t)-total_e

		enhancer_array_survived = np.zeros(len(indexes_e), bool)
		enhancer_array_interacting = np.zeros(len(indexes_e), bool)

		enhancer_array_survived[self.chrom_enh_survived-total_e] = True
		enhancer_array_interacting[interacting_enhancers_] = True

		mask_interacting_c = np.in1d(np.where(enhancer_array_survived)[0], np.where(enhancer_array_interacting)[0])

		return np.where(mask_interacting_c)[0]


	def save_to_to_file_by_line(self):
		if self.mode_of_sampler == "dirichlet_MOG":
			self.output.write(",".join(copy.deepcopy(self.c).astype(str))+"\n")
		else:
			self.output.write(",".join(copy.deepcopy(self.c[self.K:]).astype(str))+"\n")
		#print ",".join(copy.deepcopy(self.c[self.K:]).astype(str))

	def save_to_to_file_by_line_mean_and_std(self):
		np.savetxt(self.output_mean, self.mean.reshape(1,self.mean.size))
		#self.output_mean.write("\n")

		np.savetxt(self.output_std, self.standard_dev.reshape(1,self.standard_dev.size))


	def data_mean_std_normalisation(self, ts_to_norm):
		means=np.mean(ts_to_norm, axis=1)
		st_dev=np.std(ts_to_norm, axis=1)
		norm_time_series = (ts_to_norm-means[:, None]) / st_dev[:, None]

		return norm_time_series

	def add_noise(self):
		for index in len(self.X):
			self.X[index] += np.random.normal(0.0, 0.3, self.D)	
		
	def initialise_probabilities_and_lookup_arrays(self):

		self.total_n = np.histogram(self.c, bins=self.cl, density=False)[0]	
		self.probabilities = (self.total_n + self.alpha/float(self.cl))/(float(self.N) + self.alpha - 1.0)
	
	
	def calculate_distances_between_promoters_and_non_promoters(self):
		"""
		as the name suggests 
		"""
		from chrom_specific_negative_interactions import extract_TSS_coordinates, calculate_distances

		TSS_coordinates = extract_TSS_coordinates(config_variables.upstream)

		indexes_p = dict_chrom_pro_survived[self.chrom]
		indexes_e = self.chrom_enh_survived

		point_coordinates_promoter, point_coordinates_enhancer = TSS_coordinates[indexes_p], np.mean(self.enh_coords[indexes_e], axis = 1)

		distances_matrix = calculate_distances(False, point_coordinates_promoter, point_coordinates_enhancer)

		dist_from_enhancer_to_promoters = distances_matrix[len(indexes_p):, 0: len(indexes_p)]

		return dist_from_enhancer_to_promoters

		#self.dist=np.array([abs(np.mean(peak) - np.mean(self.active_promoters, axis = 1)) for peak in self.peaks])


	def allocate_to_empirical_probabilities_in_histogram(self, array_to_allocate, adequate_histogram_bins, probabilities_of_a_bin):

		init_shape = array_to_allocate.shape
		
		array_to_allocate = array_to_allocate.reshape(array_to_allocate.size)
		
		indexes_of_bins_it_falls_into = np.digitize(array_to_allocate, adequate_histogram_bins)#, right=False)
		indexes_of_nan = np.isnan(array_to_allocate)
		indexes_of_bins_it_falls_into[indexes_of_nan] = 1

		if np.min(indexes_of_bins_it_falls_into) == 0 or np.max(indexes_of_bins_it_falls_into) == len(adequate_histogram_bins):
			print 'allocation_problem', np.min(indexes_of_bins_it_falls_into) == 0, np.max(indexes_of_bins_it_falls_into) == len(adequate_histogram_bins)
			print np.where(indexes_of_bins_it_falls_into == 0)[0], array_to_allocate[np.where(indexes_of_bins_it_falls_into == 0)[0]]
					
			output_probabilities = probabilities_of_a_bin[indexes_of_bins_it_falls_into - 1]
			output_probabilities[indexes_of_nan] = np.min(probabilities_of_a_bin)*0.9
			
			return indexes_of_bins_it_falls_into.reshape(init_shape), output_probabilities.reshape(init_shape)#, indexes_of_bins_it_falls_into.reshape(init_shape)

		else:
		
			output_probabilities = probabilities_of_a_bin[indexes_of_bins_it_falls_into - 1]
			output_probabilities[indexes_of_nan] = np.min(probabilities_of_a_bin)*0.9
			return indexes_of_bins_it_falls_into.reshape(init_shape), output_probabilities.reshape(init_shape)#, indexes_of_bins_it_falls_into.reshape(init_shape)


	def distance_prob_component_saver(self):
		np.savetxt('distance_prob_component', self.distance_prob_component, delimiter='\t')   
	
	

#	def normal_component_of_X_2(self, j): # for tests only
#		"""
#		calculates the normal component of the probabilities that a non-promotor feature belongs to a promotor.
#		"""
#		
#		self.prob_Y_i = np.array([-0.5*(float(self.D)*np.log(2.0*np.pi) + np.sum(np.log(np.array(self.variance[k], dtype=np.float64))) + np.sum(np.array(self.X[j]-self.mean[k], dtype=np.float64)**2.0/np.array(self.variance[k], dtype=np.float64))) for k in range(self.K)], dtype=np.float64)
#
#		#self.prob_Y_i = np.exp(self.prob_Y_i)#-max(self.prob_Y_i))
#		self.prob_Y_i = np.exp(self.prob_Y_i-max(self.prob_Y_i))	
#		return self.prob_Y_i/sum(self.prob_Y_i)
#	
#	def normal_component_of_X_max(self): # for tests only
#
#		if self.mode_of_sampler == "distance_MOG": X = self.X
#		if self.mode_of_sampler == "dirichlet_MOG": X = self.time_series
#		
#		g = ((X[:, np.newaxis] - self.mean))**2.0/self.variance 
#
#		prob = -0.5*(np.sum(g, axis=2) + np.sum(np.log(self.variance), axis=1))
#		prob_Y_i = np.exp(prob-np.max(prob, axis = 1)[:, np.newaxis] + self.gaussian_const)
#		return prob_Y_i

	def normal_component_of_X_max_ne(self):

		if self.mode_of_sampler == "dirichlet_MOG": 
			X = self.time_series
		else: 
			X = self.X


		X, mean, variance = X[:, np.newaxis], self.mean, self.variance
		expre = "(X - mean)**2/ variance"
		g = ne.evaluate(expre)

		prob = -0.5*(np.sum(g, axis=2) + np.sum(np.log(self.variance), axis=1))
		#prob = -0.5*ne.evaluate("sum(g, 2)") +  ne.evaluate("sum(log(variance), 1)")

		prob_Y_i = np.exp(prob-np.max(prob, axis = 1)[:, np.newaxis] + self.gaussian_const)
		return prob_Y_i


#	def weighted_values(values, probabilities, size):
#		bins = np.cumsum(probabilities, axis (1))
#		return values[np.digitize(random_sample(size), bins)]

	def calculate_likelihood_of_X_max_ne(self):

		mean, variance = self.mean, self.variance

		if self.mode_of_sampler == "dirichlet_MOG": 
			X_ = self.time_series
			c_ = self.c
		else: 
			X_ = self.X
			c_ = self.c[self.K:]

		if interacting_enhancers_only_MOG:
			c_ = c_[self.interacting_enhancers]
			X_ = X_[self.interacting_enhancers]

		likelihood = -0.5*np.sum(np.log(variance[c_]*2.*np.pi) + (X_ - mean[c_])**2./variance[c_])

		np.savetxt(self.output_likelihood, [likelihood], fmt='%1.6e')

		#return likelihood



	def update_c_new(self): # c_i = j| c_-i,alpha) 
		"""
		finds posterior c 
		can be done in pararell
		"""
		#normal = self.normal_component_of_X_unlooped()
		t_new = time.time()

		if self.mode_of_sampler == "distance_prior":

			self.prob = self.distance_prob_component

		elif self.mode_of_sampler == "distance_MOG" or self.mode_of_sampler == "distance_MOG_empir_mu":

			normal = self.normal_component_of_X_max_ne()
			self.prob = normal*self.distance_prob_component		
			#print "calculate_gauss:", time.time() - t_new

		elif self.mode_of_sampler == "dirichlet_MOG":

			normal = self.normal_component_of_X_max_ne()
			self.prob = normal*self.probabilities			
			#print "calculate_gauss:", time.time() - t_new

		norm = np.sum(self.prob, axis=1)[:, np.newaxis]			

		self.prob = self.prob/norm
		t_new = time.time()
		bins = np.cumsum(self.prob, axis = 1)

		if self.mode_of_sampler == "dirichlet_MOG":

			samples = random_sample(self.N)
			self.c = self.labels[(samples[:,None] <= bins).argmax(1)]
			self.initialise_probabilities_and_lookup_arrays()

		else: 
			samples = random_sample(self.J)

			if interacting_enhancers_only_MOG:
				self.c[self.K + self.interacting_enhancers] = (self.labels[(samples[:,None] <= bins).argmax(1)])[self.interacting_enhancers]

			else:
				self.c[self.K:] = self.labels[(samples[:,None] <= bins).argmax(1)]



#	def update_mean_and_standard_deviation(self):
#		"""
#		Update mean and standrad_deviation using a conjugate-pair P(mean|standard_dev)*P(standard_dev)
#		"""
#	
#		self.matching_label_matrix[:] = self.c
#
#		matches = self.matching_label_matrix == np.arange(self.K)[:, None]
#
#		lengths_k = matches.sum(1).astype(float)
#
#		kappas_k=self.kappa_0+lengths_k
#		alphas_k=self.alpha_0+lengths_k/2.0
#
#		#emp_means_k = np.zeros((self.K, self.D)) 
#		#t_new = time.time()
#		#for k in range(self.K): emp_means_k[k] = self.time_series[matches[k]].mean(0)
#
#		emp_means_k = np.dot(matches, self.time_series)/lengths_k[:,None]
#
#		#print "calculate_means:", time.time() - t_new		
#		mus_k = (self.kappa_0 * self.mu_0 + lengths_k[:, None] * emp_means_k)/kappas_k[:, None]
#
#		#Betas_k = np.zeros((self.K, self.D))
#		t_new = time.time()
#
#		sum_squared_term = np.dot(matches, self.time_series**2) - lengths_k[:,None]*emp_means_k**2
#
#		Betas_k = self.Beta_0 + 0.5*sum_squared_term + 0.5*self.kappa_0*(emp_means_k-self.mu_0)**2.0*(lengths_k/kappas_k)[:,None]
#
#		#for k in np.arange(self.K): Betas_k[k] = self.Beta_0 + 0.5*np.sum((self.time_series[matches[k]]-emp_means_k[k])**2.0, axis=0) + 0.5*self.kappa_0*lengths_k[k]*(emp_means_k[k]-self.mu_0)**2.0/kappas_k[k]
#
#		print "calculate_betas:", time.time() - t_new
#
#		t_new = time.time()
#
#		for k in np.arange(self.K):
#
#			for n in range(self.D):	    
#
#				self.variance[k, n]=1.0/np.random.gamma(alphas_k[k], 1.0/Betas_k[k, n], 1)[0]
#				
#				self.mean[k, n] = np.random.normal(mus_k[k, n], self.variance[k, n]**0.5/(kappas_k[k]**0.5), 1)[0]
#
#		self.standard_dev = self.variance**0.5
#		print "update_means_std:", time.time() - t_new


	def update_mean_and_standard_deviation_2(self):
		"""
		Update mean and standrad_deviation using a conjugate-pair P(mean|standard_dev)*P(standard_dev)
		"""
	
		#self.matching_label_matrix[:] = self.c
 
		#matches = self.matching_label_matrix == np.arange(self.K)[:, None]


		if self.mode_of_sampler == "distance_MOG_empir_mu":
			matches = self.c[self.K:].reshape(1, self.c[self.K:].size) == self.labels[:, None]
			time_series = self.X		
		else:
			matches = self.c.reshape(1, self.c.size) == self.labels[:, None] # for dirichlet self.labels = self.K 
			time_series = self.time_series

		lengths_k = matches.sum(1).astype(float)

		kappas_k=self.kappa_0+lengths_k
		alphas_k=self.alpha_0+lengths_k/2.0

		#emp_means_k = np.zeros((self.K, self.D)) 
		#t_new = time.time()
		#for k in range(self.K): emp_means_k[k] = self.time_series[matches[k]].mean(0)

		if self.mode_of_sampler == "dirichlet_MOG" or self.mode_of_sampler == "distance_MOG_empir_mu":
			emp_means_k = np.zeros((self.cl, self.D))
			non_empty = lengths_k.astype(bool)
			emp_means_k[non_empty] = np.dot(matches, time_series)[non_empty]/lengths_k[non_empty][:,None]

		else:
			emp_means_k = np.dot(matches, time_series)/lengths_k[:,None]

	

		#print "calculate_means:", time.time() - t_new		
		mus_k = (self.kappa_0 * self.mu_0 + lengths_k[:, None] * emp_means_k)/kappas_k[:, None]

		#Betas_k = np.zeros((self.K, self.D))
		t_new = time.time()

		#matches, time_series

		sum_squared_term = np.dot(matches, time_series**2) - lengths_k[:,None]*emp_means_k**2  # devide the equation by N*1/N then terms sum_i(bar(x)*x_i) will = 2bar(x)^2 https://www.strchr.com/standard_deviation_in_one_pass

		Betas_k = self.Beta_0 + 0.5*sum_squared_term + 0.5*self.kappa_0*(emp_means_k-self.mu_0)**2.0*(lengths_k/kappas_k)[:,None]

		#sum_squared_term = np.zeros((self.K, self.D))
		#for k in np.arange(self.K): sum_squared_term[k] = 0.5*np.sum((time_series[matches[k]]-emp_means_k[k])**2.0, axis=0)
		#Betas_k = self.Beta_0 + sum_squared_term + 0.5*self.kappa_0*(emp_means_k-self.mu_0)**2.0*(lengths_k/kappas_k)[:,None]

		#for k in np.arange(self.K): Betas_k[k] = self.Beta_0 + 0.5*np.sum((self.time_series[matches[k]]-emp_means_k[k])**2.0, axis=0) + 0.5*self.kappa_0*lengths_k[k]*(emp_means_k[k]-self.mu_0)**2.0/kappas_k[k]

		#print "calculate_betas:", time.time() - t_new

		t_new = time.time()

		self.variance[:] = 1./ self.vec_gamma(alphas_k[:,None], 1.0/Betas_k)
		self.standard_dev[:] = self.variance**0.5
		self.mean[:] = self.vec_normal(mus_k, self.standard_dev/(kappas_k[:,None])**0.5)

		
		#print "update_means_std:", time.time() - t_new




#	def results_saver(self):
#
#		name = 'cluster_trace_of_c_distance_' + self.string
#
#		np.save(name, self._c_trace)
#	
#		name = 'cluster_trace_of_means_distance_' + self.string
#
#		np.save(name, self._mean_trace)
#
#		name = 'cluster_trace_of_st_dev_distance_' + self.string
#
#		np.save(name, self._standard_dev_trace)
#
#	
#	def resume(self):
#	
#		name = 'cluster_trace_of_c_distance_' + self.string
#
#		self._c_trace = np.load(name).tolist()
#
#		name = 'cluster_trace_of_means_distance_' + self.string
#
#		self._st_dev_trace = np.load(name).tolist()
#			
#		name = 'cluster_trace_of_st_dev_distance_' + self.string
#
#		self._mean_trace = np.loadtxt(name).tolist()
#
#		self.mean=self._mean_trace[-1]
#		self.standard_dev=self._standard_dev_trace[-1]
#		self.variance=self.standard_dev**2.0


def executor(arr):

	mode_of_sampler, number_of_samples, option_correl, chrom, chain_number, continue_sampling = arr

	pro_chroms, pro_coords, pro_time_series = dataset_time_series_dict[link_data_set_name_to_file_name["promoters"]["ER"]]
	enh_chroms, enh_coords, enh_time_series = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"]["ER"]]

	chrom_enh_survived = np.where((enh_chroms == chrom)*np.invert(proximal_enhancers_mask)*filtered_enhancers)[0]	

	enhancer_coordinates_chrom_survived = enh_coords[chrom_enh_survived]
	promoter_coordinates_chrom_survived = pro_coords[dict_chrom_pro_survived[chrom]]

	def data_mean_std_normalisation(ts_to_norm):
		means=np.mean(ts_to_norm, axis=1)
		st_dev=np.std(ts_to_norm, axis=1)
		norm_time_series = (ts_to_norm-means[:, None]) / st_dev[:, None]

		return norm_time_series

	concat = {}

	for feature in ["enhancers", "promoters"]:

		concat[feature] = {}	
		for datasets_name in datasets_names[option_correl]:

			ts_to_norm = dataset_time_series_dict[link_data_set_name_to_file_name[feature][datasets_name]][2]

			concat[feature][datasets_name] = data_mean_std_normalisation(ts_to_norm)#data_mean_std_normalisation(ts_to_norm)

		concat[feature]["total"] = np.column_stack([concat[feature][el] for el in datasets_names[option_correl]])

	promoter_time_series_chrom_survived = concat["promoters"]["total"][dict_chrom_pro_survived[chrom]]
	enhancer_time_series_chrom_survived = concat["enhancers"]["total"][chrom_enh_survived]

	#param = "dist"
	save_to_folder = os.getcwd() + "/MOG_results_/"
	if not os.path.exists(save_to_folder): os.makedirs(save_to_folder)

	comb = "_".join([config_variables.dict_option[el] for el in option_correl])
	kappa_0, mu_0, alpha_0, Beta_0 = config_variables.kappa_0, config_variables.mu_0, config_variables.alpha_0, config_variables.Beta_0

	#if mode_of_sampler == "distance_prior": name = 'prior_distance_trace_of_c_{0}_{1}'.format(chrom, number_of_samples)
	#elif mode_of_sampler == "distance_MOG": name = 'MOG_distance_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples)
	#elif mode_of_sampler == "dirichlet_MOG": name = 'MOG_dirichlet_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples)
	#elif mode_of_sampler == "distance_MOG_empir_mu": name = 'MOG_distance_emprirical_mu_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples)

	if mode_of_sampler == "distance_prior":	
		name = 'prior_distance_trace_of_c_{0}_{1}'.format(chrom, number_of_samples)
		output_likelihood = None	
		#name_likelihood = 'prior_distance_trace_of_likelihood_{0}_{1}'.format(chrom, number_of_samples)

	elif mode_of_sampler == "distance_MOG": 
		name = 'MOG_distance_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_std = 'MOG_distance_trace_of_std_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_mean = 'MOG_distance_trace_of_mean_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)

	elif mode_of_sampler == "dirichlet_MOG": 
		name = 'MOG_dirichlet_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_std = 'MOG_dirichlet_trace_of_std_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_mean = 'MOG_dirichlet_trace_of_mean_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)

	elif mode_of_sampler == "distance_MOG_empir_mu": 
		name = 'MOG_distance_emprirical_mu_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_std = 'MOG_distance_emprirical_mu_trace_of_std_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_mean = 'MOG_distance_emprirical_mu_trace_of_mean_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		name_likelihood = 'MOG_distance_emprirical_mu_trace_of_likelihood_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)

	if continue_sampling:

		output = open(save_to_folder + name, 'r')
		for line in output: pass
		last = map(int, line[:-1].split(','))
		output.close()

		if mode_of_sampler <> "distance_prior":
			name_likelihood__ = save_to_folder + name_likelihood
			likelihoods = np.loadtxt(name_likelihood__, dtype=float)[:-1]
			np.savetxt(name_likelihood__, likelihoods, fmt='%1.6e')

		name_2 = save_to_folder + 'MOG_distance_emprirical_mu_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, 2*number_of_samples, chain_number)
		name_likelihood_2 = save_to_folder + 'MOG_distance_emprirical_mu_trace_of_likelihood_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, 2*number_of_samples, chain_number) 

		output = open(name_2, 'w')
		if mode_of_sampler <> "distance_prior": output_likelihood = open(name_likelihood_2, 'w')

	else:	

		output = open(save_to_folder + name, 'w')
		if mode_of_sampler <> "distance_prior": output_likelihood = open(save_to_folder + name_likelihood, 'w')

		
	#output_mean = open(save_to_folder + name_mean, 'w')
	#output_std = open(save_to_folder + name_std, 'w')

	dog = finite_mix(probabilities_of_a_bin, adequate_histogram_bins, promoter_coordinates_chrom_survived, promoter_time_series_chrom_survived, enhancer_coordinates_chrom_survived, enhancer_time_series_chrom_survived,  chrom_enh_survived, enh_coords, chrom, kappa_0, mu_0, alpha_0, Beta_0, output, output_likelihood, mode_of_sampler)#output_mean, output_std, 

	if continue_sampling:
		dog.c[dog.K:] = last
		dog.update_mean_and_standard_deviation_2()

		if mode_of_sampler <> "distance_prior": dog.calculate_likelihood_of_X_max_ne()

	dog.sample(number_of_samples)


	#dog.results_saver()


	output.close()
	if mode_of_sampler <> "distance_prior": output_likelihood.close()

	#output_mean.close()
	#output_std.close()

	#if continue_sampling:
		#name = save_to_folder + 'MOG_distance_emprirical_mu_trace_of_c_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)
		#name_likelihood = save_to_folder + 'MOG_distance_emprirical_mu_trace_of_likelihood_{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(kappa_0, mu_0, alpha_0, Beta_0, chrom, comb, number_of_samples, chain_number)

		#os.system("mv {0} {1}".format(name, name_2))
		#os.system("mv {0} {1}".format(name_likelihood, name_likelihood_2))

	if test_prior:

		return dog.total_dist[1:]
	#dog.results_saver()
	#return 1

