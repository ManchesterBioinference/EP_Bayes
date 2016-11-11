def kern_scipy_gaus(x, color, xgrid, bandwidth = "scott", plot_atr = False, label=None):

	import numpy as np
	from scipy.stats import kde
	import matplotlib.pyplot as plt
 
	density = kde.gaussian_kde(x, bw_method = bandwidth)
	

	if plot_atr:
		if label == None: plt.plot(xgrid, density(xgrid), color + '-', lw=3)
		else: plt.plot(xgrid, density(xgrid), color + '-', lw=3, label = label)
		plt.fill_between(xgrid,density(xgrid), where=None, color=color, alpha = 0.2)

	return density(xgrid), xgrid

def chrom_cross_validation_distance(prior_elements, thresholds, classification_of_interactions, plot_likelihood_function = False, weights = None):
	import config_variables
	import itertools
	
	mode = config_variables.mode
	np = config_variables.np
	chroms_in_prior = config_variables.chroms_in_prior
	domain = config_variables.domain

	
	
	def calculator_distance(factor):
		try: factor = factor[0]
		except: pass
		chrom_elements = prior_elements[mode][classification_of_interactions]["distance"]["attribute_values"]
		
		likelihood_elements = np.zeros(len(chroms_in_prior))
		if domain:
			possible_distances_counts = prior_elements[mode][classification_of_interactions]["distance"]["possible_distances_counts"]	
			possible_distances_counts_mask = possible_distances_counts > 0
			weights = (1./possible_distances_counts[possible_distances_counts_mask])/np.sum(1./possible_distances_counts[possible_distances_counts_mask])


		ends = [len(chrom_elements[chrom__]) - 1 for chrom__ in chroms_in_prior]
		total_x = np.array(list(itertools.chain.from_iterable([chrom_elements[chrom__] for chrom__ in chroms_in_prior])))

		for index, chrom_ in enumerate(chroms_in_prior):

			mask = np.ones(len(total_x), bool)
			mask_chrom = np.ones(len(total_x), bool)
			mask[sum(ends[:index]):sum(ends[:index+1])] = False
			mask_chrom[mask] = False

			if domain: 
				mask = mask[possible_distances_counts_mask]
				mask_chrom = mask_chrom[possible_distances_counts_mask]
				total_x_ = total_x[possible_distances_counts_mask][:]

			total_x_ = total_x

			mask[sign*total_x_ <= 0] = False
			mask_chrom[sign*total_x_ <= 0] = False

			total_x_without_chrom = sign*total_x_[mask]
			chrom_grid = sign*total_x_[mask_chrom]

			if not(domain): f_chrom_, xgrid = kern_scipy_gaus(total_x_without_chrom, "g", chrom_grid, bandwidth = factor, plot_atr = False)
			else: f_chrom_, xgrid = kern_scipy_gaus_weighted(total_x_without_chrom, "g", chrom_grid, weights = weights[mask], bandwidth = factor, factor = None)

			likelihood_elements[index] = np.sum(np.log(f_chrom_))

		return -1*np.sum(likelihood_elements)

	from scipy.optimize import minimize_scalar
	optimum = {}

	
	for sign, positive_or_negative_side in zip([1, -1], ["positive_side", "negative_side"]):

		optimum[positive_or_negative_side] = minimize_scalar(calculator_distance, bounds = (10**-3, 1.), method='bounded', tol = 1e-5).x.reshape(1,)[0]

	#likelihood = {}
	
	#likelihood["distance"] = np.zeros((len(thresholds), 2))


	#for index__, factor in enumerate(thresholds):
	#	likelihood["distance"][index__] = calculator_distance(factor)
	#if plot_likelihood_function:
	#	import matplotlib.pyplot as plt
	#	plt.figure("likelihood")
	#	plt.plot(thresholds, likelihood["distance"][:,0])
	#	plt.plot(thresholds, likelihood["distance"][:,1])
	#	#plt.show()

	#optimum = thresholds[np.argmax(likelihood["distance"], axis = 0)].tolist()
	
	return optimum

def chrom_cross_validation_correlation(prior_elements, data_set_name, thresholds, classification_of_interactions, plot_likelihood_function = False):
	import config_variables
	import itertools
	
	mode = config_variables.mode
	np = config_variables.np
	chroms_in_prior = config_variables.chroms_in_prior
	def calculator(factor):

		try: factor = factor[0]
		except: pass
		likelihood_elements = np.zeros(len(chroms_in_prior))

		ends = [len(chrom_elements[chrom__]) - 1 for chrom__ in chroms_in_prior]
		total_x = np.array(list(itertools.chain.from_iterable([chrom_elements[chrom__] for chrom__ in chroms_in_prior])))

		for index, chrom_ in enumerate(chroms_in_prior):

			mask = np.ones(len(total_x), bool)
			mask[sum(ends[:index]):sum(ends[:index+1])] = False

			total_x_without_chrom = total_x[mask]

			f_chrom_, xgrid = kern_scipy_gaus(total_x_without_chrom, "g", chrom_elements[chrom_], bandwidth = factor, plot_atr = False)
			likelihood_elements[index] = np.sum(np.log(f_chrom_[f_chrom_>0]))

		return -1*np.sum(likelihood_elements)

	#likelihood = {}
	#likelihood[data_set_name] = np.zeros(len(thresholds))
	chrom_elements = prior_elements[mode][classification_of_interactions]["correlation"][data_set_name]["attribute_values"]

	#for index__, factor in enumerate(thresholds):
	#	likelihood[data_set_name][index__] = calculator(factor)

	#optimum = thresholds[np.argmax(likelihood[data_set_name])]
	#if plot_likelihood_function:
	#	import matplotlib.pyplot as plt
	#	plt.figure("likelihood")
	#	plt.plot(thresholds, likelihood[data_set_name])
		#plt.show()

	from scipy.optimize import minimize_scalar
	#kern = kde.gaussian_kde(x, weights=weights, bw_method = "scott")
	optimum = minimize_scalar(calculator, bounds = (10**-3, 1.), method='bounded', tol = 1e-5).x.reshape(1,)[0]



	return optimum

def kern_scipy_gaus_weighted(x, color, xgrid, weights, bandwidth = "scott", factor = None, plot_atr = False):

	import numpy as np
	from scipy.stats import kde
	import matplotlib.pyplot as plt
	import gaussian_kde_weighted_cross_validation as kde
	

	if bandwidth == "cross_validation":
		from scipy.optimize import minimize_scalar
		kern = kde.gaussian_kde(x, weights=weights, bw_method = "scott")
		bandwidth = minimize_scalar(kern.cross_validation, bounds = (10**-3, 1.), method='bounded', tol = 1e-5).x.reshape(1,)[0]
	
	density = kde.gaussian_kde(x, weights=weights, bw_method = bandwidth, bw_factor = factor)	

	if plot_atr:
		plt.plot(xgrid, density(xgrid), color + '-', lw=3)
		plt.fill_between(xgrid,density(xgrid), where=None, color=color, alpha = 0.2)

	return density(xgrid), xgrid

def cross_validation(x, **kwargs):
	import numpy as np
	from sklearn.grid_search import GridSearchCV
	from sklearn.neighbors import KernelDensity
	grid = GridSearchCV(KernelDensity(**kwargs),
                    {'bandwidth': np.linspace(10**-3, 1.0, 300)},
                    cv=20) # 20-fold cross-validation
	grid.fit(x[:, None])
	return grid.best_estimator_.bandwidth


def kern_sklearn_expon(x, color, x_grid, bandwidth = [], kernel_ = "gaussian"):

	
	import numpy as np
	import matplotlib.pyplot as plt
	from sklearn.neighbors import KernelDensity
	bandwidth = cross_validation(x)
	"""Kernel Density Estimation with Scikit-learn"""
	kde_skl = KernelDensity(bandwidth=bandwidth, kernel = kernel_)
	kde_skl.fit(x[:, np.newaxis])
	log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])	
	pdf = np.exp(log_pdf)
	plt.plot(x_grid, pdf, color=color, alpha=0.5, lw=3)
	plt.fill_between(x_grid, pdf, where=None, color=color, alpha = 0.2)
	return pdf, x_grid



def kernel_weighted_samples(x, color, x_grid, **kwargs):
	from statsmodels.nonparametric.kde import KDEUnivariate
	import matplotlib.pyplot as plt
	"""Univariate Kernel Density Estimation with Statsmodels"""
	kde = KDEUnivariate(x)
	kde.fit(**kwargs) # kernel, bw, fft, weights, gridsize, ...]
	pdf = kde.evaluate(x_grid)

	plt.plot(x_grid, pdf, color=color, alpha=0.5, lw=3)
	plt.fill_between(x_grid, pdf, where=None, color=color, alpha = 0.2)

	return pdf, x_grid
  



