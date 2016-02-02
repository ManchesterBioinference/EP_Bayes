def smooth_priors(ll, true_dist_of_total_prob, true_dist_of_total_bins, false_dist_of_total_prob, false_dist_of_total_bins, kde_mode, bandwidth = "silverman", kernel = "exponential"):

	import matplotlib.pyplot as plt
	import numpy as np

	miss_mass_true = np.diff(true_dist_of_total_bins)[[0,-1]]*true_dist_of_total_prob[[0,-1]]
	miss_mass_false = np.diff(false_dist_of_total_bins)[[0,-1]]*false_dist_of_total_prob[[0,-1]]

	import hist_sampler

	x_ = hist_sampler.histogram_sampler(true_dist_of_total_prob[1:-1], true_dist_of_total_bins[1:-1], 800000)
	y_ = hist_sampler.histogram_sampler(false_dist_of_total_prob[1:-1], false_dist_of_total_bins[1:-1], 800000)

	prob_x_2 ,bins_x_2, patch = plt.hist(abs(x_), 200, normed = True, alpha = 0.15)
	prob_x_3 ,bins_x_3, patch = plt.hist(abs(y_), bins = 600, normed = True, alpha = 0.15)

	plt.clf()

	import kern_density_est
	
	if kde_mode == "scipy_kde": true_dist_of_total_prob_smooth, true_dist_of_total_bins_smooth = kern_density_est.kern_scipy_gaus(abs(x_[abs(x_) < true_dist_of_total_bins[-2]]), 'g', np.linspace(0, true_dist_of_total_bins[-2]+50000, 1000), bandwidth) # can change the kernel to an expenential kern_density_est.kern_sklearn_expon(x, color, bandwidth, kernel = "exponential", **kwargs)

	elif kde_mode == "sklearn_kde": true_dist_of_total_prob_smooth, true_dist_of_total_bins_smooth =  kern_density_est.kern_sklearn_expon(abs(x_), 'r', np.linspace(0, true_dist_of_total_bins[-2]+50000, 3000), bandwidth, kernel_ = kernel)


	false_dist_of_total_prob_smooth, false_dist_of_total_bins_smooth = kern_density_est.kern_sklearn_expon(abs(y_), 'y', np.linspace(0, false_dist_of_total_bins[-2]+50000, 2000), bandwidth = bandwidth)

	plt.clf()

	plt.plot(true_dist_of_total_bins_smooth, true_dist_of_total_prob_smooth*(1-sum(miss_mass_true)), lw = 2,  color = 'g')
	plt.clf()	
	

	plt.plot(false_dist_of_total_bins_smooth, false_dist_of_total_prob_smooth*(1-sum(miss_mass_false)), lw = 2,  color = 'y')	

	plt.clf()	
	#plt.show()

	from scipy.stats import kde
	import bisect

	
	insert_index_2 = bisect.bisect_left(true_dist_of_total_prob_smooth[::-1]-true_dist_of_total_prob[-1], 0.)	

	true_dist_of_total_prob_smooth, true_dist_of_total_bins_smooth = true_dist_of_total_prob_smooth[0:-1*insert_index_2-1], true_dist_of_total_bins_smooth[0:-1*insert_index_2-1]


	insert_index_4 = bisect.bisect_left(false_dist_of_total_prob_smooth[::-1]-false_dist_of_total_prob[-1], 0.)

	false_dist_of_total_prob_smooth, false_dist_of_total_bins_smooth = false_dist_of_total_prob_smooth[0:-1*insert_index_4-1], false_dist_of_total_bins_smooth[0:-1*insert_index_4-1]

	prob_true_ = true_dist_of_total_prob_smooth*(1-np.sum(miss_mass_true))
	prob_false_ = false_dist_of_total_prob_smooth*(1-np.sum(miss_mass_false))
			
	bins_t_plot = np.r_[0, true_dist_of_total_bins_smooth, true_dist_of_total_bins_smooth[-1], true_dist_of_total_bins[-1]]
	prob_t_plot = np.r_[0, prob_true_, true_dist_of_total_prob[-1], true_dist_of_total_prob[-1]]

	bins_f_plot = np.r_[0, false_dist_of_total_bins_smooth, false_dist_of_total_bins_smooth[-1], false_dist_of_total_bins[-1]]
	prob_f_plot = np.r_[0, prob_false_, false_dist_of_total_prob[-1], false_dist_of_total_prob[-1]]


	plt.clf()
	plt.bar(bins_x_2[:-1], prob_x_2*(1-sum(miss_mass_true)), np.diff(bins_x_2), alpha=0.2, color = "green")
	plt.bar(bins_x_3[:-1], prob_x_3*(1-sum(miss_mass_false)), np.diff(bins_x_3), alpha=0.2, color = "yellow")
	plt.axis([bins_f_plot[1], bins_f_plot[-2]+100000, 0., max(prob_t_plot)*1.5])
	ll.savefig()
	plt.clf()

	plt.plot(bins_t_plot, prob_t_plot, 'g')
	plt.plot(bins_f_plot, prob_f_plot, 'y')	

	plt.fill_between(bins_t_plot, prob_t_plot, where=None, color='g', alpha = 0.2, lw = 2)
	plt.fill_between(bins_f_plot, prob_f_plot, where=None, color='y', alpha = 0.2, lw = 2)
	plt.axis([bins_f_plot[1], bins_f_plot[-2]+100000, 0., max(prob_t_plot)*1.5])

	ll.savefig()

	mean_diff_true = np.diff(true_dist_of_total_prob_smooth)/2.
	true_smoothed_mean_probabilities = true_dist_of_total_prob_smooth[:-1] + mean_diff_true

	mean_diff_false = np.diff(false_dist_of_total_prob_smooth)/2.
	false_smoothed_mean_probabilities = false_dist_of_total_prob_smooth[:-1] + mean_diff_false

	
	miss_mass_true_ = [true_dist_of_total_prob[-1]*abs(true_dist_of_total_bins[-1] - true_dist_of_total_bins_smooth[-1])]
	miss_mass_false_ = [false_dist_of_total_prob[-1]*abs(false_dist_of_total_bins[-1] - false_dist_of_total_bins_smooth[-1])]

	prob_true_norm = true_smoothed_mean_probabilities/np.sum(true_smoothed_mean_probabilities*np.diff(true_dist_of_total_bins_smooth))*(1.-np.sum(miss_mass_true_))
	prob_false_norm = false_smoothed_mean_probabilities/np.sum(false_smoothed_mean_probabilities*np.diff(false_dist_of_total_bins_smooth))*(1.-np.sum(miss_mass_false_))

	prob_true = np.r_[prob_true_norm, true_dist_of_total_prob[-1]]
	bins_true = np.r_[true_dist_of_total_bins_smooth, true_dist_of_total_bins[-1]]
	
	prob_false = np.r_[prob_false_norm, false_dist_of_total_prob[-1]]
	bins_false = np.r_[false_dist_of_total_bins_smooth, false_dist_of_total_bins[-1]]

	plt.clf()

	return prob_true, bins_true, prob_false, bins_false, ll


