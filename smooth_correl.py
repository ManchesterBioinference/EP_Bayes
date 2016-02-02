def smooth_correl(correl_of_true_inter_total, correl_of_false_inter_total):
	import kern_density_est
	import numpy as np

	correl_of_true_inter_total_prob_smooth, correl_of_true_inter_total_bins_smooth = kern_density_est.kern_scipy_gaus(correl_of_true_inter_total, 'g', np.linspace(-1., 1., 2000))
	correl_of_false_inter_total_prob_smooth, correl_of_false_inter_total_bins_smooth = kern_density_est.kern_scipy_gaus(np.random.choice(correl_of_false_inter_total, int(len(correl_of_false_inter_total)/50.), replace = False), 'y', np.linspace(-1., 1., 2000))


	mean_diff_true = np.diff(correl_of_true_inter_total_prob_smooth)/2.
	correl_of_true_inter_total_prob_smooth = correl_of_true_inter_total_prob_smooth[:-1] + mean_diff_true
	mean_diff_false = np.diff(correl_of_false_inter_total_prob_smooth)/2.
	correl_of_false_inter_total_prob_smooth = correl_of_false_inter_total_prob_smooth[:-1] + mean_diff_false

	extend_mass_true_first = correl_of_true_inter_total_prob_smooth[0]*np.diff(correl_of_true_inter_total_bins_smooth)[0] 
	extend_mass_true_last = correl_of_true_inter_total_prob_smooth[-1]*np.diff(correl_of_true_inter_total_bins_smooth)[-1]
	correl_of_true_inter_total_bins_smooth[0] = -1.00001
	correl_of_true_inter_total_bins_smooth[-1] = 1.00001
	correl_of_true_inter_total_prob_smooth[0] = extend_mass_true_first/np.diff(correl_of_true_inter_total_bins_smooth)[0] 
	correl_of_true_inter_total_prob_smooth[-1] = extend_mass_true_last/np.diff(correl_of_true_inter_total_bins_smooth)[-1]

	extend_mass_false_first = correl_of_false_inter_total_prob_smooth[0]*np.diff(correl_of_false_inter_total_bins_smooth)[0] 
	extend_mass_false_last = correl_of_false_inter_total_prob_smooth[-1]*np.diff(correl_of_false_inter_total_bins_smooth)[-1]
	correl_of_false_inter_total_bins_smooth[0] = -1.00001
	correl_of_false_inter_total_bins_smooth[-1] = 1.00001
	correl_of_false_inter_total_prob_smooth[0] = extend_mass_false_first/np.diff(correl_of_false_inter_total_bins_smooth)[0] 
	correl_of_false_inter_total_prob_smooth[-1] = extend_mass_false_last/np.diff(correl_of_false_inter_total_bins_smooth)[-1]				
	
	return correl_of_true_inter_total_prob_smooth, correl_of_true_inter_total_bins_smooth, correl_of_false_inter_total_prob_smooth, correl_of_false_inter_total_bins_smooth
