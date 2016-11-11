def inter_enhancer(chrom):

	negative_interactions = config_variables.negative_interactions
	from  prepare_interactions_clean import un_string

	indexes_p, indexes_e, total_p, total_e = negative_interactions.initialise_variables(chrom)[2:]

	if config_variables.disentagled_features_validation: 
		chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh_TSS[chrom]
	else:
		chr_interactions_pro_enh = config_variables.chr_interactions_dict_pro_enh[chrom]

	true_inter_pro = un_string(chr_interactions_pro_enh[:, :2]).astype(int)

	i_s_t, j_s_t = true_inter_pro[:,0], true_inter_pro[:,1]

	interacting_enhancers_ = np.unique(j_s_t)

	non_interacting_enhancers_ = np.where(np.invert(np.in1d(indexes_e, j_s_t)))[0] + total_e

	return interacting_enhancers_, non_interacting_enhancers_

proximal_enhancers_mask = config_variables.proximal_enhancers_mask
dict_chrom_pro_survived = config_variables.dict_chrom_pro_survived


def plot_double_sided_dist_smooth_histograms(x, freq, color, label):

	plt.plot(x, freq, color + '-', lw=3, label = label)
	plt.fill_between(x,freq, where=None, color=color, alpha = 0.2)


#pro_chroms, pro_coords, pro_time_series = dataset_time_series_dict[link_data_set_name_to_file_name["promoters"]["ER"]]

from matplotlib.backends.backend_pdf import PdfPages
results_folder = config_variables.results_folder

name_of_output_FDR_file = results_folder + 'interacting_vs_non_interacting'

name_of_output_FDR_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

pdf = PdfPages(name_of_output_FDR_file + ".pdf")

#plt.rcParams['xtick.labelsize'] = 15
#plt.rc('ytick', labelsize = 15)
#size_of_combination_name = 23
#size_of_y_label = 20

#plt.set_ylabel('enrichment', fontsize = size_of_y_label)


import kern_density_est as kern

for n, data_set_name in enumerate(datasets_names):
	enh_chroms, enh_coords, enh_time_series = dataset_time_series_dict[link_data_set_name_to_file_name["enhancers"][data_set_name]]

	total_interacting, total_non_interacting =  [], []

	for chrom in chrom_names[:-1]:

		filtered_enhancers = enh_time_series.sum(1) >= 30.

		interacting_enhancers_, non_interacting_enhancers_ = inter_enhancer(chrom)

		chrom_enh_survived_non_interacting = np.where(((enh_chroms == chrom)*np.invert(proximal_enhancers_mask)*filtered_enhancers)[non_interacting_enhancers_])[0]#

		chrom_enh_survived_interacting = np.where(((enh_chroms == chrom)*np.invert(proximal_enhancers_mask)*filtered_enhancers)[interacting_enhancers_])[0]

		total_interacting += enh_time_series[chrom_enh_survived_interacting].sum(1).tolist()
		total_non_interacting += enh_time_series[chrom_enh_survived_non_interacting].sum(1).tolist()


	for x, label, colour in [[np.log10(total_interacting), 'known interacting', "g"], [np.log10(total_non_interacting), 'unknown status', "r"]]:

		xgrid = np.linspace(min(x)*0.8, max(x)*1.2, num=200)

		freq__, bins__ = kern.kern_scipy_gaus(x, colour, xgrid, bandwidth = "scott", plot_atr = False, label=None)

		#print len(freq__), len(bins__)
		plt.figure(n+1, figsize=(8, 6), dpi=200)			
		plot_double_sided_dist_smooth_histograms(bins__[:-1] + np.diff(bins__)/2, freq__[:-1], colour, label)


	if data_set_name == "ER": plt.title(u'ER-\u03B1', fontsize=28)
	else: plt.title(data_set_name, fontsize=28)
	#plt.ylabel('density', fontsize=28)
	#plt.xlabel('correlation', fontsize=28)



	#plt.figure(n, figsize=(9,6))
	plt.xlabel("enrichment", fontsize = size_of_y_label)
	plt.ylabel("density", fontsize = size_of_y_label)
	#plt.title(data_name, fontsize = size_of_y_label)

	plt.rc('ytick', labelsize = 20)
	plt.rc('xtick', labelsize = 20)

	ax = plt.gca(); cur_ylim = ax.get_ylim(); ax.set_ylim([0, cur_ylim[1]])
	plt.subplots_adjust(left=0.12, right=0.95, top=0.925, bottom=0.135)
	x1,x2,y1,y2 = plt.axis(); plt.axis((x1,x2,0.,y2*1.2)); 

	tick_labels = np.arange(0.,x2, 1.0)
	string_labels = [r"$10^{%2d}$" % (i) for i in tick_labels]
	plt.xticks(tick_labels,  string_labels, fontsize=23)#["a", "b", "c", "d", "e"])#
	#plt.xlim([-8.5, 8.5])
	plt.legend(); pdf.savefig();

	#plt.hist(np.log(total_interacting), bins = 100, alpha = 0.3, normed = True, color = "g", label='known interacting')
	#plt.hist(np.log(total_non_interacting), bins = 100, alpha = 0.3, normed = True, color = "r", label='unknown status')

	#plt.legend()
	#pdf.savefig()
	
pdf.close(); plt.close("all")	


