import itertools


option_to_plot = ["SELECTIVE", "ALL"][0]

combinations = []

stuff = [1, 2, 3, 4]

for L in range(0, len(stuff)+1):
	for subset in itertools.combinations(stuff, L):
		if len(subset): combinations += [list(subset)]

if option_to_plot == "ALL": 
	selected_combinations = combinations
elif option_to_plot == "SELECTIVE":
	selected_combinations = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()


name_of_output_FDR_file = 'FDR_file_{0}_{1}_{4}_{2}_{3}_average_PolII'.format(config_variables.chroms_in_prior[0], config_variables.chroms_to_infer[0], config_variables.one_sided_or_two_sided, config_variables.use_smooth_prior_for_estimation, config_variables.number_of_bins)

name_of_output_FDR_file += "_{0}_{1}_{2}".format(config_variables.upstream, config_variables.downstream, config_variables.upstream_t_s)

if config_variables.disentagled_features_validation: 

	name_of_output_FDR_file += "_TSS" 
else:
	name_of_output_FDR_file += "_GENE"

if option_to_plot == "ALL": name_of_output_FDR_file += "_ALL"

pdf = PdfPages(name_of_output_FDR_file + "_font_size_test.pdf")


if option_to_plot == "ALL":

	f, ax = plt.subplots(1, len(selected_combinations), sharex=True, sharey=True, figsize=(50,10))
	f.subplots_adjust(left=0.035, bottom=0.1, right=0.975, top=0.925, hspace=0.1, wspace=0.1)
	plt.rcParams['xtick.labelsize'] = 28
	plt.rc('ytick', labelsize = 28)
	red_blue_yellow_cyan_marker_size = 12
	red_blue_yellow_cyan_marker_size_legend_box = 18
	legend_box_names_font_size = 22
	size_of_combination_name = 26
	size_of_y_label = 30

	#ax[0].set_title(comb, fontsize = size_of_combination_name)
	ax[0].set_ylabel('Precision', fontsize = size_of_y_label)

	import matplotlib.lines as mlines
	blue_line = mlines.Line2D([], [], color='blue', marker='^', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data + prior')
	yellow_line = mlines.Line2D([], [], color='yellow', marker='s', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='prior')
	red_line = mlines.Line2D([], [], color='red', marker='o', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data')
	pink_line = mlines.Line2D([], [], color='cyan', marker='*', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='MOG')
	handles, labels = ax[0].get_legend_handles_labels()
	ax[0].legend(handles, labels, fontsize = legend_box_names_font_size, scatterpoints = 1)


elif option_to_plot == "SELECTIVE":

	f, ax = plt.subplots(1, len(selected_combinations), sharex=True, sharey=True, figsize=(20,10))
	f.subplots_adjust(left=0.085, bottom=0.15, right=0.965, top=0.925, hspace=0.1, wspace=0.05)
	plt.rcParams['xtick.labelsize'] = 28
	plt.rc('ytick', labelsize = 28)
	red_blue_yellow_cyan_marker_size = 12
	red_blue_yellow_cyan_marker_size_legend_box = 18
	legend_box_names_font_size = 22
	size_of_combination_name = 35
	size_of_y_label = 35

	#ax[0].set_title(comb, fontsize = size_of_combination_name)
	ax[0].set_ylabel('Precision', fontsize = size_of_y_label)

	import matplotlib.lines as mlines
	blue_line = mlines.Line2D([], [], color='blue', marker='^', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data + prior')
	yellow_line = mlines.Line2D([], [], color='yellow', marker='s', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='prior')
	red_line = mlines.Line2D([], [], color='red', marker='o', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='data')
	pink_line = mlines.Line2D([], [], color='cyan', marker='*', markersize = red_blue_yellow_cyan_marker_size_legend_box, label='MOG')
	handles, labels = ax[0].get_legend_handles_labels()
	ax[0].legend(handles, labels, fontsize = legend_box_names_font_size, scatterpoints = 1)


centres_of_ticks = np.arange(4) + 0.5
ind = centres_of_ticks

for index_opt, option_ in enumerate(selected_combinations):
			
	comb = ",".join([dict_option[el] for el in option_])

	if option_ == combinations[-1]: comb = "All"

	OX = [0.1,"0.2\n TPR ",0.3, "\nMAP"]
	ax[index_opt].set_title(comb, fontsize = size_of_combination_name)

	ax[index_opt].vlines(3., 0, 1, colors=u'SlateGray', linestyles=u'dashed')
	ax[index_opt].set_xlim([0., 4.])
	ax[index_opt].set_ylim([0., 1.])

	ax[index_opt].set_xticks(centres_of_ticks)
	ax[index_opt].set_xticklabels(np.array(OX, str))

pdf.savefig()
pdf.close()

plt.show()
