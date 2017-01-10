from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

import matplotlib.pyplot as plt
import matplotlib

import config_variables

chroms_in_prior = config_variables.chroms_in_prior
chroms_to_infer = config_variables.chroms_to_infer

pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")

def plot_scikit_lda(X, y, title, normed):

    ax = plt.subplot(111)
    for label,marker,color in zip(range(1,3),('^', 's'),('blue', 'red')):

        plt.hist(x=X[:,0][y == label],
                    color=color,
                    alpha=0.5,
                    label=label_dict[label],
                    normed = normed,
					bins = 15,
					range = (min(X[:,0]), max(X[:,0])))



    plt.xlabel('projections', fontsize = 18)
    plt.ylabel('frequencies', fontsize = 18)

    leg = plt.legend(loc='upper right', fancybox=True)
    leg.get_frame().set_alpha(0.5)
    plt.title(title, y=1.04)

    # hide axis ticks
    plt.tick_params(axis="both", which="both", bottom="off", top="off",  
            labelbottom="on", left="off", right="off", labelleft="on")

    # remove axis spines
    ax.spines["top"].set_visible(False)  
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)    

    plt.grid()
    plt.tight_layout
    pdf.savefig()
    plt.close()


def calculate_single_ROC_best_True_sensitivity(probabilities_true, probabilities_false):

	_True_positives_of_threshold = []
	_False_positives_of_threshold = []

	sorted_prob_true = np.sort(probabilities_true)
	sorted_prob_false = np.sort(probabilities_false)

	sorted_thresholds = np.sort(np.unique(np.r_[probabilities_true, probabilities_false]))
	sorted_thresholds = np.unique(np.r_[sorted_thresholds, np.max(sorted_thresholds)*1.01])

	len_prob_true = len(probabilities_true)
	len_prob_false = len(probabilities_false)

	print 'len prob: ', len_prob_true, len_prob_false

	_True_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_true, sorted_thresholds)[0])

	_False_positives_of_threshold = np.cumsum(np.histogram(sorted_prob_false, sorted_thresholds)[0])

	Precision = np.array(_True_positives_of_threshold, dtype = float)/(np.array(_True_positives_of_threshold, dtype = float) + np.array(_False_positives_of_threshold, dtype = float))
	
	True_positive_Rate = np.array(_True_positives_of_threshold)/float(len_prob_true)
		
	False_positive_Rate = np.array(_False_positives_of_threshold)/float(len_prob_false)

	print 'number of thresholds', len(True_positive_Rate), len(False_positive_Rate)

	return True_positive_Rate, False_positive_Rate, Precision



a = ['promoter_enhancer_interactions', 'enhancer_enhancer_interactions']
b = ['positive_interactions', 'negative_interactions']
c = ['distance', 'correlation']
d = ['attribute_values', 'negative_side', 'prior_frequencies', 'number_in_bin_of_histogram', 'positive_side', 'prior_bins']

label_dict = {1: 'Positive', 2: 'Negative'}

score = []


#.plot(ind, probabilities_dist, alpha=1.0, color="darkviolet", marker= "s", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "distance")
#.plot(ind, probabilities_correl, alpha=1.0, color="red", marker= "o", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "data")
#.plot(ind, probabilities_dist_correl_MOG, alpha=1.0, color="cyan", marker= "*", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "LVA data+prior")
#.plot(ind, probabilities_dist_MOG, alpha=1.0, color="green", marker= "v", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "LVA prior")
#.plot(ind, probabilities_dist_correl, alpha=1.0, color="b", marker= "^", linewidth=0.0, markersize=red_blue_yellow_cyan_marker_size, label = "NB")

X = {}
y = {}

for selected_combination in selected_combinations:

	plot_label = ",".join([dict_option[el] for el in selected_combination])

	X[plot_label] = {}

	# LDA

	#X = df[[0,1,2,3]].values
	#y = df['class label'].values



	for lab, chroms_el, elements in zip(["training", "test"], [chroms_in_prior, chroms_to_infer], [prior_elements, infered_elements]):

		for feature in ["distance", "correlation"]:

			if feature == "distance":

				total_array_positive = [elements[a[0]][b[0]][c[0]][d[0]][chrom_] for chrom_ in chroms_el]
				total_array_positive = list(itertools.chain.from_iterable(total_array_positive))

				total_array_negative = [elements[a[0]][b[1]][c[0]][d[0]][chrom_] for chrom_ in chroms_el]
				total_array_negative = list(itertools.chain.from_iterable(total_array_negative))

				X[plot_label][lab] = np.zeros((len(total_array_positive) + len(total_array_negative), len(selected_combination)+1))
				X[plot_label][lab][:, 0] = np.r_[total_array_positive, total_array_negative]

				y[lab] = np.zeros(len(X[plot_label][lab])).astype(int)
				y[lab][:len(total_array_positive)] = 1
				y[lab][len(total_array_positive):] = 2
		
			elif feature == "correlation":

				for i, data_set_name in zip(range(1,len(selected_combination)+1), initiate_time_series.datasets_names[selected_combination]):

					total_array_positive = [elements[a[0]][b[0]][c[1]][data_set_name][d[0]][chrom_] for chrom_ in chroms_el]
					total_array_positive = list(itertools.chain.from_iterable(total_array_positive))

					total_array_negative = [elements[a[0]][b[1]][c[1]][data_set_name][d[0]][chrom_] for chrom_ in chroms_el]
					total_array_negative = list(itertools.chain.from_iterable(total_array_negative))

					X[plot_label][lab][:, i] = np.r_[total_array_positive, total_array_negative]


labs = ["training", "test"]
for lab in labs:

	X_trans = {}

	for selected_combination in selected_combinations:

		plot_label = ",".join(initiate_time_series.datasets_names[selected_combination])

		sklearn_lda = LDA(n_components=1)
		#X_lda_sklearn = sklearn_lda.fit(X, y)
	
		sklearn_lda.fit(X[plot_label]["training"], y["training"])
		
		X_trans = sklearn_lda.transform(X[plot_label][lab])
		
		#plot_step_lda() 
		plot_scikit_lda(X_trans, y[lab], title='{0} set, distance, {1} features'.format(lab, ", ".join(initiate_time_series.datasets_names[selected_combination])), normed = True)
		plot_scikit_lda(X_trans, y[lab], title='{0} set, distance, {1} features'.format(lab, ", ".join(initiate_time_series.datasets_names[selected_combination])), normed = False)
		score += [sklearn_lda.score(X[plot_label][lab], y[lab])]
		print sklearn_lda.predict(X[plot_label][lab]) == y[lab]
		print sklearn_lda.predict_proba(X[plot_label][lab])
		print "TRP", sum(sklearn_lda.predict(X[plot_label][lab]) == 1)/float(sklearn_lda.predict(X[plot_label][lab]).shape[0])

		#leg = plt.legend(loc='upper right', fancybox=True)
		#leg.get_frame().set_alpha(0.5)

		#pdf.savefig()
		#plt.close()

	for selected_combination in selected_combinations:

		plot_label = ",".join(initiate_time_series.datasets_names[selected_combination])

		sklearn_lda.fit(X[plot_label]["training"], y["training"])
		
		X_trans = sklearn_lda.transform(X[plot_label][lab])

		#labs = ["training", "test"]
		#lab = labs[0]
		True_positive_Rate, False_positive_Rate, Precision = calculate_single_ROC_best_True_sensitivity(X_trans[y[lab]==1] , X_trans[y[lab]==2])
		plt.plot(True_positive_Rate, Precision, label = "distance," + plot_label)
		plt.title('PR curves for {0} set'.format(lab), y=1.04, fontsize = 18)
		plt.xlabel('TPR (Recall)', fontsize = 18)
		plt.ylabel('Precision', fontsize = 18)

	leg = plt.legend(loc='upper right', fancybox=True)
	leg.get_frame().set_alpha(0.5)
	plt.axis([0,1,0,0.12])

	pdf.savefig()
	plt.close()
		

pdf.close()
#print sklearn_lda.predict(X[plot_label][lab])

