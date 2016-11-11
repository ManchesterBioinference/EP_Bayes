import os
import shlex, subprocess
import numpy as np
import re
import config_variables
from matplotlib.backends.backend_pdf import PdfPages
from pandas import *


def executor(merged_time_series_to_cluster, upstream_TSS = 0, downstream_TSS = 0, diff_bind_version = False, mode_atr = ["FIRST_TS", "SECOND_TS"][1], mode_atr2 = ["ENHANCER", "GENE", "TSS"][1], GLOBAL_OR_SURVIVED = ["survived", "global"][1], mode_of_data_sets = ["Ciiras", "Others_from_cistrom_finder"][0], sorted_mode = ["amplitude_sorted", "size_sorted"][1], dont_plot = ["ESR2", "RAD21"]):

	pwd = os.getcwd()

	hg = 'hg19'

	if mode_of_data_sets == "Ciiras": name_of_files = np.loadtxt(pwd + "/" + hg + "/list_of_files_Ciira_names_changed.txt", dtype = str) 

	elif mode_of_data_sets == "Others_from_cistrom_finder": name_of_files = np.loadtxt(pwd + "/" + hg + "/list_of_files.txt", dtype = str)

	path_to_R = pwd + "/R_scripts/"

	survived = np.loadtxt(path_to_R + '{0}_survived_indexes'.format(merged_time_series_to_cluster)).astype(int) # saved during filtering


	if diff_bind_version:
		peaks = np.loadtxt(config_variables.name_of_enhancer_file_for_overlap, dtype = str)
		indexes_of_DB_peaks = np.loadtxt(pwd + "/" + hg + "/indexes_of_DB_peaks.csv", dtype = int, skiprows = 1, usecols = (1,), delimiter = ",")
		labels = np.zeros(len(peaks), int)
		labels[indexes_of_DB_peaks] = 1
		labels = labels + 1
		labels = labels[survived]

	else:
		labels = np.loadtxt(path_to_R + '{0}_labels'.format(merged_time_series_to_cluster), str, delimiter = ",")[1:,1].astype(int) # from EP clustering

	save_to_temp_folder = pwd + "/" + hg + "/" + merged_time_series_to_cluster + "_results_temp_{0}/".format(GLOBAL_OR_SURVIVED)

	if not os.path.exists(save_to_temp_folder):
		os.makedirs(save_to_temp_folder)


	def create_GENE_file(upstream_TSS):

		data = np.loadtxt(config_variables.name_of_time_series_promoter_file_for_TSS_start, dtype = str,  delimiter = '\t')	
		plus_strand = data[:, 4] == '+'

		data[plus_strand, 1] = data[plus_strand, 1].astype(int) - upstream_TSS
		data[np.invert(plus_strand), 2] = data[np.invert(plus_strand), 2].astype(int) + upstream_TSS

		data = np.c_[data, range(len(data))]

		np.savetxt(name_of_file_for_overlap, data, fmt = '%s', delimiter = '\t')

	def create_TSS_file(upstream_TSS, downstream_TSS):

		data = np.loadtxt(config_variables.name_of_time_series_promoter_file_for_TSS_start, dtype = str,  delimiter = '\t')	
		plus_strand = data[:, 4] == '+'

		data[plus_strand, 2] = data[plus_strand, 1].astype(int) + 1
		data[np.invert(plus_strand), 1] = data[np.invert(plus_strand), 2].astype(int) - 1

		data = np.c_[data, range(len(data))]

		np.savetxt(name_of_file_for_overlap, data, fmt = '%s', delimiter = '\t')

	if mode_atr2 == "TSS":	
		name_of_file_for_overlap = save_to_temp_folder + config_variables.name_of_time_series_promoter_file_for_TSS_start[7:-3] + "_TSS_{0}_{1}".format(upstream_TSS, downstream_TSS)
		create_TSS_file(upstream_TSS, downstream_TSS)
		end_file_identifier = '{0}_{1}_{2}'.format(mode_atr2, upstream_TSS, downstream_TSS)

	if mode_atr2 == "GENE":	

		name_of_file_for_overlap = save_to_temp_folder + config_variables.name_of_time_series_promoter_file_for_TSS_start[7:-3] + "_GENE_{0}".format(upstream_TSS)
		create_GENE_file(upstream_TSS)
		end_file_identifier = '{0}_{1}'.format(mode_atr2, upstream_TSS)

	elif mode_atr2 == "ENHANCER":
		name_of_file_for_overlap = config_variables.name_of_enhancer_file_for_overlap
		end_file_identifier = '{0}'.format(mode_atr2)


	def create_enrichment_matrix():

		motif_enrichments = [[]]*len(np.loadtxt(name_of_file_for_overlap, dtype = str)) #tutaj trzeba zmienic na gena albo enhancera

		for name_of_file in name_of_files:

			name_of_file_ = pwd + "/" + hg + "/" + name_of_file	
			command_line = "windowBed -a {0} -b {1} -sw -l {2} -r {3}".format(name_of_file_, name_of_file_for_overlap, upstream_TSS, downstream_TSS) # name_of_enhancer_file_for_overlap zmienic na TSS w przypadku genow i dodac left right
			args = shlex.split(command_line)

			proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			output_raw = proc.stdout.read()
	
			if len(output_raw):
				output = np.array(map(lambda x: x.split("\t"), (output_raw).split("\n"))[:-1])

				np.savetxt(save_to_temp_folder + name_of_file[:-4] + '_overlap_{0}'.format(end_file_identifier), output, fmt = '%s', delimiter = '\t')
			
				for index_of_peak, peak_overlap in zip(output[:, -1].astype(int), output[:, :5]):
					motif_enrichments[index_of_peak] = motif_enrichments[index_of_peak] + [list(peak_overlap[[0,1,2,4]]) + [peak_overlap[3].split("_")[-1]] + [name_of_file.split("_")[1]] + [name_of_file.split("_")[0]]]


		file_1 = open(save_to_temp_folder + "enriched_peaks_{0}".format(end_file_identifier), 'w')
		peaks = np.loadtxt(name_of_file_for_overlap, dtype = str)

		for index in [ind for ind, el in enumerate(motif_enrichments) if len(el)]:
			array = motif_enrichments[index]
			for el in array:
		
				save = '\t'.join(np.r_[peaks[index], el])
				save += '\n'
				file_1.write(save)

		file_1.close()

		enriched_peaks = np.loadtxt(save_to_temp_folder + "enriched_peaks_{0}".format(end_file_identifier), str)
		legend = np.unique(enriched_peaks[:,-2])
		map_legend = {}
		for ind, el in enumerate(legend): map_legend[el] = ind

		count_matrix = np.zeros((len(motif_enrichments), len(legend)), bool)

		for el in enriched_peaks:
			count_matrix[int(el[-8]), map_legend[el[-2]]] = True

		np.save(save_to_temp_folder + "enrichment_matrix_{0}".format(end_file_identifier), count_matrix)
		return legend, count_matrix

	legend, count_matrix = create_enrichment_matrix()



	def sorts_labels():
		labels_count = np.histogram(labels, bins = range(0,max(labels)+2))[0][1:]
		sorted_counts_labels = np.argsort(labels_count)[::-1]
		sorted_counts = labels_count[sorted_counts_labels]

		sorted_labels = np.unique(labels)[sorted_counts_labels]


		def sorted_labels_func():
	
			time_series_survived = np.loadtxt(path_to_R + merged_time_series_to_cluster, dtype = np.float, delimiter = ",")

			means = []
			for ind, label in enumerate(sorted_labels):

				if mode_atr == "SECOND_TS": mean = (time_series_survived[label == labels, 8:]).mean(0)
				elif mode_atr == "FIRST_TS": mean = (time_series_survived[label == labels, :8]).mean(0)

				means += [mean]

			means = np.array(means)

			ind = np.lexsort((means[:,7],means[:,6],means[:,5],means[:,4],means[:,3],means[:,2],means[:,1],means[:,0]))

			amplitude = np.ravel(np.diff(means[:, [0,4]]))
			#amplitude = means[:, 4]/means[:, 0]
			if sorted_mode == "amplitude_sorted": ind = np.argsort(amplitude)[::-1]

			elif sorted_mode == "size_sorted": ind = np.arange(len(amplitude)).astype(int)

			return ind

		if diff_bind_version: ind_sort = [0,1]
		else: ind_sort = sorted_labels_func()

		sorted_labels = sorted_labels[ind_sort]
		sorted_counts = sorted_counts[ind_sort]

		return sorted_labels, sorted_counts, ind_sort

	sorted_labels, sorted_counts, ind_sort = sorts_labels()


	def calculates_probabilities_for_cluster():


		print count_matrix[survived].sum(0)/float(survived.shape[0])

		from scipy.stats import binom

		ps = count_matrix[survived].sum(0)/float(survived.shape[0])

		prob = np.zeros((len(np.unique(labels)), len(ps)))

		enrichments_counts = prob.astype(int)
		#sorts 	


		for index_1, label in enumerate(sorted_labels):

			n = np.sum(labels == label)
			xs = count_matrix[survived][labels == label].sum(0)
	
			for index_2, p in enumerate(ps):
				p = ps[index_2]
				x = xs[index_2]
				prob[index_1, index_2] = 1. - binom.cdf(x-1, n, p)

				enrichments_counts[index_1, index_2] = x

		np.savetxt(save_to_temp_folder + "_probabilities_of_enrichment_{0}".format(end_file_identifier), prob, delimiter = "\t", fmt = '%0.8f', header = '\t'.join(legend))
		return prob, enrichments_counts
	
	if GLOBAL_OR_SURVIVED == "survived": prob, enrichments_counts = calculates_probabilities_for_cluster()

	def calculates_probabilities_for_cluster_global():

		if mode_atr2 == "ENHANCER":

			distal_mask = np.invert(config_variables.proximal_enhancers_mask)

			print count_matrix[distal_mask].sum(0)/float(count_matrix[distal_mask].shape[0])
			
			ps = count_matrix[distal_mask].sum(0)/float(count_matrix[distal_mask].shape[0])

		elif mode_atr2 == "GENE" or mode_atr2 == "TSS":

			print count_matrix.sum(0)/float(count_matrix.shape[0])
			
			ps = count_matrix.sum(0)/float(count_matrix.shape[0])


		from scipy.stats import binom

		prob = np.zeros((len(np.unique(labels)), len(ps)))

		enrichments_counts = prob.astype(int)
		#sorts 	


		for index_1, label in enumerate(sorted_labels):

			n = np.sum(labels == label)
			xs = count_matrix[survived][labels == label].sum(0)
	
			for index_2, p in enumerate(ps):
				p = ps[index_2]
				x = xs[index_2]
				prob[index_1, index_2] = 1. - binom.cdf(x-1, n, p)

				enrichments_counts[index_1, index_2] = x

		np.savetxt(save_to_temp_folder + "_probabilities_of_enrichment_{0}".format(end_file_identifier), prob, delimiter = "\t", fmt = '%0.8f', header = '\t'.join(legend))
		return prob, enrichments_counts

	if GLOBAL_OR_SURVIVED == "global": prob, enrichments_counts = calculates_probabilities_for_cluster_global()
		


	mask_legend = np.ones_like(legend).astype(bool)
	mask_legend[np.in1d(legend, dont_plot)] = False
	

	file1 = open(save_to_temp_folder + merged_time_series_to_cluster + "_enrichment_{0}".format(end_file_identifier), "w")
	for i in range(len(prob)):
		 file1.write(','.join(legend[(prob[i] < 0.01)*mask_legend]) + "\n")

	file1.close()

	from matplotlib import pyplot as plt


	time_series_survived = np.loadtxt(path_to_R + merged_time_series_to_cluster, dtype = np.float, delimiter = ",")

	amplitude = np.zeros(len(sorted_labels))
	for ind, label in enumerate(sorted_labels):
		if mode_atr == "SECOND_TS": mean = (time_series_survived[label == labels, 8:]).mean(0) # tu trzeba to poprawic jesli chcesz dodac clustering dla geny
		elif mode_atr == "FIRST_TS": mean = (time_series_survived[label == labels, :8]).mean(0)

		#amplitude[ind] = mean[4]/mean[0]#np.diff(mean[[0,4]])
		amplitude[ind] = np.diff(mean[[0,4]])	

	idx = Index(np.unique(labels))
	df = DataFrame(np.c_[prob[:, mask_legend], sorted_counts[:,None], amplitude], index=idx, columns=np.r_[legend[mask_legend], ["Count"], ["Amplitude"]])
	vals = np.around(df.values,2)
	normal = plt.Normalize(prob[:,mask_legend].min(), prob[:,mask_legend].max())

	rise = np.zeros_like(vals).astype(bool)

	rise[:, :] = (amplitude > 0)[:,None]

	vals_enrich = np.c_[enrichments_counts[:,mask_legend], sorted_counts[:,None], (100*amplitude[:,None]).astype(int)]

	matrix_colour = plt.cm.hot(normal(vals))
	mask_encriched = np.c_[prob[:, mask_legend] < 0.01, np.ones((len(prob), 2), bool)]
	mask_encriched_2 = np.c_[prob[:, mask_legend] < 0.05, np.ones((len(prob), 2), bool)]
	mask_encriched_3 = np.c_[prob[:, mask_legend] < 0.001, np.ones((len(prob), 2), bool)]
	

	#matrix_colour[mask_encriched*rise] = np.array([0.0, 0.5019607843137255, 0.0, 0.6])
	#matrix_colour[mask_encriched*np.invert(rise)] = np.array([0.0, 0.5019607843137255, 0.0, 0.3])

	mask_depleted = np.c_[prob[:, mask_legend] > 0.99, np.ones((len(prob), 2), bool)]
	mask_depleted_2 = np.c_[prob[:, mask_legend] > 0.995, np.ones((len(prob), 2), bool)]
	mask_depleted_3 = np.c_[prob[:, mask_legend] > 0.999, np.ones((len(prob), 2), bool)]

	white = [1., 1., 1., 1.]

	mask_niether = np.invert(mask_encriched + mask_depleted)

	#matrix_colour[mask_depleted*rise] = np.array([0.768, 0.090, 0.090, 0.3])

	#matrix_colour[mask_depleted*np.invert(rise)] = np.array([0.768, 0.090, 0.090, 0.6])

	matrix_colour[mask_depleted] = [0.862745, 0.0784314, 0.235294, 0.7]
	matrix_colour[mask_depleted_2] = [0.862745, 0.0784314, 0.235294, 0.85]
	matrix_colour[mask_depleted_3] = [0.862745, 0.0784314, 0.235294, 1.]

	matrix_colour[mask_encriched] = [0.180392, 0.545098, 0.341176, .7]#[0., 1., 1., 1.]
	matrix_colour[mask_encriched_2] = [0.180392, 0.545098, 0.341176, .85]#[0., 1., 1., 1.]
	matrix_colour[mask_encriched_3] = [0.180392, 0.545098, 0.341176, 1.]#[0., 1., 1., 1.]

	matrix_colour[mask_niether] = [0.815, 0.803, 0.803, 1.]

	matrix_colour[:, -2] = white


	normal_2 = plt.Normalize(amplitude.min(), amplitude.max())

	amplitude_column = plt.cm.bwr_r(normal_2(amplitude))

	matrix_colour[:, -1] = amplitude_column

	#matrix_colour[rise[:,0], -1] = np.array([0.0, 0.5019607843137255, 0.0, 0.6])

	#matrix_colour[np.invert(rise[:,0]), -1] = np.array([0.768, 0.090, 0.090, 0.6])

	#fig = plt.figure(figsize=(12,10))
	#ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
	#the_table=plt.table(cellText=vals_enrich, rowLabels=df.index, colLabels=df.columns, 
	#                    colWidths = [0.07]*vals.shape[1], loc='center', 
	#                    cellColours=plt.get_cmap('Spectral')(normal(vals)))

	fig = plt.figure(figsize=(15,11))
	ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])

	#fig.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

	sm = plt.cm.ScalarMappable(cmap="bwr_r", norm=plt.Normalize(vmin=-1, vmax=1))
	sm._A = []
	#fig.colorbar(sm,shrink=0.25)#, ax = cbar_ax


	#cmap_r = mpl.cm.jet

	#ax1 = fig.add_axes([0.0, 0.9, 0.15])
	#norm = mpl.colors.Normalize(vmin=0, vmax=1)
	#cb1 = mpl.colorbar.ColorbarBase(ax1, cmap = cmap, norm=norm)



	the_table=plt.table(cellText=vals_enrich, rowLabels=ind_sort+1, colLabels=df.columns,
		                colWidths = [0.06]*vals.shape[1], rowLoc='right', loc='center left', 
		                cellColours=matrix_colour)

	#rowLabels=df.index
	import matplotlib.patches as mpatches


	line2a, = plt.plot([],[], label="enriched, p < 0.01", linewidth=15, color = [0.180392, 0.545098, 0.341176, 0.7])#[0., 1., 1., 1.]
	line2b, = plt.plot([],[], label="enriched, p < 0.005", linewidth=15, color = [0.180392, 0.545098, 0.341176, 0.85])
	line2c, = plt.plot([],[], label="enriched, p < 0.001", linewidth=15, color = [0.180392, 0.545098, 0.341176, 1.])

	line3a, = plt.plot([],[], label="depleted, p < 0.01", linewidth=15, color = [0.862745, 0.0784314, 0.235294, 0.7])
	line3b, = plt.plot([],[], label="depleted, p < 0.005", linewidth=15, color = [0.862745, 0.0784314, 0.235294, 0.85])
	line3c, = plt.plot([],[], label="depleted, p < 0.001", linewidth=15, color = [0.862745, 0.0784314, 0.235294, 1.])

	line1, = plt.plot([],[], label="neither", linewidth=15, color = [0.815, 0.803, 0.803, 1.])

	#line4, = plt.plot([],[], label="rises between 0-40min", linewidth=15, color = [0.0, 0.5019607843137255, 0.0, 0.6])
	#line5, = plt.plot([],[], label="drops between 0-40min", linewidth=15, color = [0.768, 0.090, 0.090, 0.6])

	#line2, = plt.plot([],[], label="enriched & rise", linewidth=15, color = [0.0, 0.5019607843137255, 0.0, 0.6])
	#line3, = plt.plot([],[], label="enriched & drops", linewidth=15, color = [0.0, 0.5019607843137255, 0.0, 0.3])

	#line4, = plt.plot([],[], label="depleted & rise", linewidth=15, color = [0.768, 0.090, 0.090, 0.3])
	#line5, = plt.plot([],[], label="depleted & drops", linewidth=15, color = [0.768, 0.090, 0.090, 0.6])

	#line1, = plt.plot([],[], label="neither", linewidth=15, color = [0.815, 0.803, 0.803, 1.])


	#fig.patch.set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)

	table_props = the_table.properties()
	table_cells = table_props['child_artists']
	for cell in table_cells: cell.set_height(1.2*cell.get_height())
	#plt.legend()

	first_legend = plt.legend(bbox_to_anchor=(0.96, 1))
	ax = plt.gca().add_artist(first_legend)

	cbaxes = fig.add_axes([0.685, 0.45, 0.025, 0.2])
	cb = fig.colorbar(sm, cax = cbaxes, ticks=[-1, 0, 1])#, shrink=2.)  
	cb.ax.set_yticklabels(['drops between 0-40min', 'stationary', 'rises between 0-40min'])


	#plt.text(2, 6, r'an equation: $E=mc^2$', fontsize=15)

	#plt.title("Transcription Factors", fontsize=20)
	#plt.ylabel('Clusters', fontsize=20)
	#plt.xlabel("distance [B]", fontsize=20)
	if diff_bind_version: name_save = '{0}TF_enrichment_{1}.pdf'.format(save_to_temp_folder, "diff_bind")
	else: name_save ='{0}TF_enrichment_{1}_{2}_{3}_{4}_0_40.pdf'.format(save_to_temp_folder, end_file_identifier, mode_atr, sorted_mode, mode_of_data_sets)

	pdf = PdfPages(name_save); pdf.savefig()
	pdf.close(); plt.close('all')

