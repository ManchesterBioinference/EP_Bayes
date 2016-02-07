import os
import shlex, subprocess
import numpy as np
import re
import config_variables
from matplotlib.backends.backend_pdf import PdfPages
from pandas import *


def executor(merged_time_series_to_cluster, diff_bind_version = False, mode_atr = ["ER", "Pol2"][0]):

	pwd = os.getcwd()
	hg = 'hg19'

	name_of_enhancer_file_for_overlap = config_variables.name_of_enhancer_file_for_overlap

	#config_variables.merged_time_series_to_cluster = merged_time_series_to_cluster

	#merged_time_series_to_cluster = "common_region_peaks_extended_less_time_points_corrected_0_indexed_unfiltered_count_concat_PolII_ER_200"

	name_of_files = np.loadtxt(pwd + "/" + hg + "/list_of_files.txt", dtype = str)

	survived = np.loadtxt('{0}_survived_indexes'.format(merged_time_series_to_cluster)).astype(int) # saved during filtering


	if diff_bind_version:
		peaks = np.loadtxt(name_of_enhancer_file_for_overlap, dtype = str)
		indexes_of_DB_peaks = np.loadtxt(pwd + "/" + hg + "/indexes_of_DB_peaks.csv", dtype = int, skiprows = 1, usecols = (1,), delimiter = ",")
		labels = np.zeros(len(peaks), int)
		labels[indexes_of_DB_peaks] = 1
		labels = labels + 1
		labels = labels[survived]

	else:
		labels = np.loadtxt('{0}_labels'.format(merged_time_series_to_cluster), str, delimiter = ",")[1:,1].astype(int) # from EP clustering

	save_to_temp_folder = pwd + "/" + hg +"/" + merged_time_series_to_cluster + "_results_temp/"
	if not os.path.exists(save_to_temp_folder):
		os.makedirs(save_to_temp_folder)

	def create_enrichment_matrix():
		motif_enrichments = [[]]*len(np.loadtxt(name_of_enhancer_file_for_overlap, dtype = str))

		for name_of_file in name_of_files:

			name_of_file_ = pwd + "/" + hg + "/" + name_of_file	
			command_line = "windowBed -a {0} -b {1} -sw -l 0 -r 0".format(name_of_file_, name_of_enhancer_file_for_overlap)
			args = shlex.split(command_line)

			proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			output_raw = proc.stdout.read()
	
			if len(output_raw):
				output = re.findall(r"[\w^+^-^//^,^.]+", output_raw)
				output = np.array(output).reshape(len(output)/9, 9)
				np.savetxt(save_to_temp_folder + name_of_file[:-4] + '_ER_peaks_overlap.txt', output, fmt = '%s', delimiter = '\t')
			
				for index_of_peak, peak_overlap in zip(output[:, 8].astype(int), output[:, :5]):
					motif_enrichments[index_of_peak] = motif_enrichments[index_of_peak] + [list(peak_overlap[[0,1,2,4]]) + [peak_overlap[3].split("_")[-1]] + [name_of_file.split("_")[1]] + [name_of_file.split("_")[0]]]


		file_1 = open(save_to_temp_folder + "enriched_peaks", 'w')
		peaks = np.loadtxt(name_of_enhancer_file_for_overlap, dtype = str)
		for index in [ind for ind, el in enumerate(motif_enrichments) if len(el)]:
			array = motif_enrichments[index]
			for el in array:
		
				save = '\t'.join(np.r_[peaks[index], el])
				save += '\n'
				file_1.write(save)

		file_1.close()

		enriched_peaks = np.loadtxt(save_to_temp_folder + "enriched_peaks", str)
		legend = np.unique(enriched_peaks[:,9])
		map_legend = {}
		for ind, el in enumerate(legend): map_legend[el] = ind

		count_matrix = np.zeros((len(motif_enrichments), len(legend)), bool)

		for el in enriched_peaks:
			count_matrix[int(el[3]), map_legend[el[-2]]] = True

		np.save("enrichment_matrix", count_matrix)
		return legend, count_matrix

	legend, count_matrix = create_enrichment_matrix()



	def sorts_labels():
		labels_count = np.histogram(labels, bins = range(0,max(labels)+2))[0][1:]
		sorted_counts_labels = np.argsort(labels_count)[::-1]
		sorted_counts = labels_count[sorted_counts_labels]

		sorted_labels = np.unique(labels)[sorted_counts_labels]


		def sorted_labels_func():
	
			time_series_survived = np.loadtxt(merged_time_series_to_cluster, dtype = np.float, delimiter = ",")

			means = []
			for ind, label in enumerate(sorted_labels):

				if mode_atr == "ER": mean = (time_series_survived[label == labels, :8]).mean(0)
				elif mode_atr == "Pol2": mean = (time_series_survived[label == labels, 8:]).mean(0)

				means += [mean]

			means = np.array(means)

			ind = np.lexsort((means[:,7],means[:,6],means[:,5],means[:,4],means[:,3],means[:,2],means[:,1],means[:,0]))
			return ind

		if diff_bind_version: ind_sort = [0,1]
		else: ind_sort = sorted_labels_func()

		sorted_labels = sorted_labels[ind_sort]
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

		np.savetxt("{0}_probabilities_of_enrichment".format(merged_time_series_to_cluster), prob, delimiter = "\t", fmt = '%0.8f', header = '\t'.join(legend))
		return prob, enrichments_counts

	prob, enrichments_counts = calculates_probabilities_for_cluster()
		


	mask_legend = np.ones_like(legend).astype(bool)
	mask_legend[2:4] = False

	#mask_legend[10:18] = False
	#mask_legend[25:27] = False

	file1 = open(save_to_temp_folder + merged_time_series_to_cluster + "_enrichment", "w")
	for i in range(len(prob)):
		 file1.write(','.join(legend[(prob[i] < 0.01)*mask_legend]) + "\n")

	file1.close()

	from matplotlib import pyplot as plt
	

	idx = Index(np.unique(labels))
	df = DataFrame(np.c_[prob[:, mask_legend], sorted_counts[ind_sort][:,None]], index=idx, columns=np.r_[legend[mask_legend], ["Count"]])
	vals = np.around(df.values,2)
	normal = plt.Normalize(prob[:,mask_legend].min()-0.3, prob[:,mask_legend].max()+0.3)

	rise = np.zeros_like(vals).astype(bool)
	time_series = np.loadtxt(merged_time_series_to_cluster, delimiter = ",")

	for ind, label in enumerate(sorted_labels):
		if mode_atr == "ER": mean = (time_series[label == labels, :8]).mean(0)
		elif mode_atr == "Pol2": mean = (time_series[label == labels, 8:]).mean(0)

		#mean = (time_series[label == labels]).mean(0)

		if mean[0] < mean[1]: rise[ind, :-1] = np.ones(len(legend[mask_legend]))

	vals_enrich = np.c_[enrichments_counts[:,mask_legend], sorted_counts[ind_sort][:,None]]

	matrix_colour = plt.cm.hot(normal(vals))
	mask_encriched = vals < 0.01
	matrix_colour[mask_encriched*rise] = np.array([0.0, 0.5019607843137255, 0.0, 0.6])

	matrix_colour[mask_encriched*np.invert(rise)] = np.array([0.0, 0.5019607843137255, 0.0, 0.3])

	mask_depleted = vals > 0.99

	white = [1., 1., 1., 1.]

	mask_niether = np.invert(mask_encriched + mask_depleted)

	matrix_colour[mask_depleted*rise] = np.array([0.768, 0.090, 0.090, 0.3])

	matrix_colour[mask_depleted*np.invert(rise)] = np.array([0.768, 0.090, 0.090, 0.6])

	matrix_colour[mask_niether] = [0.815, 0.803, 0.803, 1.]

	matrix_colour[:, -1] = white


	#fig = plt.figure(figsize=(12,10))
	#ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
	#the_table=plt.table(cellText=vals_enrich, rowLabels=df.index, colLabels=df.columns, 
	#                    colWidths = [0.07]*vals.shape[1], loc='center', 
	#                    cellColours=plt.get_cmap('Spectral')(normal(vals)))

	fig = plt.figure(figsize=(15,11))
	ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
	the_table=plt.table(cellText=vals_enrich, rowLabels=df.index, colLabels=df.columns, 
		                colWidths = [0.06]*vals.shape[1], rowLoc='right', loc='center left', 
		                cellColours=matrix_colour)


	import matplotlib.patches as mpatches


	line2, = plt.plot([],[], label="enriched & rise", linewidth=15, color = [0.0, 0.5019607843137255, 0.0, 0.6])
	line3, = plt.plot([],[], label="enriched & drops", linewidth=15, color = [0.0, 0.5019607843137255, 0.0, 0.3])

	line4, = plt.plot([],[], label="depleted & rise", linewidth=15, color = [0.768, 0.090, 0.090, 0.3])
	line5, = plt.plot([],[], label="depleted & drops", linewidth=15, color = [0.768, 0.090, 0.090, 0.6])

	line1, = plt.plot([],[], label="neither", linewidth=15, color = [0.815, 0.803, 0.803, 1.])


	#fig.patch.set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)

	table_props = the_table.properties()
	table_cells = table_props['child_artists']
	for cell in table_cells: cell.set_height(1.2*cell.get_height())
	#plt.legend()
	first_legend = plt.legend(bbox_to_anchor=(0.9, 1))
	ax = plt.gca().add_artist(first_legend)

	#plt.text(2, 6, r'an equation: $E=mc^2$', fontsize=15)

	#plt.title("Transcription Factors", fontsize=20)
	#plt.ylabel('Clusters', fontsize=20)
	#plt.xlabel("distance [B]", fontsize=20)
	if diff_bind_version: name_save = '{0}TF_enrichment_{1}.pdf'.format(save_to_temp_folder, "diff_bind")
	else: name_save ='{0}TF_enrichment_{1}.pdf'.format(save_to_temp_folder, mode_atr)

	pdf = PdfPages(name_save); pdf.savefig()
	pdf.close(); plt.show()



#count_matrix[survived][labels == 1].sum(0)/float(count_matrix[survived].sum(0)[0])


#for i in $( ls ); do
#	echo item: $i
#done
