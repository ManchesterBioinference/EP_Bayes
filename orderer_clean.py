def executor(dataset_name, chrom_names):
	import numpy as np
	import itertools

	dataset = np.loadtxt(dataset_name, delimiter = '\t', dtype = str)
	mask_ordered = list(itertools.chain.from_iterable([np.where(dataset[:,0] == chrom)[0] for chrom in chrom_names]))
	dataset_ordered = dataset[mask_ordered]
	np.savetxt(dataset_name + '_ordered', dataset_ordered, fmt='%s', delimiter='\t')
