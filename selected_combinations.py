def selected_combinations(PR_CURVES):
	import itertools
	import numpy as np
	stuff = [1, 2, 3, 4]
	combinations = []
	for L in range(0, len(stuff)+1):
		for subset in itertools.combinations(stuff, L):
			if len(subset): combinations += [list(subset)]

	if PR_CURVES == "ALL": 
		selected_combinations_ = combinations
	elif PR_CURVES == "SELECTIVE":
		selected_combinations_ = np.array(combinations)[[0, 2, 5, 10, 14]].tolist()

	else: selected_combinations_ = "something_wrong"

	return combinations, selected_combinations_
