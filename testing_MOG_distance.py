test_2 = probabilities_for_promoters_of_interacting_enhancers["MOG_dist"]["chr2"] distance prior only

match_MAP["MOG_dist"]["chr2"]

prob and match_MAP["dist"]["chr2"] is equivalent to NB




test_2.argmax(0)[np.invert(match_MAP["MOG_dist"]["chr2"])] prior_only

prob.argmax(0)[np.invert(match_MAP["dist"]["chr2"])]




gives_sum_of_uncommon_elements = NB_match - NB_match * prior_match

own_by_NB = NB_match - NB_match * prior_match

own_by_prior = prior_match - NB_match * prior_match
