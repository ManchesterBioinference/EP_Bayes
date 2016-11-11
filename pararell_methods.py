import numpy as np
import numexpr as ne
def pararell_calc(smaller_c_trace):
	dim = np.shape(smaller_c_trace)
	
	return ( smaller_c_trace.reshape(dim[0], dim[1], 1) == smaller_c_trace.reshape(dim[0], 1, dim[1]) ).sum(0)

def pararell_calc_ne(smaller_c_trace):
	dim = np.shape(smaller_c_trace)

	A = smaller_c_trace.reshape(dim[0], dim[1], 1)
	B = smaller_c_trace.reshape(dim[0], 1, dim[1])
	
	return ( ne.evaluate("A == B") ).sum(0)


		#t_new = time.time()
		#total_matrix = np.zeros((_c_trace.shape[1], _c_trace.shape[1]), dtype = int)

		#for im, trace in enumerate(_c_trace):
		#	print im
		#	match_matrix = trace.reshape(len(trace),1) == trace.reshape(1,len(trace)) # meaning: copy vector vertically n times,  copy vector horizontally n times. 
		#	total_matrix = total_matrix + match_matrix.astype(int)

		#print time.time() - t_new


		#t_new = time.time()

		#total_matrix_2 = np.zeros((_c_trace.shape[1], _c_trace.shape[1]), dtype = int)
		#dim = _c_trace.shape
		#pack = 40
		#incr = int(dim[0]/pack)

		#for i in range(pack): total_matrix_2 = total_matrix_2 + ( _c_trace.reshape(dim[0], dim[1], 1)[i*incr:(i+1)*incr] == _c_trace.reshape(dim[0], 1, dim[1])[i*incr: (i+1)*incr] ).sum(0) # meaning: copy vector vertically n times,  copy vector horizontally n times. 
		#missing = dim[0] - (i+1)*incr
		#if missing: total_matrix_2 = total_matrix_2 + ( _c_trace.reshape(dim[0], dim[1], 1)[(i+1)*incr: dim[0]] == _c_trace.reshape(dim[0], 1, dim[1])[(i+1)*incr: dim[0]] ).sum(0)

		#print time.time() - t_new

		#print "all:", (total_matrix == total_matrix_2).all()
		


