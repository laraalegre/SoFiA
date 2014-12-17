# distutils: language = c++

# cython: boundscheck = False
# cython: wraparound = False
# cython: cdivision = True

cimport numpy as np
import numpy as np

cdef extern from "CNHI.h":
     
     cdef void CNHI_find_sources(float * data, int * mask, int verbose, int metric[9], float p_req, int tr_min, int tr_max, int flag_val, bint median_test, float q_req)

def find_sources(data, mask, p_req = 0.05, min_scale = 4, max_scale = -1, verbose = 1, median_test = 1, q_req = 3.8, flag_val = 1):
     
     return _find_sources(data.astype(np.single, copy = False), mask.astype(np.intc, copy = False), p_req, min_scale, max_scale, verbose, flag_val, median_test, q_req)

cdef _find_sources(np.ndarray[dtype = float, ndim = 3] data, np.ndarray[dtype = int, ndim = 3] mask, float p_req, int min_scale, int max_scale, int verbose, int flag_val, int median_test, float q_req):
     
     # define variables
     cdef bint c_median_test
     cdef int metric[9]

     # set the boolean flag, c_median_test, that sets whether detections must have a median greater than the median of the remaining data
     c_median_test = 0
     if (median_test > 0):
          c_median_test = 1
    
     # set the metric used to parse the datacube
     metric[0] = data.shape[2]
     metric[1] = data.shape[1]
     metric[2] = data.shape[0]
     metric[3] = 0
     metric[4] = 0
     metric[5] = 0
     metric[6] = data.shape[2]
     metric[7] = data.shape[1]
     metric[8] = data.shape[0]

     # call the function
     CNHI_find_sources(<float *> data.data, <int *> mask.data, verbose, metric, <float> p_req, min_scale, max_scale, <int> flag_val, c_median_test, <float> q_req)

     return mask

