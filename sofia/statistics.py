import os
import numpy as np
import ctypes as ct



# =========================
# Load C statistics library
# =========================

_stat = ct.CDLL(os.environ["SOFIA_MODULE_PATH"] + "/sofia/_statistics.so")



# ============================================
# Declare types of arguments and return values
# ============================================

# Summation
# ---------
_stat.sum.argtypes = [ct.POINTER(ct.c_double), ct.c_size_t]
_stat.sum.restype = ct.c_double

# Moment map
# ----------
_stat.moment.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_size_t, ct.c_size_t, ct.c_uint, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)]
_stat.moment.restype = ct.POINTER(ct.c_double)

# Memory de-allocation
# --------------------
_stat.free_memory.argtypes = [ct.POINTER(ct.c_double)]
_stat.free_memory.restype = None



# ========================
# Python function wrappers
# ========================

# Summation
# ---------
def nansum(data):
	global _stat
	
	# Prepare arguments
	ptr_data = data.ctypes.data_as(ct.POINTER(ct.c_double))
	n = ct.c_size_t(data.size)
	
	# Call sum function
	return _stat.sum(ptr_data, n)


# Moment maps
# -----------
def moment(data, mom=0, mom0=None, mom1=None):
	global _stat
	
	# Prepare arguments
	ptr_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	nx = ct.c_size_t(data.shape[0])
	ny = ct.c_size_t(data.shape[1])
	nz = ct.c_size_t(data.shape[2])
	mm = ct.c_uint(mom)
	if mom0 is not None: ptr_mom0 = mom0.ctypes.data_as(ct.POINTER(ct.c_double))
	else: ptr_mom0 = None
	if mom1 is not None: ptr_mom1 = mom1.ctypes.data_as(ct.POINTER(ct.c_double))
	else: ptr_mom1 = None
	
	# Call moment function (returns pointer)
	ptr_moment_map = _stat.moment(ptr_data, nx, ny, nz, mm, ptr_mom0, ptr_mom1)
	
	# Copy data into new NumPy array
	moment_map = np.array(np.fromiter(ptr_moment_map, dtype=np.float64, count=data.shape[0]*data.shape[1])).reshape((data.shape[0], data.shape[1]))
	
	# De-allocate memory
	_stat.free_memory(ptr_moment_map)
	
	return moment_map
