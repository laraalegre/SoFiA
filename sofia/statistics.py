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

# Check for NaN
# -------------
_stat.check_nan.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_uint]
_stat.check_nan.restype = ct.c_uint

# Replace NaN
# -----------
_stat.replace_nan.argtypes = [ct.POINTER(ct.c_float), ct.POINTER(ct.c_float), ct.c_size_t, ct.c_float, ct.c_uint, ct.c_uint]
_stat.replace_nan.restype = None

# Set mask based on threshold
# ---------------------------
_stat.set_mask.argtypes = [ct.POINTER(ct.c_ubyte), ct.POINTER(ct.c_float), ct.c_size_t, ct.c_float, ct.c_uint]
_stat.set_mask.restype = None

# Standard deviation
# ------------------
_stat.standard_deviation.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_size_t, ct.c_uint, ct.c_uint]
_stat.standard_deviation.restype = ct.c_double

# Summation
# ---------
_stat.sum.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t]
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

# Check for NaN
# -------------
def check_nan(data):
	global _stat
	
	# Prepare arguments
	ptr_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	n        = ct.c_size_t(data.size)
	bo       = ct.c_uint(get_byte_order(data))
	
	# Call sum function
	return _stat.check_nan(ptr_data, n, bo)


# Replace NaN
# -----------
def replace_nan(data, mask, value):
	global _stat
	
	# Prepare arguments
	ptr_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	ptr_mask = mask.ctypes.data_as(ct.POINTER(ct.c_float))
	n        = ct.c_size_t(data.size)
	val      = ct.c_float(value)
	bo_data  = ct.c_uint(get_byte_order(data))
	bo_mask  = ct.c_uint(get_byte_order(mask))
	
	# Call sum function
	return _stat.replace_nan(ptr_data, ptr_mask, n, val, bo_data, bo_mask)


# Set mask based on threshold
# ---------------------------
def set_mask(mask, data, threshold):
	global _stat
	
	# Prepare arguments
	ptr_mask = mask.ctypes.data_as(ct.POINTER(ct.c_ubyte))
	ptr_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	n        = ct.c_size_t(data.size)
	thresh   = ct.c_float(threshold)
	bo       = ct.c_uint(get_byte_order(data))
	
	# Call sum function
	return _stat.set_mask(ptr_mask, ptr_data, n, thresh, bo)


# Standard deviation
# ------------------
def standard_deviation(data, flux_range, cadence):
	global _stat
	
	# Define flux range and method values
	flux_ranges = {
		"all": 0,
		"negative": 1,
		"positive": 2}
	
	# Prepare arguments
	arg_data    = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size    = ct.c_size_t(data.size)
	arg_cadence = ct.c_size_t(cadence)
	arg_range   = ct.c_uint(flux_ranges[flux_range])
	arg_bo      = ct.c_uint(get_byte_order(data))
	
	# Call C function	
	return _stat.standard_deviation(arg_data, arg_size, arg_cadence, arg_range, arg_bo)


# Summation
# ---------
def nansum(data):
	global _stat
	
	# Prepare arguments
	ptr_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
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


# Determine byte order of data
# ----------------------------

def get_byte_order(data):
	byte_order = 0
	if(data.dtype.byteorder == ">"):
		byte_order = 1
	return byte_order
