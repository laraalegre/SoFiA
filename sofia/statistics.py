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
_stat.check_nan.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t]
_stat.check_nan.restype = ct.c_uint

# Set mask based on threshold
# ---------------------------
_stat.set_mask.argtypes = [ct.POINTER(ct.c_ubyte), ct.POINTER(ct.c_float), ct.c_size_t, ct.c_float]
_stat.set_mask.restype = None

# Standard deviation
# ------------------
_stat.stddev.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_size_t, ct.c_int, ct.c_float]
_stat.stddev.restype = ct.c_double

# Median
# ------
_stat.median.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_uint]
_stat.median.restype = ct.c_float

# Median absolute deviation
# -------------------------
_stat.mad.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_float]
_stat.mad.restype = ct.c_float

# Summation
# ---------
_stat.sum.argtypes = [ct.POINTER(ct.c_float), ct.c_size_t, ct.c_uint]
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
	arg_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size = ct.c_size_t(data.size)
	
	# Call C function
	return _stat.check_nan(arg_data, arg_size)


# Set mask based on threshold
# ---------------------------
def set_mask(mask, data, threshold):
	global _stat
	
	# Prepare arguments
	arg_mask      = mask.ctypes.data_as(ct.POINTER(ct.c_ubyte))
	arg_data      = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size      = ct.c_size_t(data.size)
	arg_threshold = ct.c_float(threshold)
	
	# Call C function
	return _stat.set_mask(arg_mask, arg_data, arg_size, arg_threshold)


# Standard deviation
# ------------------
def stddev(data, flux_range="all", cadence=1, value=np.nan):
	global _stat
	
	# Define flux range and method values
	flux_ranges = {
		"all": 0,
		"negative": -1,
		"positive": 1}
	
	# Prepare arguments
	arg_data    = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size    = ct.c_size_t(data.size)
	arg_cadence = ct.c_size_t(cadence)
	arg_range   = ct.c_int(flux_ranges[flux_range])
	arg_value   = ct.c_float(value)
	
	# Call C function
	return _stat.stddev(arg_data, arg_size, arg_cadence, arg_range, arg_value)


# Median
# ------
def median(data, flux_range="all", cadence=1):
	global _stat
	
	#if flux_range == "negative":
	#	data = np.array(data[data < 0][0::cadence], dtype=data.dtype)
	#elif flux_range == "positive":
	#	data = np.array(data[data > 0][0::cadence], dtype=data.dtype)
	#else:
	#	data = np.array(data[~np.isnan(data)][0::cadence], dtype=data.dtype)
	
	# Prepare arguments
	arg_data    = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size    = ct.c_size_t(data.size)
	arg_approx  = ct.c_uint(0)
	
	# Call C function
	return _stat.median(arg_data, arg_size, arg_approx)


# Median absolute deviation
# -------------------------
def mad(data, flux_range="all", cadence=1, value=np.nan):
	global _stat
	
	if flux_range == "negative":
		data = np.array(data[data < 0][0::cadence], dtype=data.dtype)
	elif flux_range == "positive":
		data = np.array(data[data > 0][0::cadence], dtype=data.dtype)
	else:
		data = np.array(data[~np.isnan(data)][0::cadence], dtype=data.dtype)
	
	# Prepare arguments
	arg_data  = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size  = ct.c_size_t(data.size)
	arg_value = ct.c_float(value)
	
	# Call C function
	return _stat.mad(arg_data, arg_size, arg_value)


# Summation
# ---------
def sum(data):
	global _stat
	
	# Prepare arguments
	arg_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size = ct.c_size_t(data.size)
	arg_mean = ct.c_uint(0)
	
	# Call C function
	return _stat.sum(arg_data, arg_size, arg_mean)


# Mean
# ----
def mean(data):
	global _stat
	
	# Prepare arguments
	arg_data = data.ctypes.data_as(ct.POINTER(ct.c_float))
	arg_size = ct.c_size_t(data.size)
	arg_mean = ct.c_uint(1)
	
	# Call C function
	return _stat.sum(arg_data, arg_size, arg_mean)


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
	
	# Call C function (returns pointer)
	ptr_moment_map = _stat.moment(ptr_data, nx, ny, nz, mm, ptr_mom0, ptr_mom1)
	
	# Copy data into new NumPy array
	moment_map = np.array(np.fromiter(ptr_moment_map, dtype=np.float64, count=data.shape[0]*data.shape[1])).reshape((data.shape[0], data.shape[1]))
	
	# Release memory
	_stat.free_memory(ptr_moment_map)
	
	return moment_map


# Determine byte order of data
# ----------------------------

def get_byte_order(data):
	return (data.dtype.byteorder == ">")
