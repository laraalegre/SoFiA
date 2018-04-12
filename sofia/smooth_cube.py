#! /usr/bin/env python

import numpy as np
from scipy import ndimage
from sofia import error as err


def smooth(indata, kernel, edgeMode, kernelX, kernelY, kernelZ):
	"""
	Smooth a data cube with the specified kernel type and size.
	
	Arguments:
	  indata:       The input data cube.
	  kernel:       The smoothing kernel; "gaussian", "boxcar" or "median".
	  edgeMode:     Determines how borders are handled; "reflect", "constant", "nearest", "mirror" or "wrap".
	  kernelX/Y/Z:  Size of the kernel (standard deviation in the case of a Gaussian kernel).
	
	Returns:
	  Smoothed copy of the data cube.
	"""
	
	err.message("Smoothing data cube")
	
	# Sanity checks of user input
	err.ensure(
		kernel in {"gaussian", "boxcar", "median"},
		"Smoothing failed. Illegal smoothing type: '" + str(kernel) + "'.")
	err.ensure(
		edgeMode in {"reflect", "constant", "nearest", "mirror", "wrap"},
		"Smoothing failed. Illegal edge mode: '" + str(edgeMode) + "'.")
	err.ensure(
		kernelX or kernelY or kernelZ,
		"Smoothing failed. All smoothing kernels are zero.")
	err.ensure(
		kernel != "median" or (kernelX and kernelY and kernelZ),
		"Smoothing failed. Cannot determine median for kernel size of zero.")
	
	# Print some information
	err.message("  Kernel type: " + str(kernel).title())
	err.message("  Kernel size: [" + str(kernelX) + ", " +  str(kernelY) + ", " + str(kernelZ) + "]")
	err.message("  Edge mode:   " + str(edgeMode))
	
	# Create copy of input cube to be smoothed
	outdata = np.copy(indata)
	
	# Remove NaNs (and INFs) if necessary
	found_nan = np.isnan(indata).sum()
	if found_nan: outdata = np.nan_to_num(outdata, copy=False)
	
	# Smooth with the selected kernel
	if kernel == "gaussian":
		outdata = ndimage.filters.gaussian_filter(outdata, sigma=(kernelZ, kernelX, kernelY), mode=edgeMode)
	elif kernel == "boxcar":
		outdata = ndimage.filters.uniform_filter(outdata, size=(kernelZ, kernelX, kernelY), mode=edgeMode)
	else: # kernel == "median"
		outdata = ndimage.filters.median_filter(outdata, size=(kernelZ, kernelX, kernelY), mode=edgeMode)
	
	# Put NaNs back in if necessary
	if found_nan: outdata[np.isnan(indata)] = np.nan
	
	return outdata
