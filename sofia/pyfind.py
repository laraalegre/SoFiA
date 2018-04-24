#! /usr/bin/env python

import math
from time import time
import numpy as np
from scipy import ndimage
from sofia.functions import GetRMS
from sofia import error as err


def SCfinder_mem(cube, header, t0, kernels=[[0, 0, 0, "b"],], threshold=3.5, sizeFilter=0, maskScaleXY=2.0, maskScaleZ=2.0, kernelUnit="pixel", edgeMode="constant", rmsMode="negative", fluxRange="all", verbose=0):
	# Define a few constants
	FWHM_CONST    = 2.0 * math.sqrt(2.0 * math.log(2.0))   # Conversion between sigma and FWHM of Gaussian function
	MAX_PIX_CONST = 1e+6                                   # Maximum number of pixels for noise calculation; sampling is set accordingly
	
	# Create binary mask array
	msk = np.zeros(cube.shape, np.bool)
	found_nan = np.isnan(cube).sum()
	
	# Set sampling sampleRms for rms measurement
	sampleRms = max(1, int((float(np.array(cube.shape).prod()) / MAX_PIX_CONST)**(1.0 / min(3, len(cube.shape)))))
	
	# Measure noise in original cube with sampling "sampleRms"
	rms = GetRMS(cube, rmsMode=rmsMode, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=verbose, sample=sampleRms)
	
	# Loop over all kernels
	for kernel in kernels:
		[kx, ky, kz, kt] = kernel
		if verbose:
			err.linebreak()
			err.print_progress_time(t0)
			err.message("    Filter %s %s %s %s ..." % (kx, ky, kz, kt))
		if kernelUnit == "world" or kernelUnit == "w":
			if verbose: err.message("    Converting filter size to pixels ...")
			kx = abs(float(kx) / header["CDELT1"])
			ky = abs(float(ky) / header["CDELT2"])
			kz = abs(float(kz) / header["CDELT3"])
		if kt == "b":
			if kz != int(math.ceil(kz)) and verbose: err.warning("Rounding width of boxcar z kernel to next integer.")
			kz = int(math.ceil(kz))
		
		# Create a copy of the original cube:
		smoothedCube = np.copy(cube)
		
		# Replace all NaNs with zero (and INFs with a finite number):
		if found_nan: smoothedCube = np.nan_to_num(smoothedCube)
		
		smoothedCube[(smoothedCube > 0) & (msk > 0)] = +maskScaleXY * rms
		smoothedCube[(smoothedCube < 0) & (msk > 0)] = -maskScaleXY * rms
		
		# Spatial smoothing:
		if kx + ky:
			smoothedCube = ndimage.filters.gaussian_filter(smoothedCube, [0, ky / FWHM_CONST, kx / FWHM_CONST], mode=edgeMode)
		
		# Spectral smoothing:
		if kz:
			if   kt == "b": smoothedCube = ndimage.filters.uniform_filter1d(smoothedCube, kz, axis=0, mode=edgeMode)
			elif kt == "g": smoothedCube = ndimage.filters.gaussian_filter1d(smoothedCube, kz / FWHM_CONST, axis=0, mode=edgeMode)
		
		# Re-insert the NaNs (but not the INFs) taken out earlier:
		if found_nan: smoothedCube[np.isnan(cube)] = np.nan
		
		# Calculate the RMS of the smoothed cube:
		smoothedrms = GetRMS(smoothedCube, rmsMode=rmsMode, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=verbose, sample=sampleRms)
		
		# Get rid of the NaNs a second time:
		#if found_nan: smoothedCube = np.nan_to_num(smoothedCube)
		# NOTE: This should not be necessary because any comparison with NaN will always yield False.
		#       Hence, NaN pixels will never be included in the mask below.
		
		# Add pixels above threshold to mask by setting bit 1:
		with np.errstate(invalid="ignore"):
			msk = np.bitwise_or(msk, np.greater_equal(np.absolute(smoothedCube), threshold * smoothedrms))
		
		# Delete smoothed cube again:
		del(smoothedCube)
	return msk
