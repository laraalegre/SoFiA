#! /usr/bin/env python

import math
import numpy as np
from scipy import ndimage
from sofia.functions import GetRMS
from sofia import error as err
from sofia.sigma_cube import sigma_scale



# ==========================================
# FUNCTION: Implementation of the S+C finder
# ==========================================

def SCfinder_mem(cube, mask, header, t0, kernels=[[0, 0, 0, "b"],], threshold=3.5, sizeFilter=0, maskScaleXY=2.0, maskScaleZ=2.0, kernelUnit="pixel", edgeMode="constant", rmsMode="negative", fluxRange="all", verbose=0, scaleEachSCkernel=False, scaleX=False, scaleY=False, scaleZ=True, edgeX=0, edgeY=0, edgeZ=0, method="1d2d", windowSpatial=20, windowSpectral=20, gridSpatial=0, gridSpectral=0, interpolation="none"):
	# Define a few constants
	FWHM_CONST    = 2.0 * math.sqrt(2.0 * math.log(2.0))   # Conversion between sigma and FWHM of Gaussian function
	MAX_PIX_CONST = 1.0e+6                                 # Maximum number of pixels for noise calculation; sampling is set accordingly
	
	# Check for NaN in cube
	found_nan = np.isnan(cube).any()
	
	# Set sampling sampleRms for rms measurement
	sampleRms = max(1, int((float(cube.size) / MAX_PIX_CONST)**(1.0 / min(3, len(cube.shape)))))
	
	# Measure noise in original cube with sampling "sampleRms"
	rms = GetRMS(cube, rmsMode=rmsMode, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=verbose, sample=sampleRms)
	
	# Loop over all kernels
	for kernel in kernels:
		[kx, ky, kz, kt] = kernel
		if verbose:
			err.linebreak()
			err.print_progress_time(t0)
			err.message("    Filter {0:} {1:} {2:} {3:} ...".format(kx, ky, kz, kt))
		if kernelUnit == "world" or kernelUnit == "w":
			if verbose: err.message("    Converting filter size to pixels ...")
			kx = abs(float(kx) / header["CDELT1"])
			ky = abs(float(ky) / header["CDELT2"])
			kz = abs(float(kz) / header["CDELT3"])
		if kt == "b":
			if kz != int(math.ceil(kz)) and verbose: err.warning("Rounding width of boxcar z kernel to next integer.")
			kz = int(math.ceil(kz))
		
		# Create a copy of the original cube
		cube_smooth = np.copy(cube)
		
		# Replace all NaNs with zero
		if found_nan:
			cube_smooth[np.isnan(cube)] = 0.0
		
		cube_smooth[(cube_smooth > 0) & (mask > 0)] = maskScaleXY * rms
		cube_smooth[(cube_smooth < 0) & (mask > 0)] = -maskScaleXY * rms
		
		# Spatial smoothing
		if kx + ky:
			cube_smooth = ndimage.filters.gaussian_filter(cube_smooth, [0, ky / FWHM_CONST, kx / FWHM_CONST], mode=edgeMode)
		
		# Spectral smoothing
		if kz:
			if   kt == "b": cube_smooth = ndimage.filters.uniform_filter1d(cube_smooth, kz, axis=0, mode=edgeMode)
			elif kt == "g": cube_smooth = ndimage.filters.gaussian_filter1d(cube_smooth, kz / FWHM_CONST, axis=0, mode=edgeMode)
		
		# Re-insert the NaNs taken out earlier
		if found_nan:
			cube_smooth[np.isnan(cube)] = np.nan
		
		# Per-kernel noise normalisation (Time consuming!)
		if scaleEachSCkernel:
			cube_smooth, noise_smooth = sigma_scale(cube_smooth, scaleX=scaleX, scaleY=scaleY, scaleZ=scaleZ, edgeX=edgeX, edgeY=edgeY, edgeZ=edgeZ, statistic=rmsMode, fluxRange=fluxRange, method=method, windowSpatial=windowSpatial, windowSpectral=windowSpectral, gridSpatial=gridSpatial, gridSpectral=gridSpectral, interpolation=interpolation)

		# Calculate the RMS of the smoothed (possibly normalised) cube
		rms_smooth = GetRMS(cube_smooth, rmsMode=rmsMode, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=verbose, sample=sampleRms)

		# Add pixels above threshold to mask by setting bit 1
		err.message("    Applying +/- {0:} sigma detection threshold".format(threshold))
		with np.errstate(invalid="ignore"):
			mask |= (np.absolute(cube_smooth) >= threshold * rms_smooth)
			#mask = np.bitwise_or(mask, np.greater_equal(np.absolute(cube_smooth), threshold * rms_smooth))
		
		# Delete smoothed cube again
		del cube_smooth
	return
