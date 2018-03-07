#! /usr/bin/env python

import sys
import math
import numpy as np
import scipy as sp
from distutils.version import StrictVersion, LooseVersion

# Check numpy and scipy version numbers for the nanmedian function import
if LooseVersion(np.__version__) >= LooseVersion("1.9.0"):
	from numpy import nanmedian
elif LooseVersion(sp.__version__) < LooseVersion("0.15.0"):
	from scipy.stats import nanmedian
else:
	from scipy import nanmedian


def Gaussian(x, A, sigma):
	return A * np.exp(-x**2 / (2.0 * sigma**2))


def GetRMS(cube, rmsMode="negative", zoomx=1, zoomy=1, zoomz=1, verbose=0, min_hist_peak=0.05, sample=1, fluxRange="all"):
	"""
	Description of arguments
	------------------------
	   fluxRange  Define which part of the data are to be used in the noise measurement.
	              Allowed values:
	                'negative'  Use only pixels with negative flux.
	                'positive'  Use only pixels with positive flux.
	                'all'       Use both positive and negative (i.e. all) pixels.
	"""
	
	# Check input for sanity
	if fluxRange != "all" and fluxRange != "positive" and fluxRange != "negative":
		sys.stderr.write("WARNING: Illegal value of fluxRange = '" + str(fluxRange) + "'.\n")
		sys.stderr.write("         Using default value of 'all' instead.\n")
		fluxRange = "all"
	if rmsMode != "std" and rmsMode != "mad" and rmsMode != "negative":
		sys.stderr.write("WARNING: Illegal value of rmsMode = '" + str(rmsMode) + "'.\n")
		sys.stderr.write("         Using default value of 'mad' instead.\n")
		rmsMode = "mad"
	
	# Ensure that we have a 3D cube
	if len(cube.shape) == 2: cube = np.array([cube])
	
	x0, x1 = int(math.ceil((1 - 1.0 / zoomx) * cube.shape[2] / 2)), int(math.floor((1 + 1.0 / zoomx) * cube.shape[2] / 2)) + 1
	y0, y1 = int(math.ceil((1 - 1.0 / zoomy) * cube.shape[1] / 2)), int(math.floor((1 + 1.0 / zoomy) * cube.shape[1] / 2)) + 1
	z0, z1 = int(math.ceil((1 - 1.0 / zoomz) * cube.shape[0] / 2)), int(math.floor((1 + 1.0 / zoomz) * cube.shape[0] / 2)) + 1
	if verbose: sys.stdout.write("    Estimating rms on subcube (x,y,z zoom = %.0f,%.0f,%.0f) ..." % (zoomx, zoomy, zoomz))
	if verbose: sys.stdout.write("    Estimating rms on subcube sampling every %i voxels ..." % (sample))
	if verbose: sys.stdout.write("    ... Subcube shape is " + str(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample].shape) + " ...")
	
	
	# Check if only negative or positive pixels are to be used:
	if fluxRange == "negative":
		subCube = cube[cube[z0:z1:sample, y0:y1:sample, x0:x1:sample] < 0]
		subCube[::2] *= -1   # Flip the sign of every other element
	elif fluxRange == "positive":
		subCube = cube[cube[z0:z1:sample, y0:y1:sample, x0:x1:sample] > 0]
		subCube[::2] *= -1   # Flip the sign of every other element
	
	
	# GAUSSIAN FIT TO NEGATIVE FLUXES
	if rmsMode == "negative":
		nrbins = max(100, int(math.ceil(float(cube.size) / 1e+5)))
		cubemin = np.nanmin(cube)
		bins = np.arange(cubemin, abs(cubemin) / nrbins - 1e-12, abs(cubemin) / nrbins)
		fluxval = (bins[:-1] + bins[1:]) / 2
		rmshisto = np.histogram(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][~np.isnan(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])], bins=bins)[0]
		
		nrsummedbins = 0
		while rmshisto[-nrsummedbins-1:].sum() < min_hist_peak * rmshisto.sum(): nrsummedbins += 1
		
		if nrsummedbins:
			if verbose: sys.stdout.write("    ... adjusting bin size to get a fraction of voxels in central bin >= " + str(min_hist_peak))
			nrbins /= (nrsummedbins + 1)
			bins = np.arange(cubemin, abs(cubemin) / nrbins - 1e-12, abs(cubemin) / nrbins)
			fluxval = (bins[:-1] + bins[1:]) / 2
			rmshisto = np.histogram(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][~np.isnan(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])], bins=bins)[0]
		
		rms = abs(sp.optimize.curve_fit(Gaussian, fluxval, rmshisto, p0=[rmshisto.max(), -fluxval[rmshisto < rmshisto.max() / 2].max() * 2 / 2.355])[0][1])
	
	# MEDIAN ABSOLUTE DEVIATION
	elif rmsMode == "mad":
		if fluxRange == "all":
			rms = 1.4826 * nanmedian(abs(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample] - nanmedian(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample], axis=None)), axis=None)
		else:
			rms = 1.4826 * nanmedian(abs(subCube - nanmedian(subCube, axis=None)), axis=None)
	
	# STANDARD DEVIATION
	elif rmsMode == "std":
		if fluxRange == "all":
			rms = np.nanstd(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample], axis=None, dtype=np.float64)
		else:
			rms = np.nanstd(subCube, axis=None, dtype=np.float64)
	
	if verbose: sys.stdout.write("    ... %s rms = %.2e (data units)" % (rmsMode, rms))
	
	return rms
