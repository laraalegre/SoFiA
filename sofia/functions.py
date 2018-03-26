#! /usr/bin/env python

import sys
import math
import numpy as np
import scipy as sp
from distutils.version import StrictVersion, LooseVersion
from sofia import error as err

# Check numpy and scipy version numbers for the nanmedian function import
if LooseVersion(np.__version__) >= LooseVersion("1.9.0"):
	from numpy import nanmedian
elif LooseVersion(sp.__version__) < LooseVersion("0.15.0"):
	from scipy.stats import nanmedian
else:
	from scipy import nanmedian


# ------------------------------
# Gaussian function centred at 0
# ------------------------------
def Gaussian(x, A, sigma):
	err.ensure(sigma != 0, "Invalid width of 0 encountered in Gaussian function.")
	return A * np.exp(-x**2 / (2.0 * sigma**2))


# --------------------------------------------------------
# Second moment of data sample y(x) with first moment of 0
# --------------------------------------------------------
def moment2(x, y):
	err.ensure(x.size == y.size, "Incompatible array sizes encountered in moment calculation.")
	return np.sqrt(np.sum(x * x * y) / np.sum(y))


# ------------------------------------------------
# Standard deviation of data sample with mean of 0
# ------------------------------------------------
def standard_deviation(x):
	err.ensure(x.size > 0, "Array size of 0 encountered in calculation of std. dev.")
	return np.sqrt(np.sum(x * x) / x.size)


# --------------------------------------------------------------------
# Standard deviation of data sample with mean of 0, accounting for NaN
# --------------------------------------------------------------------
def nan_standard_deviation(x):
	y = x[~np.isnan(x)]
	err.ensure(y.size > 0, "Array size of 0 encountered in calculation of std. dev.")
	return np.sqrt(np.sum(y * y) / y.size)


# -----------------------------------
# Function to measure RMS noise level
# -----------------------------------
def GetRMS(cube, rmsMode="negative", fluxRange="all", zoomx=1, zoomy=1, zoomz=1, verbose=0, min_hist_peak=0.05, sample=1, twoPass=True):
	"""
	Description of arguments
	------------------------
	rmsMode    Select which algorithm should be used for calculating the noise.
	           Allowed options:
	             'std'       Standard deviation about 0.
	             'mad'       Median absolute deviation about 0.
	             'moment'    2nd moment of flux histogram, assuming a 1st moment of 0.
	             'gauss'     Width of Gaussian fitted to flux histogram, assuming a centroid of 0.
	             'negative'  Width of Gaussian fitted to negative side of the flux histogram,
	                         again assuming a centroid of 0. This is a legacy option and may be
	                         removed from SoFiA in the future.
	fluxRange  Define which part of the data are to be used in the noise measurement.
	           Allowed options:
	             'negative'  Use only pixels with negative flux.
	             'positive'  Use only pixels with positive flux.
	             'all'       Use both positive and negative (i.e. all) pixels.
	verbose    Print additional progress messages if set to True.
	twoPass    Run a second pass of MAD and STD, this time with a clip level of 5 times
	           the RMS from the first pass.
	"""
	
	# Check input for sanity
	if fluxRange != "all" and fluxRange != "positive" and fluxRange != "negative":
		sys.stderr.write("WARNING: Illegal value of fluxRange = '" + str(fluxRange) + "'.\n")
		sys.stderr.write("         Using default value of 'all' instead.\n")
		fluxRange = "all"
	if rmsMode != "std" and rmsMode != "mad" and rmsMode != "negative" and rmsMode != "gauss" and rmsMode != "moment":
		sys.stderr.write("WARNING: Illegal value of rmsMode = '" + str(rmsMode) + "'.\n")
		sys.stderr.write("         Using default value of 'mad' instead.\n")
		rmsMode = "mad"
	
	# Ensure that we have a 3D cube
	if len(cube.shape) == 2: cube = np.array([cube])
	
	x0, x1 = int(math.ceil((1 - 1.0 / zoomx) * cube.shape[2] / 2)), int(math.floor((1 + 1.0 / zoomx) * cube.shape[2] / 2)) + 1
	y0, y1 = int(math.ceil((1 - 1.0 / zoomy) * cube.shape[1] / 2)), int(math.floor((1 + 1.0 / zoomy) * cube.shape[1] / 2)) + 1
	z0, z1 = int(math.ceil((1 - 1.0 / zoomz) * cube.shape[0] / 2)), int(math.floor((1 + 1.0 / zoomz) * cube.shape[0] / 2)) + 1
	err.print_info("    Estimating rms on subcube (x,y,z zoom = %.0f,%.0f,%.0f) ..." % (zoomx, zoomy, zoomz), verbose)
	err.print_info("    Estimating rms on subcube sampling every %i voxels ..." % (sample), verbose)
	err.print_info("    ... Subcube shape is " + str(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample].shape) + " ...", verbose)
	
	
	# Check if only negative or positive pixels are to be used:
	if fluxRange == "negative":
		with np.errstate(invalid="ignore"):
			halfCube = cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][cube[z0:z1:sample, y0:y1:sample, x0:x1:sample] < 0]
		err.ensure(halfCube.size, "Cannot measure noise from negative flux values.\nNo negative fluxes found in data cube.")
	elif fluxRange == "positive":
		with np.errstate(invalid="ignore"):
			halfCube = cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][cube[z0:z1:sample, y0:y1:sample, x0:x1:sample] > 0]
		err.ensure(halfCube.size, "Cannot measure noise from positive flux values.\nNo positive fluxes found in data cube.")
	# NOTE: The purpose of the with... statement is to temporarily disable certain warnings, as otherwise the
	#       Python interpreter would print a warning whenever a value of NaN is compared to 0. The comparison
	#       is defined to yield False, which conveniently removes NaNs by default without having to do that
	#       manually in a separate step, but the associated warning message is unfortunately a nuisance.
	
	
	# GAUSSIAN FIT TO NEGATIVE FLUXES
	if rmsMode == "negative":
		nrbins = max(100, int(math.ceil(float(cube.size) / 1e+5)))
		
		cubemin = np.nanmin(cube)
		err.ensure(cubemin < 0, "Cannot estimate noise from Gaussian fit to negative flux\nhistogram; no negative fluxes found in data cube.")
		
		bins = np.arange(cubemin, abs(cubemin) / nrbins - 1e-12, abs(cubemin) / nrbins)
		fluxval = (bins[:-1] + bins[1:]) / 2
		rmshisto = np.histogram(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][~np.isnan(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])], bins=bins)[0]
		
		nrsummedbins = 0
		while rmshisto[-nrsummedbins-1:].sum() < min_hist_peak * rmshisto.sum(): nrsummedbins += 1
		
		if nrsummedbins:
			if verbose: sys.stdout.write("    ... adjusting bin size to get a fraction of voxels in central bin >= " + str(min_hist_peak) + "\n")
			nrbins /= (nrsummedbins + 1)
			bins = np.arange(cubemin, abs(cubemin) / nrbins - 1e-12, abs(cubemin) / nrbins)
			fluxval = (bins[:-1] + bins[1:]) / 2.0
			rmshisto = np.histogram(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][~np.isnan(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])], bins=bins)[0]
		
		rms = abs(sp.optimize.curve_fit(Gaussian, fluxval, rmshisto, p0=[rmshisto.max(), -fluxval[rmshisto < rmshisto.max() / 2.0].max() * 2.0 / 2.355])[0][1])
	
	# GAUSSIAN FIT TO FLUX HISTOGRAM / SECOND MOMENT OF FLUX HISTOGRAM
	elif rmsMode == "gauss" or rmsMode == "moment":
		nBins = 100
		dataMin = float(np.nanmin(cube))
		dataMax = float(np.nanmax(cube))
		err.ensure(dataMin < dataMax, "Maximum not greater than minimum. Cannot determine noise level.")
		
		if fluxRange == "negative":
			# Set upper limit to 0
			err.ensure(dataMin < 0.0, "Minimum > 0. Cannot determine noise level for negative pixels.")
			dataMax = 0.0
		elif fluxRange == "positive":
			# Set lower limit to 0
			err.ensure(dataMax > 0.0, "Maximum < 0. Cannot determine noise level for positive pixels.")
			dataMin = 0.0
		else:
			# Select the smallest of the two for both limits
			err.ensure(dataMin < 0.0 and dataMax > 0.0, "Noise values not scattered around 0. Cannot measure noise level.")
			dataMin = -min(abs(dataMin), abs(dataMax))
			dataMax =  min(abs(dataMin), abs(dataMax))
		
		binWidth = (dataMax - dataMin) / float(nBins)
		bins = np.arange(dataMin, dataMax, binWidth)
		binCtr = (bins[:-1] + bins[1:]) / 2.0
		hist = np.histogram(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][~np.isnan(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])], bins=bins)[0]
		
		# Calculate 2nd moment
		mom2 = moment2(binCtr, hist)
		
		# Adjust bin size if necessary
		counter = 0
		while mom2 < 5.0 * binWidth and counter < 2:
			counter += 1
			err.print_info("Increasing number of bins by factor of " + str(int(20.0 * binWidth / mom2)) + " for Gaussian fit.")
			nBins = int(nBins * 20.0 * binWidth / mom2)
			binWidth = (dataMax - dataMin) / float(nBins)
			binCtr = (bins[:-1] + bins[1:]) / 2.0
			hist = np.histogram(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][~np.isnan(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])], bins=bins)[0]
			mom2 = moment2(binCtr, hist)
		
		# Carry out Gaussian fitting if requested
		if rmsMode == "gauss": rms = abs(sp.optimize.curve_fit(Gaussian, binCtr, hist, p0=[hist.max(), mom2])[0][1])
		else: rms = mom2
	
	# MEDIAN ABSOLUTE DEVIATION
	elif rmsMode == "mad":
		if fluxRange == "all":
			# NOTE: Here we assume that the median of the data is zero!
			rms = 1.4826 * nanmedian(abs(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample]), axis=None)
			if twoPass:
				err.print_info("Repeating noise estimation with 5-sigma clip.", verbose)
				with np.errstate(invalid="ignore"):
					rms = 1.4826 * nanmedian(abs(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][abs(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample]) < 5.0 * rms]), axis=None)
		else:
			# NOTE: Here we assume that the median of the data is zero! There are no more NaNs in halfCube.
			rms = 1.4826 * np.median(abs(halfCube), axis=None)
			if twoPass:
				err.print_info("Repeating noise estimation with 5-sigma clip.", verbose)
				rms = 1.4826 * np.median(abs(halfCube[abs(halfCube) < 5.0 * rms]), axis=None)
	
	# STANDARD DEVIATION
	elif rmsMode == "std":
		if fluxRange == "all":
			# NOTE: Here we assume that the mean of the data is zero!
			rms = nan_standard_deviation(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample])
			if twoPass:
				err.print_info("Repeating noise estimation with 5-sigma clip.", verbose)
				with np.errstate(invalid="ignore"):
					rms = nan_standard_deviation(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample][abs(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample]) < 5.0 * rms])
		else:
			# NOTE: Here we assume that the mean of the data is zero! There are no more NaNs in halfCube.
			rms = standard_deviation(halfCube)
			if twoPass:
				err.print_info("Repeating noise estimation with 5-sigma clip.", verbose)
				rms = standard_deviation(halfCube[abs(halfCube) < 5.0 * rms])
	
	err.print_info("    ... %s rms = %.2e (data units)" % (rmsMode, rms), verbose)
	
	return rms


# Copy of previous, full MAD and STD calculations
# -----------------------------------------------
# MAD all: rms = 1.4826 * nanmedian(abs(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample] - nanmedian(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample], axis=None)), axis=None)
# MAD pos/neg: rms = 1.4826 * nanmedian(abs(halfCube - nanmedian(halfCube, axis=None)), axis=None)
# STD all: rms = np.nanstd(cube[z0:z1:sample, y0:y1:sample, x0:x1:sample], axis=None, dtype=np.float64)
# STD pos/neg: rms = np.nanstd(halfCube, axis=None, dtype=np.float64)