#! /usr/bin/env python

import astropy.io.fits as pf
import numpy as np
import math as mt
import scipy
from scipy import optimize
import scipy.ndimage as nd
from sys import argv
import string
from os import path
import sys
from functions import *
from time import time


def GaussianNoise(F, N0, s0):
	return N0 * np.exp(-F**2 / 2 / s0**2)


def SizeFilter(mskt, sfx, sfy, sfz, sbx, sby, sbz, zt, sizeFilter, edgeMode='constant', verbose=0):
	mskt = nd.filters.gaussian_filter(mskt, [0, mt.sqrt(sfy**2 + sby**2) / 2.355, mt.sqrt(sfx**2 + sbx**2) / 2.355], mode=edgeMode)
	if zt == 'b': mskt = nd.filters.uniform_filter1d(mskt, max(sbz, sfz), axis=0, mode=edgeMode)
	elif zt == 'g': mskt = nd.filters.gaussian_filter1d(mskt, max(sbz, sfz / 2.355), axis=0, mode=edgeMode)
	mskt[mskt< sizeFilter] = 0
	mskt[mskt>=sizeFilter] = 1
	return mskt


def SortKernels(kernels):
	# Sorting kernels
	uniquesky = []
	velsmooth = []
	velfshape = []
	
	for jj in np.array(kernels):
		if list(jj[:2].astype(float)) not in uniquesky: uniquesky.append(list(jj[:2].astype(float)))
	uniquesky = [kk[1] for kk in sorted([(float(jj[0]), jj) for jj in uniquesky])]
	
	for jj in uniquesky:
		velsmooth.append([])
		velfshape.append([])
		for ii in np.array(kernels):
			if list(ii[:2].astype(float)) == jj:
				velsmooth[-1].append(int(ii[2]))
				velfshape[-1].append(ii[3])
	return uniquesky, velsmooth, velfshape


def SCfinder_mem(cube, header, t0, kernels=[[0, 0, 0, 'b'],], threshold=3.5, sizeFilter=0, maskScaleXY=2.0, maskScaleZ=2.0, kernelUnit='pixel', edgeMode='constant', rmsMode='negative', verbose=0):
	# Create binary mask array
	msk = np.zeros(cube.shape, np.bool)
	found_nan = np.isnan(cube).sum()
	
	# Set sampling sampleRms for rms measurement
	maxNrVox = 1e+6 # maximum nr of voxels for noise calculation; sampling is set accordingly
	sampleRms = max(1, int((float(np.array(cube.shape).prod()) / maxNrVox)**(1.0 / min(3, len(cube.shape)))))
	
	# Measure noise in original cube with sampling "sampleRms"
	rms = GetRMS(cube, rmsMode=rmsMode, zoomx=1, zoomy=1, zoomz=1, verbose=verbose, sample=sampleRms)
	
	# Loop over all kernels
	for jj in kernels:
		[kx, ky, kz, kt] = jj
		if verbose:
			print '\n--- %.3f seconds since start' % (time() - t0)
			print '    Filter %s %s %s %s ...' % (kx, ky, kz, kt)
		if kernelUnit == 'world' or kernelUnit == 'w':
			if verbose: print '    Converting filter size to pixels ...'
			kx = abs(float(kx) / header['CDELT1'])
			ky = abs(float(ky) / header['CDELT2'])
			kz = abs(float(kz) / header['CDELT3'])
		if kt == 'b':
			if kz != int(mt.ceil(kz)) and verbose: sys.stderr.write('    WARNING: Rounding width of boxcar z kernel to next integer')
			kz = int(mt.ceil(kz))
		sys.stdout.flush()
		
		# Create a copy of the original cube:
		smoothedcube = cube * 1.0
		
		# Replace all NaNs with zero (and INFs with a finite number):
		if found_nan: smoothedcube = np.nan_to_num(smoothedcube)
		
		smoothedcube[(smoothedcube > 0) & (msk > 0)] = +maskScaleXY * rms
		smoothedcube[(smoothedcube < 0) & (msk > 0)] = -maskScaleXY * rms
		
		# Spatial smoothing:
		if kx + ky: smoothedcube = nd.filters.gaussian_filter(smoothedcube, [0, ky / 2.355, kx / 2.355], mode=edgeMode)
		
		# Spectral smoothing:
		if kz:
			if kt == 'b': smoothedcube = nd.filters.uniform_filter1d(smoothedcube, kz, axis=0, mode=edgeMode)
			elif kt == 'g': smoothedcube = nd.filters.gaussian_filter1d(smoothedcube, kz / 2.355, axis=0, mode=edgeMode)
		
		# Re-insert the NaNs (but not the INFs) taken out earlier:
		if found_nan: smoothedcube[np.isnan(cube)] = np.nan
		
		# Calculate the RMS of the smoothed cube:
		smoothedrms = GetRMS(smoothedcube, rmsMode=rmsMode, zoomx=1, zoomy=1, zoomz=1, verbose=verbose, sample=sampleRms)
		
		# Get rid of the NaNs a second time:
		if found_nan: smoothedcube = np.nan_to_num(smoothedcube)
		
		# Add pixels above threshold to mask by setting bit 1:
		msk = np.bitwise_or(msk, np.greater_equal(np.absolute(smoothedcube), threshold * smoothedrms), msk)
		
		# Delete smoothed cube again:
		del(smoothedcube)
	return msk
