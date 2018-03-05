#! /usr/bin/env python

import numpy as np
from .functions import *
import sys

# Function to read in a cube and scale it by the RMS.
# This script is useful to correct for variation in noise as function of frequency, noisy edges of cubes and channels with strong RFI.

def sigma_scale(cube, scaleX=False, scaleY=False, scaleZ=True, edgeX=0, edgeY=0, edgeZ=0, statistic="mad"):
	verbose = 0
	
	# Sigma scaling only works for 3D cubes, as it is mainly designed to correct for differences in frequency.
	
	# cube:       the input cube
	# edgeX,Y,Z:  the edges of the cube that are excluded from the noise calculation (default: 0,0,0)
	# statistic:  "mad" or "std" (default: "mad")
	
	sys.stdout.write('Generating noise-scaled data cube:\n')
	
	if statistic == 'mad': sys.stdout.write('Applying Median Absolute Deviation (MAD) statistic.\n')
	if statistic == 'std': sys.stdout.write('Applying Standard Deviation (STD) statistic.\n')
	if statistic == 'negative': sys.stdout.write('Applying Negative statistic.\n')
	sys.stdout.flush()
	
	# Check the dimensions of the cube (could be obtained from header information)
	dimensions = np.shape(cube)
	
	# Define the range over which statistics are calculated
	z1 = edgeZ
	z2 = dimensions[0] - edgeZ
	y1 = edgeY
	y2 = dimensions[1] - edgeY
	x1 = edgeX
	x2 = dimensions[2] - edgeX
	
	if z1 >= z2 or y1 >= y2 or x1 >= x2:
		sys.stderr.write("ERROR: Edge size exceeds cube size for at least one axis.\n")
		sys.exit(1)
	
	if scaleZ:
		for i in range(dimensions[0]):
			if not np.all(np.isnan(cube[i, y1:y2, x1:x2])):
				rms = GetRMS(cube[i, y1:y2, x1:x2], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
				if rms > 0: cube[i, :, :] /= rms
	
	if scaleY:
		for i in range(dimensions[1]):
			if not np.all(np.isnan(cube[z1:z2, i, x1:x2])):
				rms = GetRMS(cube[z1:z2, i, x1:x2], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
				if rms > 0: cube[:, i, :] /= rms
	
	if scaleX:
		for i in range(dimensions[2]):
			if not np.all(np.isnan(cube[z1:z2, y1:y2, i])):
				rms = GetRMS(cube[z1:z2, y1:y2, i], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
				if rms > 0: cube[:, :, i] /= rms
	
	sys.stdout.write('Noise-scaled data cube generated.\n\n')
	sys.stdout.flush()
	return cube
