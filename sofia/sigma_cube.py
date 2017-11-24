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
	
	if scaleZ:
		z_rms = np.zeros(dimensions[0])
		for i in range(len(z_rms)):
			if np.all(np.isnan(cube[i,y1:y2,x1:x2])):
				z_rms[i] = 0
				sys.stderr.write("WARNING: No valid data found in plane at z = " + str(i + 1) + "\n")
			else:
				z_rms[i] = GetRMS(cube[i, y1:y2, x1:x2], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
		# Scale the cube by the rms
		for i in range(len(z_rms)):
			if z_rms[i] > 0: cube[i,:,:] = cube[i,:,:] / z_rms[i]
	
	if scaleY:
		y_rms = np.zeros(dimensions[1])
		for i in range(len(y_rms)):
			if np.all(np.isnan(cube[z1:z2,i,x1:x2])):
				y_rms[i] = 0
				sys.stderr.write("WARNING: No valid data found in plane at y = " + str(i + 1) + "\n")
			else:
				y_rms[i] = GetRMS(cube[z1:z2, i, x1:x2], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
		# Scale the cube by the rms
		for i in range(len(y_rms)):
			if y_rms[i] > 0: cube[:,i,:] = cube[:,i,:] / y_rms[i]
	
	if scaleX:
		x_rms = np.zeros(dimensions[2])
		for i in range(len(x_rms)):
			if np.all(np.isnan(cube[z1:z2,y1:y2,i])):
				x_rms[i] = 0
				sys.stderr.write("WARNING: No valid data found in plane at x = " + str(i + 1) + "\n")
			else:
				x_rms[i] = GetRMS(cube[z1:z2, y1:y2, i], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
		# Scale the cube by the rms
		for i in range(len(x_rms)):
			if x_rms[i] > 0: cube[:,:,i] = cube[:,:,i] / x_rms[i]
	
	sys.stdout.write('Noise-scaled data cube generated.\n\n')
	sys.stdout.flush()
	return cube
