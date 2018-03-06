#! /usr/bin/env python

import numpy as np
from .functions import *
import sys

# Function to read in a cube and scale it by the RMS.
# This script is useful to correct for variation in noise as function of frequency, noisy edges of cubes and channels with strong RFI.

def sigma_scale(cube, scaleX=False, scaleY=False, scaleZ=True, edgeX=0, edgeY=0, edgeZ=0, statistic="mad"):
	verbose = 0
	method_global = True # WARNING This is currently hard-coded to turn local RMS on/off!
	halfwidth = 10
	cadence = 10
	
	# Sigma scaling only works for 3D cubes, as it is mainly designed to correct for differences in frequency.
	
	# cube:       the input cube
	# edgeX,Y,Z:  the edges of the cube that are excluded from the noise calculation (default: 0,0,0)
	# statistic:  "mad" or "std" (default: "mad")
	
	sys.stdout.write("Generating noise-scaled data cube:\n")
	
	if statistic == "mad": sys.stdout.write("Applying Median Absolute Deviation (MAD) statistic.\n")
	if statistic == "std": sys.stdout.write("Applying Standard Deviation (STD) statistic.\n")
	if statistic == "negative": sys.stdout.write("Applying Negative statistic.\n")
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
	
	if method_global:
		# Global noise measurement of entire 2-D plane (fast and memory-friendly)
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
	else:
		# Local noise measurement within running window (slow and memory-intensive)
		
		# Create empty cube (filled with 0) to hold noise values
		rms_cube = np.zeros(cube.shape)
		
		# Iterate over data cube
		it = np.nditer(cube, flags=["multi_index"])
		counter = 0
		rms = 0.0
		while not it.finished:
			if counter % 1000 == 0:
				sys.stdout.write(" Progress: " + str(100 * counter / cube.size) + "%\r")
				sys.stdout.flush()
			
			if counter % cadence == 0:
				z1 = it.multi_index[0] - halfwidth
				z2 = it.multi_index[0] + halfwidth + 1
				y1 = it.multi_index[1] - halfwidth
				y2 = it.multi_index[1] + halfwidth + 1
				x1 = it.multi_index[2] - halfwidth
				x2 = it.multi_index[2] + halfwidth + 1
			
				if z1 < 0: z1 = 0
				if z2 > dimensions[0]: z2 = dimensions[0]
				if y1 < 0: y1 = 0
				if y2 > dimensions[1]: y2 = dimensions[1]
				if x1 < 0: x1 = 0
				if x2 > dimensions[2]: x2 = dimensions[2]
				
				if not np.all(np.isnan(cube[z1:z2, y1:y2, x1:x2])):
					rms = GetRMS(cube[z1:z2, y1:y2, x1:x2], rmsMode=statistic, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
			
			rms_cube[it.multi_index] = rms
			
			counter += 1
			it.iternext()
		
		# Divide data cube by local RMS cube
		cube = rms_cube
		
		# Delete the RMS cube to release its memory
		del rms_cube
	
	sys.stdout.write("Noise-scaled data cube generated.\n\n")
	sys.stdout.flush()
	
	return cube
