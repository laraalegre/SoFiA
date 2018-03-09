#! /usr/bin/env python

import numpy as np
from .functions import *
import sys
import math
from sofia import error as err

"""
Function to read in a cube and scale it by the RMS. This is useful to correct for variation in noise as
a function of frequency, noisy edges of cubes and channels with strong RFI.

Parameters
----------
  cube:            The input cube.
  scaleX,Y,Z:      True or False to choose which axis should be scaled.
  edgeX,Y,Z:       The edges of the cube that are excluded from the noise calculation (default: 0,0,0).
  statistic:       Which algorithm to use: "negative", "mad" or "std" (default: "mad").
  method:          "local"   Measure local rms on specified grid.
                   "global"  Measure RMS in entire planes perpendicular to the axis being scaled.
  windowSpatial:   Spatial window size over which to measure local RMS. Must be even.
  windowSpectral:  Spectral window size over which to measure local RMS. Must be even.
  gridSpatial:     Size of each spatial grid cell for local RMS measurement. Must be even.
  gridSpectral:    Size of each spectral grid cell for local RMS measurement. Must be even.
"""

def sigma_scale(cube, scaleX=False, scaleY=False, scaleZ=True, edgeX=0, edgeY=0, edgeZ=0, statistic="mad", fluxRange="all", method="global", windowSpatial=20, windowSpectral=20, gridSpatial=10, gridSpectral=10):
	# Print some informational messages
	err.print_info("Generating noise-scaled data cube:")
	err.print_info("  Selecting " + str(method) + " noise measurement method.")
	
	if statistic == "mad": err.print_info("  Applying Median Absolute Deviation (MAD) to " + str(fluxRange) + " pixels.")
	if statistic == "std": err.print_info("  Applying Standard Deviation (STD) to " + str(fluxRange) + " pixels.")
	if statistic == "negative": err.print_info("  Applying Gaussian fit to negative pixels.")
	sys.stdout.flush()
	
	# Check the dimensions of the cube (could be obtained from header information)
	dimensions = np.shape(cube)
	
	# LOCAL noise measurement within running window (slower and less memory-friendly)
	if method == "local":
		# Ensure that window and grid sizes are integers greater than 0 and divisible by 2
		windowSpatial  = int(max(windowSpatial, 2))
		windowSpectral = int(max(windowSpectral, 2))
		gridSpatial    = int(max(gridSpatial, 2))
		gridSpectral   = int(max(gridSpectral, 2))
		windowSpatial  += windowSpatial % 2
		windowSpectral += windowSpectral % 2
		gridSpatial    += gridSpatial % 2
		gridSpectral   += gridSpectral % 2
		
		err.print_info("  Using grid size of [" + str(gridSpatial) + ", " + str(gridSpectral) + "]")
		err.print_info("  and window size of [" + str(windowSpatial) + ", " + str(windowSpectral) + "].")
		
		# Create empty cube (filled with 0) to hold noise values
		rms_cube = np.zeros(cube.shape)
		
		# Generate grid points to be used
		gridPointsZ = np.arange((dimensions[0] - gridSpectral * (int(math.ceil(float(dimensions[0]) / float(gridSpectral))) - 1)) / 2, dimensions[0], gridSpectral)
		gridPointsY = np.arange((dimensions[1] - gridSpatial  * (int(math.ceil(float(dimensions[1]) / float(gridSpatial)))  - 1)) / 2, dimensions[1], gridSpatial)
		gridPointsX = np.arange((dimensions[2] - gridSpatial  * (int(math.ceil(float(dimensions[2]) / float(gridSpatial)))  - 1)) / 2, dimensions[2], gridSpatial)
		
		# Divide grid and window sizes by 2 to get radii
		gridSpatial /= 2
		gridSpectral /= 2
		windowSpatial /= 2
		windowSpectral /= 2
		
		# Determine RMS across window centred on grid cell
		for z in gridPointsZ:
			for y in gridPointsY:
				for x in gridPointsX:
					grid = (max(0, z - gridSpectral), min(dimensions[0], z + gridSpectral), max(0, y - gridSpatial), min(dimensions[1], y + gridSpatial), max(0, x - gridSpatial), min(dimensions[2], x + gridSpatial))
					
					window = (max(0, z - windowSpectral), min(dimensions[0], z + windowSpectral), max(0, y - windowSpatial), min(dimensions[1], y + windowSpatial), max(0, x - windowSpatial), min(dimensions[2], x + windowSpatial))
					
					if not np.all(np.isnan(cube[window[0]:window[1], window[2]:window[3], window[4]:window[5]])):
						rms_cube[grid[0]:grid[1], grid[2]:grid[3], grid[4]:grid[5]] = GetRMS(cube[window[0]:window[1], window[2]:window[3], window[4]:window[5]], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
		
		# Divide data cube by RMS cube
		rms_cube[rms_cube <= 0] = np.nan
		cube /= rms_cube
		
		# Delete the RMS cube again to release its memory
		del rms_cube
	
	# GLOBAL noise measurement on entire 2D plane (faster and more memory-friendly)
	else:
		# Define the range over which statistics are calculated
		z1 = int(edgeZ)
		z2 = int(dimensions[0] - edgeZ)
		y1 = int(edgeY)
		y2 = int(dimensions[1] - edgeY)
		x1 = int(edgeX)
		x2 = int(dimensions[2] - edgeX)
		
		# Make sure edges don't exceed cube size
		err.ensure(z1 < z2 and y1 < y2 and x1 < x2, "Edge size exceeds cube size for at least one axis.")
		
		# Measure noise across 2D planes and scale cube accordingly
		if scaleZ:
			for i in range(dimensions[0]):
				if not np.all(np.isnan(cube[i, y1:y2, x1:x2])):
					rms = GetRMS(cube[i, y1:y2, x1:x2], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
					if rms > 0: cube[i, :, :] /= rms
		
		if scaleY:
			for i in range(dimensions[1]):
				if not np.all(np.isnan(cube[z1:z2, i, x1:x2])):
					rms = GetRMS(cube[z1:z2, i, x1:x2], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
					if rms > 0: cube[:, i, :] /= rms
		
		if scaleX:
			for i in range(dimensions[2]):
				if not np.all(np.isnan(cube[z1:z2, y1:y2, i])):
					rms = GetRMS(cube[z1:z2, y1:y2, i], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
					if rms > 0: cube[:, :, i] /= rms
	
	err.print_info("Noise-scaled data cube generated.\n")
	
	return cube
