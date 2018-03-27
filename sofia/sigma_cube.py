#! /usr/bin/env python

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
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
  interpolation:   Interpolate values between grid points? Can be "none", "linear" or "cubic"
"""

def sigma_scale(cube, scaleX=False, scaleY=False, scaleZ=True, edgeX=0, edgeY=0, edgeZ=0, statistic="mad", fluxRange="all", method="global", windowSpatial=20, windowSpectral=20, gridSpatial=0, gridSpectral=0, interpolation="none"):
	# Print some informational messages
	err.print_info("Generating noise-scaled data cube:")
	err.print_info("  Selecting " + str(method) + " noise measurement method.")
	
	if statistic == "mad": err.print_info("  Applying Median Absolute Deviation (MAD) to " + str(fluxRange) + " pixels.")
	if statistic == "std": err.print_info("  Applying Standard Deviation (STD) to " + str(fluxRange) + " pixels.")
	if statistic == "negative": err.print_info("  Applying Gaussian fit to negative pixels.")
	sys.stdout.flush()
	
	# Check the dimensions of the cube (could be obtained from header information)
	dimensions = np.shape(cube)
	
	# Create empty cube (filled with 0) to hold noise values
	rms_cube = np.full(cube.shape, np.nan)
	
	# LOCAL noise measurement within running window (slower and less memory-friendly)
	if method == "local":
		# Make window sizes integers >= 1
		windowSpatial = max(int(windowSpatial), 1)
		windowSpectral = max(int(windowSpectral), 1)
		
		# Ensure that window sizes are odd
		windowSpatial += (1 - windowSpatial % 2)
		windowSpectral += (1 - windowSpectral % 2)
		
		# Set grid sizes to half the window sizes if undefined
		if not gridSpatial: gridSpatial = windowSpatial // 2
		if not gridSpectral: gridSpectral = windowSpectral // 2
		
		# Make grid sizes integers >= 1
		gridSpatial = max(int(gridSpatial), 1)
		gridSpectral = max(int(gridSpectral), 1)
		
		# Ensure that grid sizes are odd
		gridSpatial += (1 - gridSpatial % 2)
		gridSpectral += (1 - gridSpectral % 2)
		
		# Print grid and window sizes adopted
		err.print_info("  Using grid size of [" + str(gridSpatial) + ", " + str(gridSpectral) + "]")
		err.print_info("  and window size of [" + str(windowSpatial) + ", " + str(windowSpectral) + "].")
		
		# Generate grid points to be used
		gridPointsZ = np.arange((dimensions[0] - gridSpectral * (int(math.ceil(float(dimensions[0]) / float(gridSpectral))) - 1)) // 2, dimensions[0], gridSpectral)
		gridPointsY = np.arange((dimensions[1] - gridSpatial  * (int(math.ceil(float(dimensions[1]) / float(gridSpatial)))  - 1)) // 2, dimensions[1], gridSpatial)
		gridPointsX = np.arange((dimensions[2] - gridSpatial  * (int(math.ceil(float(dimensions[2]) / float(gridSpatial)))  - 1)) // 2, dimensions[2], gridSpatial)
		
		# Divide grid and window sizes by 2 to get radii
		radiusGridSpatial = gridSpatial // 2
		radiusGridSpectral = gridSpectral // 2
		radiusWindowSpatial = windowSpatial // 2
		radiusWindowSpectral = windowSpectral // 2
		
		# Determine RMS across window centred on grid cell
		for z in gridPointsZ:
			for y in gridPointsY:
				for x in gridPointsX:
					grid = (max(0, z - radiusGridSpectral), min(dimensions[0], z + radiusGridSpectral + 1), max(0, y - radiusGridSpatial), min(dimensions[1], y + radiusGridSpatial + 1), max(0, x - radiusGridSpatial), min(dimensions[2], x + radiusGridSpatial + 1))
					
					window = (max(0, z - radiusWindowSpectral), min(dimensions[0], z + radiusWindowSpectral + 1), max(0, y - radiusWindowSpatial), min(dimensions[1], y + radiusWindowSpatial + 1), max(0, x - radiusWindowSpatial), min(dimensions[2], x + radiusWindowSpatial + 1))
					
					if not np.all(np.isnan(cube[window[0]:window[1], window[2]:window[3], window[4]:window[5]])):
						if interpolation == "linear" or interpolation == "cubic":
							# Write value into grid point for later interpolation
							rms_cube[z, y, x] = GetRMS(cube[window[0]:window[1], window[2]:window[3], window[4]:window[5]], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
						else:
							# Fill entire grid cell
							rms_cube[grid[0]:grid[1], grid[2]:grid[3], grid[4]:grid[5]] = GetRMS(cube[window[0]:window[1], window[2]:window[3], window[4]:window[5]], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
		
		# Carry out interpolation if requested, taking NaNs into account
		if interpolation == "linear" or interpolation == "cubic":
			err.print_info("Interpolating in between grid points.")
			
			# First across each spatial plane
			if gridSpatial > 1:
				for z in gridPointsZ:
					for y in gridPointsY:
						data_values   = rms_cube[z, y, gridPointsX]
						not_nan = np.logical_not(np.isnan(data_values))
						if any(not_nan):
							interp_coords = np.arange(0, dimensions[2])
							if interpolation == "cubic":
								spline = InterpolatedUnivariateSpline(gridPointsX[not_nan], data_values[not_nan])
								rms_cube[z, y, 0:dimensions[2]] = spline(interp_coords)
							else:
								interp_values = np.interp(interp_coords, gridPointsX[not_nan], data_values[not_nan])
								rms_cube[z, y, 0:dimensions[2]] = interp_values
					for x in range(dimensions[2]):
						data_values   = rms_cube[z, gridPointsY, x]
						not_nan = np.logical_not(np.isnan(data_values))
						if any(not_nan):
							interp_coords = np.arange(0, dimensions[1])
							if interpolation == "cubic":
								spline = InterpolatedUnivariateSpline(gridPointsY[not_nan], data_values[not_nan])
								rms_cube[z, 0:dimensions[1], x] = spline(interp_coords)
							else:
								interp_values = np.interp(interp_coords, gridPointsY[not_nan], data_values[not_nan])
								rms_cube[z, 0:dimensions[1], x] = interp_values
					# Alternative option: 2-D spatial interpolation using SciPy's interp2d
					#from scipy.interpolate import interp2d
					#xx, yy = np.meshgrid(gridPointsX, gridPointsY)
					#data_values = rms_cube[z, yy, xx]
					#f = interp2d(gridPointsX, gridPointsY, data_values, kind="cubic")
					#interp_coords_x = np.arange(0, dimensions[2])
					#interp_coords_y = np.arange(0, dimensions[1])
					#rms_cube[z, :, :] = f(interp_coords_x, interp_coords_y)
			
			# Then along the spectral axis
			if gridSpectral > 1:
				for y in range(dimensions[1]):
					for x in range(dimensions[2]):
						data_values   = rms_cube[gridPointsZ, y, x]
						not_nan = np.logical_not(np.isnan(data_values))
						if any(not_nan):
							interp_coords = np.arange(0, dimensions[0])
							if interpolation == "cubic":
								spline = InterpolatedUnivariateSpline(gridPointsZ[not_nan], data_values[not_nan])
								rms_cube[0:dimensions[0], y, x] = spline(interp_coords)
							else:
								interp_values = np.interp(interp_coords, gridPointsZ[not_nan], data_values[not_nan])
								rms_cube[0:dimensions[0], y, x] = interp_values
		
		# Replace any invalid RMS values with NaN
		with np.errstate(invalid="ignore"):
			rms_cube[rms_cube <= 0] = np.nan
		
		# Divide data cube by RMS cube
		cube /= rms_cube
		
		# Delete the RMS cube again to release its memory
		#del rms_cube
	
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
					if rms > 0:
						rms_cube[i, :, :] = rms
						cube[i, :, :] /= rms
		
		if scaleY:
			for i in range(dimensions[1]):
				if not np.all(np.isnan(cube[z1:z2, i, x1:x2])):
					rms = GetRMS(cube[z1:z2, i, x1:x2], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
					if rms > 0:
						rms_cube[:, i, :] = rms
						cube[:, i, :] /= rms
		
		if scaleX:
			for i in range(dimensions[2]):
				if not np.all(np.isnan(cube[z1:z2, y1:y2, i])):
					rms = GetRMS(cube[z1:z2, y1:y2, i], rmsMode=statistic, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=0)
					if rms > 0:
						rms_cube[:, :, i] = rms
						cube[:, :, i] /= rms
	
	err.print_info("Noise-scaled data cube generated.\n")
	
	return cube, rms_cube
