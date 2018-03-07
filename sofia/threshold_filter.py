#! /usr/bin/env python

# import default python libraries
import numpy as np
import os
from .functions import *

# Run a simple threshold filter and write out mask:

def filter(mask, cube, header, clipMethod, threshold, rmsMode, fluxRange, verbose):
	if clipMethod == 'relative':
		# Determine the clip level
		# Measure noise in original cube
		rms = GetRMS(cube, rmsMode=rmsMode, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
		print ('Estimated rms = ' + str(rms))
		clip = threshold * rms
	if clipMethod == 'absolute':
		clip = threshold
	print ('Using clip threshold: ' + str(clip))

	# Check whether there are NaNs:
	nan_mask = np.isnan(cube)
	found_nan = nan_mask.sum()
	
	# Set NaNs to zero (and INFs to a finite number) if necessary:
	if found_nan: cube = np.nan_to_num(cube)
	
	# Run the threshold finder, setting bit 1 of the mask for |cube| > clip:
	np.bitwise_or(mask, np.greater_equal(np.absolute(cube), clip), mask)
	
	# Put NaNs (but not INFs) back into the data cube:
	if found_nan: cube[nan_mask] = np.nan
	
	return
