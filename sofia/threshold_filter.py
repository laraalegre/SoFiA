#! /usr/bin/env python

# import default python libraries
import numpy as np
from sofia.functions import GetRMS
from sofia import error as err

# Run a simple threshold filter and write out mask:

def filter(mask, cube, header, clipMethod, threshold, rmsMode, fluxRange, verbose):
	err.message("Running threshold finder.")
	
	# Sanity checks of user input
	err.ensure(
		clipMethod in {"absolute", "relative"},
		"Threshold finder failed. Illegal clip method: '" + str(clipMethod) + "'.")
	err.ensure(
		rmsMode in {"std", "mad", "gauss", "negative"},
		"Threshold finder failed. Illegal RMS mode: '" + str(rmsMode) + "'.")
	err.ensure(
		fluxRange in {"positive", "negative", "all"},
		"Threshold finder failed. Illegal flux range: '" + str(fluxRange) + "'.")
	
	# Scale threshold by RMS if requested
	if clipMethod == "relative":
		threshold *= GetRMS(cube, rmsMode=rmsMode, fluxRange=fluxRange, zoomx=1, zoomy=1, zoomz=1, verbose=verbose)
	
	# Print some information and check sign of threshold
	err.message("  Using threshold of " + str(threshold) + ".")
	err.ensure(threshold >= 0.0, "Threshold finder failed. Threshold value is negative.")
	
	# Run the threshold finder, setting bit 1 of the mask for |cube| >= |threshold|:
	np.bitwise_or(mask, np.greater_equal(np.absolute(cube), threshold), out=mask)
	
	return
