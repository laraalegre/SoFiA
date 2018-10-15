#! /usr/bin/env python

import math
import numpy as np
from scipy.ndimage.morphology import binary_dilation
from sofia.functions import GetRMS
from sofia import error as err

"""
Function to automatically mask poor regions of the cube. At the moment it only works on LOS's. It masks LOS's with anomalously high RMS.

Parameters
----------
  cube:            The input cube.
  mskthr:          Threshold for masking LOS's with anomalously high RMS. A LOS is masked if its RMS is mskthr*STD above the typical RMS of the cube,\
                   where STD is the width of the distribution of RMS values across all LOS's.
  dilside:         Side of the box structuring element used to propagate the flagging to neighbouring LOS's.
"""

def maskerrors(cube,mskthr,dilside):
	rms0 = GetRMS(cube, rmsMode="mad", fluxRange="negative", zoomx=1, zoomy=1, zoomz=1, verbose=0, twoPass=True)

	# COMPACT ...
	#los_rms = np.nanstd(cube,axis=0)
	#los_rms = 1.4826*np.nanmedian(abs(cube),axis=0)

	# ... OR NOT in order to use the TwoPass option of GetRMS
	los_rms=np.zeros((cube.shape[-2],cube.shape[-1]))
	for xx in range(cube.shape[-1]):
		for yy in range(cube.shape[-2]): los_rms[yy,xx]=GetRMS(cube[:,yy:yy+1,xx:xx+1],rmsMode="mad", fluxRange="all", twoPass=True)

	# Mask all LOS's whose RMS is > mskthr*STD above rms0 (optionally extended to neighbouring LOS's using binary dilation with a box structuring element)
	los_rms_disp=np.nanstd(los_rms)
	los_rms=(los_rms<rms0+mskthr*los_rms_disp)
	los_rms=binary_dilation(~los_rms,structure=np.ones((dilside,dilside)))
	cube[:,los_rms]=np.nan

	return cube
