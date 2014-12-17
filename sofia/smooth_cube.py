#! /usr/bin/env python

from __future__ import division  # make all division floating point. Otherwise is integer division
import string                    # for split, strip string functions
import os                        # for unix/linux system calls (eg 'ls')
import math                      # for maths commands (fabs)
import sys                       # system (sys.exit)
import re                        # for regular expressions
import random                    # for random numbers
import astropy.io.fits as pyfits
import scipy
import numpy
import ast

from numpy import *
from scipy import ndimage



def smooth(indata, kernel, edgeMode, kernelX, kernelY, kernelZ):
	# indata: the input data
	# type: the smoothing type, gaussian (g), boxcar (b), or median (m)
	# edgeMOde: the "mode" parameter which determines how borders are handled: reflect, constant, nearest, mirror or wrap
	
	print 'Start smoothing cube'
	
	if (kernelX + kernelY + kernelZ) == 0:
		sys.stderr.write("WARNING: All the smoothing kernels are set to zero; no smoothing is applied.\n")
		outdata = indata
	else:
		found_nan=numpy.isnan(indata).sum()
		outdata = indata*1
		if found_nan: outdata=numpy.nan_to_num(outdata)
		if type == "g" or kernel == "gaussian":
			print 'The smoothing type is: Gaussian'
			print 'The standard deviation for the Gaussian kernel (x,y,z) is: ', kernelX, kernelY, kernelZ
			print 'The edge is handled as:', edgeMode
			outdata=ndimage.filters.gaussian_filter(indata, sigma=(kernelZ,kernelX,kernelY), mode=edgeMode)
		elif type == "b" or kernel == "boxcar":
			print 'The smoothing type is: Boxcar'
			print 'The size of the uniform filter is (x,y,z) is: ', kernelX, kernelY, kernelZ
			print 'The edge is handled as:', edgeMode
			outdata=ndimage.filters.uniform_filter(indata, size=(kernelZ,kernelX,kernelY), mode=edgeMode)
		elif type == "m" or kernel == "median":
			if kernelX == 0 or kernelY == 0 or kernelZ == 0:
				sys.stderr.write("WARNING: cannot determine median over kernel length zero; no smoothing is applied.\n")
				outdata = indata
			else:
				print 'The smoothing type is: Median'
				print 'The size of the filter (x,y,z) is: ', kernelX, kernelY, kernelZ
				print 'The edge is handled as:', edgeMode
				outdata=ndimage.filters.median_filter(indata, size=(kernelZ,kernelX,kernelY), mode=edgeMode)
		else:
			sys.stderr.write("ERROR: Smoothing type not recognised.\n")
			sys.stderr.write("       The program continues without applying any smoothing.\n")
			#outdata = indata
		if found_nan: outdata[numpy.isnan(indata)]=numpy.nan

	return outdata

