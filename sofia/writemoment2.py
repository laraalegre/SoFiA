#!/usr/bin/env python
import os
import numpy as np
import scipy.constants
from scipy import interpolate
import astropy.io.fits as pyfits
from sofia import global_settings as glob
from sofia import version
from sofia import error as err



# ==========================
# FUNCTION: Regrid data cube
# ==========================

def regridMaskedChannels(datacube,maskcube,header):
	maskcubeFlt = maskcube.astype("float")
	maskcubeFlt[maskcube > 1] = 1.0
	
	err.message("Regridding...")
	z = (np.arange(1.0, header["naxis3"] + 1) - header["CRPIX3"]) * header["CDELT3"] + header["CRVAL3"]
	
	if header["CTYPE3"] == "VELO-HEL":
		pixscale = (1.0 - header["CRVAL3"] / scipy.constants.c) / (1.0 - z / scipy.constants.c)
	else:
		err.warning("Cannot convert 3rd axis coordinates to frequency.\nIgnoring the effect of CELLSCAL = 1/F.")
		pixscale = np.ones((header["naxis3"]))
	
	x0 = header["crpix1"] - 1
	y0 = header["crpix2"] - 1
	xs = np.arange(datacube.shape[2], dtype=float) - x0
	ys = np.arange(datacube.shape[1], dtype=float) - y0
	
	for zz in range(datacube.shape[0]):
		regrid_channel = interpolate.RectBivariateSpline(ys * pixscale[zz], xs * pixscale[zz], datacube[zz])
		datacube[zz] = regrid_channel(ys, xs)
		regrid_channel_mask = interpolate.RectBivariateSpline(ys * pixscale[zz], xs * pixscale[zz], maskcubeFlt[zz])
		maskcubeFlt[zz] = regrid_channel_mask(ys, xs)
	
	datacube[abs(maskcubeFlt) <= abs(np.nanmin(maskcubeFlt))] = 0.0
	del maskcubeFlt
	
	return datacube



# ===========================================
# FUNCTION: Create nr-of-chan and moment maps
# ===========================================

def writeMoments(datacube, maskcube, filename, debug, header, compress, domom0, domom1, flagOverwrite):
	# Disable moment 2 for now
	domom2 = False
	
	# Exit if nothing to be done
	if not (domom0 or domom1 or domom2):
		err.warning("No moment maps requested; skipping moment map generation.")
		return
	
	# ---------------------------
	# Number of detected channels
	# ---------------------------
	nrdetchan = (maskcube > 0).sum(axis=0)
	if np.nanmax(nrdetchan) < 65535:
		nrdetchan = nrdetchan.astype("int16")
	else:
		nrdetchan = nrdetchan.astype("int32")
	
	hdu = pyfits.PrimaryHDU(data=nrdetchan, header=header)
	hdu.header["BUNIT"] = "Nchan"
	hdu.header["DATAMIN"] = np.nanmin(nrdetchan)
	hdu.header["DATAMAX"] = np.nanmax(nrdetchan)
	hdu.header["ORIGIN"] = version.getVersion(full=True)
	del(hdu.header["CRPIX3"])
	del(hdu.header["CRVAL3"])
	del(hdu.header["CDELT3"])
	del(hdu.header["CTYPE3"])
	
	name = str(filename) + "_nrch.fits"
	if compress: name += ".gz"
	
	# Check for overwrite flag
	if not flagOverwrite and os.path.exists(name):
		err.error("Output file exists: " + str(name) + ".", fatal=False)
	else:
		hdu.writeto(name, output_verify="warn", clobber=True)
	
	# ----------------------
	# Moment 0, 1 and 2 maps
	# ----------------------
	# WARNING: The generation of moment maps will mask the copy of the data cube held
	#          in memory by SoFiA. If you wish to use the original data cube after
	#          this point, please reload it first!
	datacube[maskcube == 0] = 0
	
	# Regrid cube if necessary
	if "CELLSCAL" in header and header["CELLSCAL"] == "1/F":
		err.warning(
			"CELLSCAL keyword with value of 1/F found.\n"
			"Will regrid masked cube before making moment images.")
		datacube = regridMaskedChannels(datacube, maskcube, header)
	
	# ALERT: Why are we doing this?
	#datacube = np.array(datacube, dtype=np.single)
	
	# Extract relevant WCS parameters
	if "CTYPE3" in header and "CDELT3" in header and "CRPIX3" in header and "CRVAL3" in header:
		width = header["CDELT3"]
		chan0 = header["CRPIX3"]
		freq0 = header["CRVAL3"]
		mom_scale_factor = 1.0
		
		# Velocity
		if glob.check_header_keywords(glob.KEYWORDS_VELO, header["CTYPE3"]):
			if not "CUNIT3" in header or header["CUNIT3"].lower() == "m/s":
				# Assuming m/s and converting to km/s
				mom_scale_factor = 1.0e-3
				unit_spec = "km/s"
			elif header["CUNIT3"].lower() == "km/s":
				# Assuming km/s
				unit_spec = "km/s"
			else:
				# Working with whatever velocity units the cube has
				unit_spec = str(header["CUNIT3"])
		# Frequency
		elif glob.check_header_keywords(glob.KEYWORDS_FREQ, header["CTYPE3"]):
			if not "CUNIT3" in header or header["CUNIT3"].lower() == "hz":
				# Assuming Hz
				unit_spec = "Hz"
			elif header["CUNIT3"].lower() == "khz":
				# Assuming kHz and converting to Hz
				mom_scale_factor = 1.0e+3
				unit_spec = "Hz"
			else:
				# Working with whatever frequency units the cube has
				unit_spec = str(header["CUNIT3"])
	else:
		err.warning("Axis descriptors for kinematic axis missing from header.\nMoment maps will not be scaled!")
		width = 1.0
		chan0 = 0.0
		freq0 = 0.0
		mom_scale_factor = 1.0
		unit_spec = "chan"
	
	# Prepare moment maps
	# NOTE: We are making use of NumPy's array broadcasting rules here to avoid
	#       having to cast array sizes to the full 3-D data cube size!
	moments = [None, None, None]
	with np.errstate(invalid="ignore"):
		if domom0 or domom1 or domom2:
			# Definition of moment 0
			moments[0] = np.nansum(datacube, axis=0)
		
		if domom1 or domom2:
			# Definition of moment 1
			velArr = ((np.arange(datacube.shape[0]) + 1.0 - chan0) * width + freq0).reshape((datacube.shape[0], 1, 1))
			moments[1] = np.divide(np.nansum(velArr * datacube, axis=0), moments[0])
		
		if domom2:
			# Definition of moment 2
			velArr = velArr - moments[1]
			moments[2] = np.sqrt(np.divide(np.nansum(velArr * velArr * datacube, axis=0), moments[0]))
	
	# Multiply moment 0 by channel width
	moments[0] *= abs(width)
	
	# Set up unit strings
	if "BUNIT" in header:
		unit_flux = str(header["BUNIT"])
		# Correct for common misspellings of "Jy[/beam]"
		if unit_flux.lower() == "jy":
			unit_flux = "Jy." + unit_spec
		elif unit_flux.lower() == "jy/beam":
			unit_flux = "Jy/beam." + unit_spec
		else:
			unit_flux += "." + unit_spec
	else:
		unit_flux = ""
	unit_mom = [unit_flux, unit_spec, unit_spec]
	
	# Writing moment maps to disk
	for i in range(3):
		if moments[i] is not None:
			err.message("  Writing moment {0:d} image.".format(i))
			moments[i] *= mom_scale_factor
			
			hdu = pyfits.PrimaryHDU(data=moments[i], header=header)
			hdu.header["BUNIT"] = unit_mom[i]
			hdu.header["DATAMIN"] = np.nanmin(moments[i])
			hdu.header["DATAMAX"] = np.nanmax(moments[i])
			hdu.header["ORIGIN"] = version.getVersion(full=True)
			hdu.header["CELLSCAL"] = "CONSTANT"
			del(hdu.header["CRPIX3"])
			del(hdu.header["CRVAL3"])
			del(hdu.header["CDELT3"])
			del(hdu.header["CTYPE3"])
			
			if debug:
				hdu.writeto(str(filename) + "_mom{0:d}.debug.fits".format(i), output_verify="warn", clobber=True)
			else:
				name = str(filename) + "_mom{0:d}.fits".format(i)
				if compress: name += ".gz"
				
				# Check for overwrite flag
				if not flagOverwrite and os.path.exists(name):
					err.error("Output file exists: " + str(name) + ".", fatal=False)
				else:
					hdu.writeto(name, output_verify="warn", clobber=True)
	
	return
