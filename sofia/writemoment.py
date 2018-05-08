#!/usr/bin/env python

import os
import numpy as np
import astropy.io.fits as pyfits
from sofia import statistics as stat
from sofia import functions as func
from sofia import error as err
from sofia import __version_full__ as sofia_version_full
from sofia import __astropy_arg_overwrite__



# ======================================
# FUNCTION: Create Nchan and moment maps
# ======================================

def writeMoments(datacube, maskcube, filename, debug, header, compress, write_mom, flagOverwrite):
	# Exit if nothing is to be done
	if not any(write_mom):
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
	hdu.header["ORIGIN"] = sofia_version_full
	func.delete_header(hdu.header, "CTYPE3")
	func.delete_header(hdu.header, "CRPIX3")
	func.delete_header(hdu.header, "CRVAL3")
	func.delete_header(hdu.header, "CDELT3")
	
	name = str(filename) + "_nrch.fits"
	if compress: name += ".gz"
	
	# Check for overwrite flag
	if func.check_overwrite(name, flagOverwrite):
		hdu.writeto(name, output_verify="warn", **__astropy_arg_overwrite__)
	
	
	# ----------------------
	# Moment 0, 1 and 2 maps
	# ----------------------
	# WARNING: The generation of moment maps will mask the copy of the data cube held
	#          in memory by SoFiA. If you wish to use the original data cube after
	#          this point, please reload it first!
	datacube[maskcube == 0] = 0
	
	# Regrid cube if necessary
	# ------------------------
	if "CELLSCAL" in header and header["CELLSCAL"] == "1/F":
		err.warning(
			"CELLSCAL keyword with value of 1/F found.\n"
			"Will regrid data cube before creating moment images.")
		datacube = func.regridMaskedChannels(datacube, maskcube, header)
	
	# ALERT: Why are we doing this?
	#datacube = np.array(datacube, dtype=np.single)
	
	# Extract relevant WCS parameters
	# -------------------------------
	if func.check_wcs_info(header):
		width = header["CDELT3"]
		chan0 = header["CRPIX3"]
		freq0 = header["CRVAL3"]
		mom_scale_factor = 1.0
		
		# Velocity
		if func.check_header_keywords(func.KEYWORDS_VELO, header["CTYPE3"]):
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
		elif func.check_header_keywords(func.KEYWORDS_FREQ, header["CTYPE3"]):
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
		err.warning("Axis descriptors missing from FITS file header.\nMoment maps will not be scaled!")
		width = 1.0
		chan0 = 0.0
		freq0 = 0.0
		mom_scale_factor = 1.0
		unit_spec = "chan"
	
	# Calculate moment maps
	# ---------------------
	moments = [None, None, None]
	# ALERT: The order of axes and memory layout differ, so transposition
	#        and memory layout change are necessary! This will need to be
	#        fixed in some way.
	datacube = np.transpose(datacube, axes=[1, 2, 0]).astype(np.float32, order="C")
	with np.errstate(invalid="ignore"):
		if any(write_mom):
			# Definition of moment 0
			#moments[0] = np.nansum(datacube, axis=0)
			moments[0] = stat.moment(datacube, mom=0)
		
		if write_mom[1] or write_mom[2]:
			# Definition of moment 1
			#velArr = ((np.arange(datacube.shape[0]) + 1.0 - chan0) * width + freq0).reshape((datacube.shape[0], 1, 1))
			#moments[1] = np.divide(np.nansum(velArr * datacube, axis=0), moments[0])
			moments[1] = stat.moment(datacube, mom=1, mom0=moments[0])
		
		if write_mom[2]:
			# Definition of moment 2
			#velArr = velArr - moments[1]
			#moments[2] = np.sqrt(np.divide(np.nansum(velArr * velArr * datacube, axis=0), moments[0]))
			moments[2] = stat.moment(datacube, mom=2, mom0=moments[0], mom1=moments[1])
	datacube = np.transpose(datacube, axes=[2, 0, 1])
	
	# Convert moments to physical units
	# ---------------------------------
	if moments[0] is not None: moments[0] *= abs(width)
	if moments[1] is not None: moments[1] = (moments[1] + 1.0 - chan0) * width + freq0
	if moments[2] is not None: moments[2] *= abs(width)
	
	# Set up unit strings
	# -------------------
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
		err.warning("Cannot determine flux unit; BUNIT missing from header.")
		unit_flux = ""
	unit_mom = [unit_flux, unit_spec, unit_spec]
	
	# Writing moment maps to disk
	# ---------------------------
	for i in range(3):
		if write_mom[i] and moments[i] is not None:
			err.message("Writing moment {0:d} image.".format(i))
			moments[i] *= mom_scale_factor
			
			hdu = pyfits.PrimaryHDU(data=moments[i], header=header)
			hdu.header["BUNIT"] = unit_mom[i]
			hdu.header["DATAMIN"] = np.nanmin(moments[i])
			hdu.header["DATAMAX"] = np.nanmax(moments[i])
			hdu.header["ORIGIN"] = sofia_version_full
			hdu.header["CELLSCAL"] = "CONSTANT"
			func.delete_header(hdu.header, "CRPIX3")
			func.delete_header(hdu.header, "CRVAL3")
			func.delete_header(hdu.header, "CDELT3")
			func.delete_header(hdu.header, "CTYPE3")
			
			if debug:
				hdu.writeto(str(filename) + "_mom{0:d}.debug.fits".format(i), output_verify="warn", **__astropy_arg_overwrite__)
			else:
				name = str(filename) + "_mom{0:d}.fits".format(i)
				if compress: name += ".gz"
				
				# Check for overwrite flag
				if func.check_overwrite(name, flagOverwrite):
					hdu.writeto(name, output_verify="warn", **__astropy_arg_overwrite__)
	
	return
