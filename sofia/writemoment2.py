#!/usr/bin/env python
import os
import numpy as np
import scipy.constants
from scipy import interpolate
import astropy.io.fits as pyfits
from sofia import global_settings as glob
from sofia import version
from sofia import error as err



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


def writeMoments(datacube, maskcube, filename, debug, header, compress, domom0, domom1, flagOverwrite):
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
	
	name = "%s_nrch.fits" % filename
	if compress: name += ".gz"
	
	# Check for overwrite flag
	if not flagOverwrite and os.path.exists(name):
		err.error("Output file exists: " + str(name) + ".", fatal=False)
	else:
		hdu.writeto(name, output_verify="warn", clobber=True)
	
	# WARNING: The generation of moment maps will mask the copy of the data cube held
	#          in memory by SoFiA. If you wish to use the original data cube after
	#          this point, please reload it first!
	datacube[maskcube == 0] = 0
	
	if "CELLSCAL" in header and header["CELLSCAL"] == "1/F":
		err.warning(
			"CELLSCAL keyword with value of 1/F found.\n"
			"Will regrid masked cube before making moment images.")
		datacube = regridMaskedChannels(datacube, maskcube, header)
	
	datacube = np.array(datacube, dtype=np.single)
	
	# --------------
	# Moment 0 image
	# --------------
	if domom0 or domom1:
		# Calculate moment 0
		m0 = np.nansum(datacube, axis=0)
	
	if domom0:
		err.message("Writing moment-0") # in units of header["bunit"]*header["CDELT3"]
		
		# Velocity
		if glob.check_values(glob.KEYWORDS_VELO, header["CTYPE3"]):
			if not "CUNIT3" in header or header["CUNIT3"].lower() == "m/s":
				# Converting (assumed) m/s to km/s
				dkms = abs(header["CDELT3"]) / 1e+3
				scalemom12 = 1.0 / 1e+3
				bunitExt = ".km/s"
			elif header["CUNIT3"].lower() == "km/s":
				# Working in km/s
				dkms = abs(header["CDELT3"])
				scalemom12 = 1.0
				bunitExt = ".km/s"
			else:
				# Working with whatever units the cube has
				dkms = abs(header["CDELT3"])
				scalemom12 = 1.0
				bunitExt = "." + header["CUNIT3"]
		
		# Frequency
		elif glob.check_values(glob.KEYWORDS_FREQ, header["CTYPE3"]):
			if not "CUNIT3" in header or header["CUNIT3"].lower() == "hz":
				# Using (or assuming) Hz
				dkms = abs(header["CDELT3"])
				scalemom12 = 1.0
				bunitExt = ".Hz"
			elif header["CUNIT3"].lower() == "khz":
				# Converting kHz to Hz
				dkms = abs(header["CDELT3"]) * 1e+3
				scalemom12 = 1e+3
				bunitExt = ".Hz"
			else:
				# Working with whatever frequency units the cube has
				dkms = abs(header["CDELT3"])
				scalemom12 = 1.0
				bunitExt = "." + header["CUNIT3"]
		
		hdu = pyfits.PrimaryHDU(data=m0*dkms, header=header)
		hdu.header["BUNIT"] += bunitExt
		hdu.header["DATAMIN"] = np.nanmin(m0 * dkms)
		hdu.header["DATAMAX"] = np.nanmax(m0 * dkms)
		hdu.header["ORIGIN"] = version.getVersion(full=True)
		del(hdu.header["CRPIX3"])
		del(hdu.header["CRVAL3"])
		del(hdu.header["CDELT3"])
		del(hdu.header["CTYPE3"])
		hdu.header["CELLSCAL"] = "constant"
		
		if debug:
			hdu.writeto("%s_mom0.debug.fits" % filename, output_verify="warn", clobber=True)
		else:
			name = "%s_mom0.fits" % filename
			if compress: name += ".gz"
			
			# Check for overwrite flag
			if not flagOverwrite and os.path.exists(name):
				err.error("Output file exists: " + str(name) + ".", fatal=False)
			else:
				hdu.writeto(name, output_verify="warn", clobber=True)
	
	# --------------
	# Moment 1 image
	# --------------
	if domom1:
		err.message("Writing moment-1")
		
		# Calculate moment 1
		velArr = ((np.arange(datacube.shape[0]) + 1.0 - header["CRPIX3"]) * header["CDELT3"] + header["CRVAL3"]).reshape((datacube.shape[0], 1, 1))
		with np.errstate(invalid="ignore"):
			m1 = np.divide(np.nansum(velArr * datacube, axis=0), m0)
		
		# Velocity
		if glob.check_values(glob.KEYWORDS_VELO, header["CTYPE3"]):
			if not "CUNIT3" in header:
				m1 /= 1e+3 # Assuming m/s
				bunitExt = "km/s"
			elif header["CUNIT3"].lower() == "km/s":
				bunitExt = "km/s"
			else:
				bunitExt = header["CUNIT3"]
		
		# Frequency
		elif glob.check_values(glob.KEYWORDS_FREQ, header["CTYPE3"]):
			if not "CUNIT3" in header or header["CUNIT3"].lower() == "hz":
				bunitExt = "Hz"
			else:
				bunitExt = header["CUNIT3"]
			dkms = 1.0 # No scaling, avoids crashing
		
		hdu = pyfits.PrimaryHDU(data=m1, header=header)
		hdu.header["BUNIT"] = bunitExt
		hdu.header["DATAMIN"] = np.nanmin(m1)
		hdu.header["DATAMAX"] = np.nanmax(m1)
		hdu.header["ORIGIN"] = version.getVersion(full=True)
		del(hdu.header["CRPIX3"])
		del(hdu.header["CRVAL3"])
		del(hdu.header["CDELT3"])
		del(hdu.header["CTYPE3"])
		hdu.header["CELLSCAL"] = "constant"
		
		if debug:
			hdu.writeto("%s_mom1.debug.fits" % filename, output_verify="warn", clobber=True)
		else:
			name = "%s_mom1.fits" % filename
			if compress: name += ".gz"
			
			# Check for overwrite flag
			if not flagOverwrite and os.path.exists(name):
				err.error("Output file exists: " + str(name) + ".", fatal=False)
			else:
				hdu.writeto(name, output_verify="warn", clobber=True)
