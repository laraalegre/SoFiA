#!/usr/bin/python
# -*- coding: utf-8 -*-

# WCS coordinates
# GIPSY header repaired following this description:
# http://www.astro.rug.nl/software/kapteyn/spectralbackground.html#a-recipe-for-modification-of-nmap-gipsy-fits-data

import math
import numpy as np
import scipy.constants
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Longitude, Latitude
from sofia import error as err
import sys


# -------------------------------------
# ---- FUNCTION TO IMPLEMENT SGN() ----
# -------------------------------------

def math_sgn(value):
	# Note: This does not handle special cases such as +-0, +-inf, etc.,
	#       but that's not really needed here anyway.
	if(value < 0): return -1
	return 1


# --------------------------------------------
# ---- FUNCTION TO FIX GIPSY FILE HEADERS ----
# --------------------------------------------

def fix_gipsy_header(header_orig):
	# GIPSY keys for spectral axis
	key_opt = ["FREQ-OHEL","FREQ-OLSR"]
	key_rad = ["FREQ-RHEL","FREQ-RLSR"]
	header = header_orig.copy()
	naxis = header["NAXIS"]
	
	for i in range(1, naxis + 1):
		ctype = header["CTYPE%d" % i]
		if ctype in key_opt + key_rad:
			axis = i
			# Read reference velocity - from VELR or DRVAL
			try:
				if "VELR" in header:
					vel = header["VELR"]
				elif "DRVAL%d" % axis in header:
					vel = header["VELR"]
					unit = header["DUNIT%d" % axis]
					if unit.lower() == "km/s":
						vel = vel * 1000.0
					elif unit.lower() != "m/s":
						break
			except:
				err.warning("Problem with reference velocity.")
				break
			
			# Convert reference frequency to Hz
			try:
				freq  = header["CRVAL%d" % axis]
				dfreq = header["CDELT%d" % axis]
				unit  = header["CUNIT%d" % axis]
				freqUnits = ["hz", "khz", "mhz", "ghz"]
				j = freqUnits.index(unit.lower())
				freq  *= 10**j
				dfreq *= 10**j
			except:
				err.warning("Problem with reference frequency.")
				break
			
			# Need rest frequency for conversion
			try:
				freq0Names = ["FREQ0", "FREQR", "RESTFRQ"]
				for key in freq0Names:
					try:
						freq0 = header[key]
						#foundFreq0 = 1
					except:
						pass
				header["RESTFRQ"] = freq0
				#foundFreq0
			except:
				err.warning("Rest frequency not found.")
				break
			
			# Calculate reference frequency in the barycentric system
			if ctype in key_opt:
				freqB = freq0 / (1.0 + vel / scipy.constants.c)
			else:
				freqB = freq0 / (1.0 - vel / scipy.constants.c)
			
			# Calculate topocentric velocity
			velT = scipy.constants.c * ((freqB**2 - freq**2) / (freqB**2 + freq**2))
			dfreqB = dfreq * math.sqrt((scipy.constants.c - velT) / (scipy.constants.c + velT))
			header["CTYPE%d" % axis] = "FREQ"
			header["CUNIT%d" % axis] = "Hz"
			header["CRVAL%d" % axis] = freqB
			header["CDELT%d" % axis] = dfreqB
			## GIPSY headers seem to contain the unit "DEGREE" for RA/Dec
			## WCS lib does not like that
			for key in header:
				if "CUNIT" in key and header[key] == "DEGREE":
					header[key] = "deg"
			err.message("Header repaired successfully.")
			
			return header


def add_wcs_coordinates(objects, catParNames, catParFormt, catParUnits, Parameters):
	try:
		hdulist = fits.open(Parameters["import"]["inFile"])
		header = hdulist[0].header
		hdulist.close()
		
		# Fix headers where "per second" is written "/S" instead of "/s"
		# (assuming they mean "per second" and not "per Siemens").
		if "cunit3" in header and "/S" in header["cunit3"]:
			err.warning("Converting '/S' to '/s' in CUNIT3.")
			header["cunit3"] = header["cunit3"].replace("/S","/s")
		
		# Check if there is a Nmap/GIPSY FITS header keyword value present
		gipsyKey = [k for k in ["FREQ-OHEL", "FREQ-OLSR", "FREQ-RHEL", "FREQ-RLSR"] if (k in [header[key] for key in header if ("CTYPE" in key)])]
		if gipsyKey:
			err.message("GIPSY header found. Trying to convert to FITS standard.")
			from astropy.wcs import Wcsprm
			header = fix_gipsy_header(header)
			wcsin = Wcsprm(str(header))
			wcsin.sptr("VOPT-F2W")
			objects = np.concatenate((objects, wcsin.p2s(objects[:, catParNames.index("x"):catParNames.index("x") + 3], 0)["world"]), axis=1)
			catParUnits = tuple(list(catParUnits) + [str(cc).replace(" ", "") for cc in wcsin.cunit])
			catParNames = tuple(list(catParNames) + [(cc.split("--")[0]).lower() for cc in wcsin.ctype])
			catParFormt = tuple(list(catParFormt) + ["%15.7e", "%15.7e", "%15.7e"])

		else:
			# Constrain the RA axis reference value CRVAL_ to be between 0 and 360 deg
			rafound = 0
			for kk in range(header["naxis"]):
				if header["ctype1"][:2] == "RA":
					rafound = 1
					break
			if rafound:
				if header["crval%i" % (kk + 1)] < 0:
					err.warning("Adding 360 deg to RA reference value.")
					header["crval%i" % (kk + 1)] += 360
				elif header["crval%i" % (kk + 1)] > 360:
					err.warning("Subtracting 360 deg from RA reference value.")
					header["crval%i" % (kk + 1)] -= 360

			if header['naxis']==2:
				wcsin = wcs.WCS(header, naxis=[wcs.WCSSUB_CELESTIAL])
				xy = objects[:, catParNames.index("x"):catParNames.index("x") + 2].astype(float)
				objects = np.concatenate((objects, wcsin.wcs_pix2world(xy, 0)), axis=1)
				catParUnits = tuple(list(catParUnits) + [str(cc).replace(" ", "") for cc in wcsin.wcs.cunit])
				catParNames = tuple(list(catParNames) + [(cc.split("--")[0]).lower() for cc in wcsin.wcs.ctype])
				catParFormt = tuple(list(catParFormt) + ["%15.7e", "%15.7e"])
			else:
				wcsin = wcs.WCS(header, naxis=[wcs.WCSSUB_CELESTIAL, wcs.WCSSUB_SPECTRAL])
				xyz = objects[:, catParNames.index("x"):catParNames.index("x") + 3].astype(float)
				if "cellscal" in header and header["cellscal"] == "1/F":
					err.warning(
						"CELLSCAL keyword with value of 1/F found.\n"
						"Will account for varying pixel scale in WCS coordinate calculation.")
					x0, y0 = header["crpix1"] - 1, header["crpix2"] - 1
					# Will calculate the pixscale factor of each channel as:
					# pixscale = ref_frequency / frequency
					if header["ctype3"] == "VELO-HEL":
						pixscale = (1 - header["crval3"] / scipy.constants.c) / (1 - (((xyz[:, 2] + 1) - header["crpix3"]) * header["cdelt3"] + header["crval3"]) / scipy.constants.c)
					else:
						err.warning("Cannot convert 3rd axis coordinates to frequency. Ignoring the effect of CELLSCAL = 1/F.")
						pixscale = 1.0
					xyz[:, 0] = (xyz[:, 0] - x0) * pixscale + x0
					xyz[:, 1] = (xyz[:, 1] - y0) * pixscale + y0
				objects = np.concatenate((objects, wcsin.wcs_pix2world(xyz, 0)), axis=1)
				catParUnits = tuple(list(catParUnits) + [str(cc).replace(" ", "") for cc in wcsin.wcs.cunit])
				catParNames = tuple(list(catParNames) + [(cc.split("--")[0]).lower() for cc in wcsin.wcs.ctype])
				catParFormt = tuple(list(catParFormt) + ["%15.7e", "%15.7e", "%15.7e"])
		err.message("WCS coordinates added to catalogue.")

		# Create IAU-compliant source name:
		# WARNING: This currently assumes a regular, ≥ 2-dim. data cube where the first two axes are longitude and latitude.
		n_src = objects.shape[0]
		n_par = objects.shape[1]
		
		iau_names = np.empty([n_src, 1], dtype=object)
		
		if header["ctype1"][:4] == "RA--":
			# Equatorial coordinates; try to figure out equinox:
			iau_coord = "equ"
			if "equinox" in header:
				if int(header["equinox"]) >= 2000: iau_equinox = "J"
				else: iau_equinox = "B"
			elif "epoch" in header:
				# Assume that EPOCH has been abused to record the equinox:
				if int(header["epoch"]) >= 2000: iau_equinox = "J"
				else: iau_equinox = "B"
			else:
				# Equinox undefined:
				iau_equinox = "X"
		elif header["ctype1"][:4] == "GLON":
			# Galactic coordinates:
			iau_coord = "gal"
			iau_equinox = "G"
		else:
			# Unsupported coordinate system:
			iau_coord = ""
			iau_equinox = ""
		
		for src in range(n_src):
			lon = objects[src][n_par - 2] if header['naxis']==2 else objects[src][n_par - 3]
			lat = objects[src][n_par - 1] if header['naxis']==2 else objects[src][n_par - 2]
			
			if iau_coord == "equ":
				ra = Longitude(lon, unit=u.deg)
				dec = Latitude(lat, unit=u.deg)
				iau_pos = ra.to_string(unit=u.h, decimal=False, sep="", precision=2, alwayssign=False, pad=True, fields=3)
				iau_pos += dec.to_string(unit=u.deg, decimal=False, sep="", precision=1, alwayssign=True, pad=True, fields=3)
			else:
				iau_pos = "{0:08.4f}".format(lon)
				if lat < 0.0: iau_pos += "-"
				else: iau_pos += "+"
				iau_pos += "{0:07.4f}".format(abs(lat))
			
			iau_names[src][0] = "SoFiA " + iau_equinox + iau_pos
		
		objects = np.concatenate((objects, iau_names), axis = 1)
		catParUnits = tuple(list(catParUnits) + ["-"])
		catParNames = tuple(list(catParNames) + ["name"])
		catParFormt = tuple(list(catParFormt) + ["%30s"])
	except:
		err.warning("WCS conversion of parameters failed with the following error:")
		err.warning("  {0:}".format(sys.exc_info()))
		err.warning("  at line {0:} of wcs_coordinates.py".format(sys.exc_info()[2].tb_lineno))

	return (objects, catParNames, catParFormt, catParUnits)
