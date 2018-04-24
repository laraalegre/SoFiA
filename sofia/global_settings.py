#! /usr/bin/env python

# The purpose of this file is to define global constants and functions
# that are needed across modules.



# =========================================
# SETTINGS: FITS header keywords and values
# =========================================

KEYWORDS_VELO = ["VOPT", "VRAD", "VELO", "FELO"]  # NOTE: "FELO" is not a valid FITS coordinate type!
KEYWORDS_FREQ = ["FREQ"]



# ===============================================================
# FUNCTION: Check if keyword contains any from the list of values
# ===============================================================

def check_header_keywords(values, keyword):
	for value in values:
		if value in keyword.upper(): return True
	return False



# =====================================================
# FUNCTION: Delete header keywords with existence check
# =====================================================

def delete_header(header, keyword):
	if keyword in header:
		del header[keyword]
		return True
	return False



# =============================================
# FUNCTION: Delete 3rd-axis WCS header elements
# =============================================

def delete_3rd_axis(header):
	delete_header(header, "CTYPE3")
	delete_header(header, "CRPIX3")
	delete_header(header, "CRVAL3")
	delete_header(header, "CDELT3")
	delete_header(header, "CUNIT3")
	return



# ============================================================
# FUNCTION: Ensure that basic WCS axis descriptors are present
# ============================================================

def check_wcs_info(header, spatial=False):
	if spatial: return "CTYPE3" in header and "CDELT3" in header and "CRPIX3" in header and "CRVAL3" in header \
		and "CTYPE2" in header and "CDELT2" in header and "CRPIX2" in header and "CRVAL2" in header \
		and "CTYPE1" in header and "CDELT1" in header and "CRPIX1" in header and "CRVAL1" in header
	return "CTYPE3" in header and "CDELT3" in header and "CRPIX3" in header and "CRVAL3" in header



# ==========================
# FUNCTION: Regrid data cube
# ==========================

def regridMaskedChannels(datacube, maskcube, header):
	import numpy as np
	import scipy.constants
	from scipy import interpolate
	from sofia import error as err
	
	if not check_wcs_info(header, spatial=True):
		err.warning("Axis descriptors missing from FITS file header.\nIgnoring the effect of CELLSCAL = 1/F.")
		return datacube
	
	maskcubeFlt = maskcube.astype("float")
	maskcubeFlt[maskcube > 1] = 1.0
	
	err.message("Regridding...")
	z = (np.arange(1.0, header["NAXIS3"] + 1) - header["CRPIX3"]) * header["CDELT3"] + header["CRVAL3"]
	
	if check_header_keywords(KEYWORDS_VELO, header["CTYPE3"]):
		pixscale = (1.0 - header["CRVAL3"] / scipy.constants.c) / (1.0 - z / scipy.constants.c)
	elif check_header_keywords(KEYWORDS_FREQ, header["CTYPE3"]):
		pixscale = header["CRVAL3"] / z
	else:
		err.warning("Cannot convert 3rd axis coordinates to frequency.\nIgnoring the effect of CELLSCAL = 1/F.")
		pixscale = np.ones((header["NAXIS3"]))
	
	x0 = header["CRPIX1"] - 1
	y0 = header["CRPIX2"] - 1
	xs = np.arange(datacube.shape[2], dtype=float) - x0
	ys = np.arange(datacube.shape[1], dtype=float) - y0
	
	for zz in range(datacube.shape[0]):
		regrid_channel = interpolate.RectBivariateSpline(ys * pixscale[zz], xs * pixscale[zz], datacube[zz])
		datacube[zz] = regrid_channel(ys, xs)
		regrid_channel_mask = interpolate.RectBivariateSpline(ys * pixscale[zz], xs * pixscale[zz], maskcubeFlt[zz])
		maskcubeFlt[zz] = regrid_channel_mask(ys, xs)
	
	datacube[abs(maskcubeFlt) <= abs(np.nanmin(maskcubeFlt))] = np.nan
	del maskcubeFlt
	
	return datacube
