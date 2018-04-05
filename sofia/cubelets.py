#! /usr/bin/env python
import numpy as np
from astropy.io import fits
import os
import sys
from sofia import writemoment
import math
from scipy import interpolate
from scipy.ndimage import map_coordinates
import scipy.constants
from .version import *
from sofia import error as err


# -----------------------------------------------------
# FUNCTION: Regrid data cube for moment map calculation
# -----------------------------------------------------

def regridMaskedChannels(datacube, maskcube, header):
	maskcubeFlt = maskcube.astype("float")
	maskcubeFlt[maskcubeFlt > 1] = 1
	
	err.message("Regridding...")
	z = (np.arange(1., header["naxis3"] + 1) - header["CRPIX3"]) * header["CDELT3"] + header["CRVAL3"]
	if "vopt" in header["CTYPE3"].lower() or "vrad" in header["CTYPE3"].lower() or "velo" in header["CTYPE3"].lower() or "felo" in header["CTYPE3"].lower():
		pixscale=(1 - header["CRVAL3"] / scipy.constants.c) / (1 - z / scipy.constants.c)
		# WARNING: Strictly correct only for the radio velocity definition!
	elif "freq" in header["CTYPE3"].lower():
		pixscale = header["CRVAL3"] / z
	else:
		err.warning("Cannot convert 3rd axis coordinates to frequency. Will ignore the effect of CELLSCAL = 1/F.")
		pixscale = np.ones((header["naxis3"]))
	
	x0, y0 = header["crpix1"] - 1, header["crpix2"] - 1
	xs = np.arange(datacube.shape[2], dtype=float) - x0
	ys = np.arange(datacube.shape[1], dtype=float) - y0
	
	for zz in range(datacube.shape[0]):
		regrid_channel = interpolate.RectBivariateSpline(ys * pixscale[zz], xs * pixscale[zz], datacube[zz])
		datacube[zz] = regrid_channel(ys, xs)
		regrid_channel_mask = interpolate.RectBivariateSpline(ys * pixscale[zz], xs * pixscale[zz], maskcubeFlt[zz])
		maskcubeFlt[zz] = regrid_channel_mask(ys, xs)
	
	datacube[abs(maskcubeFlt) <= abs(maskcubeFlt.min())] = np.nan
	del maskcubeFlt
	return datacube


# ------------------------------------------------------
# FUNCTION: Create various data products for each source
# ------------------------------------------------------

def writeSubcube(cube, header, mask, objects, cathead, outroot, outputDir, compress, flagOverwrite):
	# Strip path variable to get the file name and the directory separately
	splitroot = outroot.split("/")
	cubename  = splitroot[-1]
	#if len(splitroot) > 1:
	#	outputDir = "/".join(splitroot[:-1]) + "/objects/"
	#else:
	#	outputDir = "./objects/"
	
	# Check if output directory exists and create it if not
	if not os.path.exists(outputDir):
		os.system("mkdir " + outputDir)
	
	# Copy of header for manipulation
	headerCubelets = header.copy()
	
	# Read all important information (central pixels & values, increments) from the header
	dX    = headerCubelets["CDELT1"]
	dY    = headerCubelets["CDELT2"]
	dZ    = headerCubelets["CDELT3"]
	cValX = headerCubelets["CRVAL1"]
	cValY = headerCubelets["CRVAL2"]
	cValZ = headerCubelets["CRVAL3"]
	cPixX = headerCubelets["CRPIX1"] - 1
	cPixY = headerCubelets["CRPIX2"] - 1
	cPixZ = headerCubelets["CRPIX3"] - 1
	cubeDim = cube.shape
	
	for obj in objects:
		# Centres and bounding boxes
		Xc = obj[cathead == "x"][0]
		Yc = obj[cathead == "y"][0]
		Zc = obj[cathead == "z"][0]
		Xmin = obj[cathead == "x_min"][0]
		Ymin = obj[cathead == "y_min"][0]
		Zmin = obj[cathead == "z_min"][0]
		Xmax = obj[cathead == "x_max"][0]
		Ymax = obj[cathead == "y_max"][0]
		Zmax = obj[cathead == "z_max"][0]
		
		# If centre of mass estimation is wrong replace by geometric centre
		if Xc < 0 or Xc > cubeDim[2] - 1: Xc = obj[cathead == "x_geo"][0]
		if Yc < 0 or Yc > cubeDim[1] - 1: Yc = obj[cathead == "y_geo"][0]
		if Zc < 0 or Zc > cubeDim[0] - 1: Zc = obj[cathead == "z_geo"][0]
		
		cPixXNew = int(Xc)
		cPixYNew = int(Yc)
		cPixZNew = int(Zc)
		
		# Largest distance of source limits from the centre
		maxX = 2 * max(abs(cPixXNew - Xmin), abs(cPixXNew - Xmax))
		maxY = 2 * max(abs(cPixYNew - Ymin), abs(cPixYNew - Ymax))
		maxZ = 2 * max(abs(cPixZNew - Zmin), abs(cPixZNew - Zmax))
		
		# Calculate the new bounding box for the mass centred cube
		XminNew = cPixXNew - maxX
		if XminNew < 0: XminNew = 0
		YminNew = cPixYNew - maxY
		if YminNew < 0: YminNew = 0
		ZminNew = cPixZNew - maxZ
		if ZminNew < 0: ZminNew = 0
		XmaxNew = cPixXNew + maxX
		if XmaxNew > cubeDim[2] - 1: XmaxNew = cubeDim[2] - 1
		YmaxNew = cPixYNew + maxY
		if YmaxNew > cubeDim[1] - 1: YmaxNew = cubeDim[1] - 1
		ZmaxNew = cPixZNew + maxZ
		if ZmaxNew > cubeDim[0] - 1: ZmaxNew = cubeDim[0] - 1
		
		# Calculate the centre with respect to the cutout cube
		cPixXCut = cPixX - XminNew
		cPixYCut = cPixY - YminNew
		cPixZCut = cPixZ - ZminNew
		
		# Update header keywords:
		headerCubelets["CRPIX1"] = cPixXCut + 1
		headerCubelets["CRPIX2"] = cPixYCut + 1
		headerCubelets["CRPIX3"] = cPixZCut + 1
		
		# Extract the cubelet
		[ZminNew, ZmaxNew, YminNew, YmaxNew, XminNew, XmaxNew] = map(int, [ZminNew, ZmaxNew, YminNew, YmaxNew, XminNew, XmaxNew])
		subcube = cube[ZminNew:ZmaxNew + 1, YminNew:YmaxNew + 1, XminNew:XmaxNew + 1]
		
		# Update header keywords:
		headerCubelets["NAXIS1"] = subcube.shape[2]
		headerCubelets["NAXIS2"] = subcube.shape[1]
		headerCubelets["NAXIS3"] = subcube.shape[0]
		
		headerCubelets["ORIGIN"] = getVersion(full=True)
		
		# Write the cubelet
		hdu = fits.PrimaryHDU(data=subcube, header=headerCubelets)
		hdulist = fits.HDUList([hdu])
		name = outputDir + cubename + "_" + str(int(obj[0])) + ".fits"
		if compress: name += ".gz"
		
		# Check for overwrite flag:
		if check_overwrite(name, flagOverwrite): hdulist.writeto(name, output_verify="warn", clobber=True)
		
		hdulist.close()
		
		
		# -------------------------
		# Position-velocity diagram
		# -------------------------
		
		if "kin_pa" in cathead:
			pv_sampling = 10
			pv_r = np.arange(-max(subcube.shape[1:]), max(subcube.shape[1:]) - 1 + 1.0 / pv_sampling, 1.0 / pv_sampling)
			pv_y = Yc - float(YminNew) + pv_r * math.cos(float(obj[cathead == "kin_pa"][0]) / 180 * math.pi)
			pv_x = Xc - float(XminNew) - pv_r * math.sin(float(obj[cathead == "kin_pa"][0]) / 180 * math.pi)
			pv_x, pv_y = pv_x[(pv_x >= 0) * (pv_x <= subcube.shape[2] - 1)], pv_y[(pv_x >= 0) * (pv_x <= subcube.shape[2] - 1)]
			pv_x, pv_y = pv_x[(pv_y >= 0) * (pv_y <= subcube.shape[1] - 1)], pv_y[(pv_y >= 0) * (pv_y <= subcube.shape[1] - 1)]
			pv_x.resize((1, pv_x.shape[0]))
			pv_y.resize((pv_x.shape))
			pv_coords = np.concatenate((pv_y, pv_x), axis=0)
			pv_array=[]
			for jj in range(subcube.shape[0]):
				plane = map_coordinates(subcube[jj], pv_coords)
				plane = [plane[ii::pv_sampling] for ii in range(pv_sampling)]
				plane = np.array([ii[:plane[-1].shape[0]] for ii in plane])
				pv_array.append(plane.mean(axis=0))
			pv_array = np.array(pv_array)
			hdu = fits.PrimaryHDU(data=pv_array, header=headerCubelets)
			hdulist = fits.HDUList([hdu])
			hdulist[0].header["CTYPE1"] = "PV--DIST"
			hdulist[0].header["CDELT1"] = hdulist[0].header["CDELT2"]
			hdulist[0].header["CRVAL1"] = 0
			hdulist[0].header["CRPIX1"] = pv_array.shape[1] / 2
			hdulist[0].header["CTYPE2"] = hdulist[0].header["CTYPE3"]
			hdulist[0].header["CDELT2"] = hdulist[0].header["CDELT3"]
			hdulist[0].header["CRVAL2"] = hdulist[0].header["CRVAL3"]
			hdulist[0].header["CRPIX2"] = hdulist[0].header["CRPIX3"]
			hdulist[0].header["ORIGIN"] = getVersion(full=True)
			delete_3rd_axis(hdulist[0].header)
			name = outputDir + cubename + "_" + str(int(obj[0])) + "_pv.fits"
			if compress: name += ".gz"
			
			# Check for overwrite flag:
			if check_overwrite(name, flagOverwrite): hdulist.writeto(name,output_verify="warn", clobber=True)
			hdulist.close()
		
		# Remove all other sources from the mask
		submask = mask[ZminNew:ZmaxNew + 1, YminNew:YmaxNew + 1, XminNew:XmaxNew + 1].astype("int")
		submask[submask != obj[0]] = 0
		submask[submask == obj[0]] = 1
		
		# Write mask
		hdu = fits.PrimaryHDU(data=submask.astype("int16"), header=headerCubelets)
		hdu.header["BUNIT"] = "Source-ID"
		hdu.header["DATAMIN"] = submask.min()
		hdu.header["DATAMAX"] = submask.max()
		hdu.header["ORIGIN"] = getVersion(full=True)
		hdulist = fits.HDUList([hdu])
		name = outputDir + cubename + "_" + str(int(obj[0])) + "_mask.fits"
		if compress: name += ".gz"
		
		# Check for overwrite flag:
		if check_overwrite(name, flagOverwrite): hdulist.writeto(name, output_verify="warn", clobber=True)
		hdulist.close()
		
		# Units of moment images
		# Velocity
		if "vopt" in headerCubelets["CTYPE3"].lower() or "vrad" in headerCubelets["CTYPE3"].lower() or "velo" in headerCubelets["CTYPE3"].lower() or "felo" in headerCubelets["CTYPE3"].lower():
			if not "CUNIT3" in headerCubelets or headerCubelets["CUNIT3"].lower() == "m/s":
				# Converting m/s to km/s
				dkms = abs(headerCubelets["CDELT3"]) / 1e+3
				scalemom12 = 1.0 / 1e+3
				bunitExt = ".km/s"
			elif headerCubelets["CUNIT3"].lower() == "km/s":
				dkms = abs(headerCubelets["CDELT3"])
				scalemom12 = 1.0
				bunitExt = ".km/s"
			else:
				# Working with whatever units the cube has
				dkms = abs(headerCubelets["CDELT3"])
				scalemom12 = 1.0
				bunitExt = "." + headerCubelets["CUNIT3"]
		# Frequency
		elif "freq" in headerCubelets["CTYPE3"].lower():
			if not "CUNIT3" in headerCubelets or headerCubelets["CUNIT3"].lower() == "hz":
				dkms = abs(headerCubelets["CDELT3"])
				scalemom12 = 1.0
				bunitExt = ".Hz"
			elif headerCubelets["CUNIT3"].lower() == "khz":
				# Converting kHz to Hz
				dkms = abs(headerCubelets["CDELT3"]) * 1e+3
				scalemom12 = 1e+3
				bunitExt = ".Hz"
			else:
				# Working with whatever units the cube has
				dkms = abs(headerCubelets["CDELT3"])
				scalemom12 = 1.0
				bunitExt = "." + headerCubelets["CUNIT3"]
		# Other
		else:
			# Working with whatever units the cube has
			dkms = abs(headerCubelets["CDELT3"])
			scalemom12 = 1.0
			if not "CUNIT3" in headerCubelets: bunitExt = ".std_unit_" + headerCubelets["CTYPE3"]
			else: bunitExt = "." + headerCubelets["CUNIT3"]
		
		# Make copy of subcube and regrid
		subcubeCopy = subcube.copy()
		subcubeCopy[submask == 0] = 0
		if "cellscal" in headerCubelets:
			if headerCubelets["cellscal"] == "1/F": subcubeCopy = regridMaskedChannels(subcubeCopy, submask, headerCubelets)
		
		
		# ------------------
		# Moments 0, 1 and 2
		# ------------------
		
		moments = [None, None, None]
		with np.errstate(invalid="ignore"):
			moments[0] = np.nansum(subcubeCopy, axis=0)
			tmp = ((np.arange(subcubeCopy.shape[0]).reshape((subcubeCopy.shape[0], 1, 1)) * np.ones(subcubeCopy.shape) - headerCubelets["CRPIX3"] + 1) * headerCubelets["CDELT3"] + headerCubelets["CRVAL3"]) * scalemom12
			moments[1] = np.divide(np.nansum(tmp * subcubeCopy, axis=0), moments[0])
			tmp -= moments[1]
			moments[2] = np.sqrt(np.divide(np.nansum(tmp * tmp * subcubeCopy, axis=0), moments[0]))
		
		moments[0] *= dkms
		units = [headerCubelets["BUNIT"] + bunitExt, bunitExt[1:], bunitExt[1:]]
		
		for i in range(3):
			hdu = fits.PrimaryHDU(data=moments[i], header=headerCubelets)
			delete_3rd_axis(hdu.header)
			hdu.header["BUNIT"]   = units[i]
			hdu.header["DATAMIN"] = np.nanmin(moments[i])
			hdu.header["DATAMAX"] = np.nanmax(moments[i])
			hdu.header["ORIGIN"]  = getVersion(full=True)
			filename = outputDir + cubename + "_" + str(int(obj[0])) + "_mom" + str(i) + ".fits"
			if compress: filename += ".gz"
			if check_overwrite(filename, flagOverwrite): hdu.writeto(filename, output_verify="warn", clobber=True)
		
		
		# -------------------
		# Integrated spectrum
		# -------------------
		spec = np.nansum(subcubeCopy, axis=(1, 2))
		nPix = np.sum(~np.isnan(subcubeCopy), axis=(1, 2))
		
		name = outputDir + cubename + "_" + str(int(obj[0])) + "_spec.txt"
		if compress: name += ".gz"
		
		# Check for overwrite flag:
		if check_overwrite(name, flagOverwrite):
			if compress:
				import gzip
				f = gzip.open(name, "wb")
			else:
				f = open(name, "w")
			
			f.write("# Integrated source spectrum\n")
			f.write("# Creator: %s\n#\n" % getVersion(full=True))
			f.write("# Description of columns:\n")
			f.write("# - Chan      Channel number.\n")
			f.write("# - Spectral  Associated value of the spectral coordinate according to\n")
			f.write("#             the WCS information in the FITS file header.\n")
			f.write("# - Sum       Sum of flux values of all spatial pixels covered by the\n")
			f.write("#             source in that channel. Note that this has not yet been\n")
			f.write("#             divided by the beam solid angle! If your data cube is in\n")
			f.write("#             Jy/beam, you will have to manually divide by the beam\n")
			f.write("#             size which, for Gaussian beams, is given as\n")
			f.write("#               PI * a * b / (4 * ln(2))\n")
			f.write("#             where a and b are the major and minor axis of the beam in\n")
			f.write("#             units of pixels.\n")
			f.write("# - Npix      Number of spatial pixels covered by the source in that\n")
			f.write("#             channel. This can be used to determine the statistical\n")
			f.write("#             uncertainty of the summed flux value. Again, this has\n")
			f.write("#             not yet been corrected for any potential spatial correla-\n")
			f.write("#             tion of pixels due to the beam solid angle!\n#\n")
			f.write("# Chan        Spectral             Sum    Npix\n")
			f.write("# --------------------------------------------\n")
			
			for i in range(0,len(spec)):
				xspec = cValZ + (i + float(ZminNew) - cPixZ) * dZ
				f.write("%6d %15.6e %15.6e %7d\n" % (i + ZminNew, xspec, spec[i], nPix[i]))
			
			f.close()


# -------------------------------------
# FUNCTION: Check file overwrite status
# -------------------------------------

def check_overwrite(filename, flagOverwrite, fatal=False):
	if not flagOverwrite and os.path.exists(filename):
		err.error("Output file exists: " + filename + ".", fatal=fatal)
		return False
	return True


# -----------------------------------------
# FUNCTION: Delete header element if exists
# -----------------------------------------

def delete_header(header, element):
	if element in header:
		del(header[element])
		return True
	return False


# ---------------------------------------------
# FUNCTION: Delete 3rd-axis WCS header elements
# ---------------------------------------------

def delete_3rd_axis(header):
	delete_header(header, "CRPIX3")
	delete_header(header, "CRVAL3")
	delete_header(header, "CDELT3")
	delete_header(header, "CTYPE3")
	delete_header(header, "CUNIT3")
	return
