#!/usr/bin/python

import sys
import numpy as np
import scipy.ndimage as nd
from scipy.optimize import curve_fit
from scipy.special import erf
from sofia import cparametrizer as cp
from sofia import error as err


# ==============================
# FUNCTION: Define Busy Function
# ==============================

def BusyFunction(x, a, b1, b2, c, w, xe, xp):
	return (a / 4.0) * (erf(b1 * (w + x - xe)) + 1) * (erf(b2 * (w - x + xe)) + 1) * (c * np.absolute(x - xp) * np.absolute(x - xp) + 1)



# =======================
# FUNCTION: Mask dilation
# =======================

def dilate(cube, mask, objects, cathead, Parameters):
	dilateThreshold = Parameters["parameters"]["dilateThreshold"]
	dilatePixMax = Parameters["parameters"]["dilatePixMax"]
	dilateChan = Parameters["parameters"]["dilateChan"]
	
	# Stops dilating when (flux_new - flux_old) / flux_new < dilateThreshold
	for mm in range(1, mask.max() + 1):
		obj = objects[mm - 1]
		xmin = max(0, obj[list(cathead).index("x_min")] - dilatePixMax)
		xmax = min(cube.shape[2] - 1, obj[list(cathead).index("x_max")] + dilatePixMax)
		ymin = max(0, obj[list(cathead).index("y_min")] - dilatePixMax)
		ymax = min(cube.shape[1] - 1, obj[list(cathead).index("y_max")] + dilatePixMax)
		zmin = max(0, obj[list(cathead).index("z_min")] - dilateChan)
		zmax = min(cube.shape[0] - 1, obj[list(cathead).index("z_max")] + dilateChan)
		
		[zmin, zmax, ymin, ymax, xmin, xmax] = map(int, [zmin, zmax, ymin, ymax, xmin, xmax])
		
		objcube = cube[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1].copy()
		objmask = mask[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1].copy()
		allmask = mask[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1].copy()
		otherobjs = (allmask > 0) * (allmask != mm)
		
		if (otherobjs).sum():
			# Ensure that objects != mm within dilatePixMax, dilateChan are not
			# included in the flux growth calculation
			err.warning("Object {0:d} has possible overlapping objects within {1:d} pix, {2:d} chan.".format(mm, dilatePixMax, dilateChan))
			objcube[(allmask > 0) * (allmask != mm)] = 0
		
		fluxes = []
		
		for dil in range(dilatePixMax + 1):
			dd = dil * 2 + 1
			dilstruct = (np.sqrt(((np.indices((dd, dd)) - dil)**2).sum(axis=0)) <= dil).astype(int)
			dilstruct.resize((1, dilstruct.shape[0], dilstruct.shape[1]))
			dilstruct = dilstruct.repeat(dilateChan * 2 + 1, axis=0)
			fluxes.append(objcube[nd.morphology.binary_dilation(objmask==mm, structure=dilstruct)].sum())
			if dil > 0 and (fluxes[-1] - fluxes[-2]) / fluxes[-1] < dilateThreshold:
				dil -= 1
				break
		
		# Pick the best dilation kernel for current object and update mask
		err.message("Mask dilation of source {0:d} by {1:d} px and {2:d} chan.".format(mm, dil, dilateChan))
		dd = dil * 2 + 1
		dilstruct = (np.sqrt(((np.indices((dd, dd)) - dil)**2).sum(axis=0)) <= dil).astype(int)
		dilstruct.resize((1, dilstruct.shape[0], dilstruct.shape[1]))
		dilstruct = dilstruct.repeat(dilateChan * 2 + 1, axis=0)
		
		# Only grow the mask of object mm even when other objects are present in objmask
		objmask[nd.morphology.binary_dilation(objmask==mm, structure=dilstruct).astype(int) == 1] = mm
		
		# Put back in objmask objects != mm that may have been inside objmask before 
		# dilation or may have been temporarily replaced by the dilated object mm
		if (otherobjs).sum():
			objmask[otherobjs] = allmask[otherobjs]
		mask[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1] = objmask
		
		del objcube
		del objmask
		del allmask
		del otherobjs
	return mask



# ==========================================
# FUNCTION: Call C++ parameterisation module
# ==========================================

def parametrise(cube, mask, objects, cathead, catformt, catparunits, Parameters, dunits):
	cathead = np.array(cathead)
	objects = np.array(objects)
	initcatalog = cp.PySourceCatalog()
	
	for obj in objects:
		newSource = cp.PySource()
		newSource.ID = obj[cathead == "id"]
		newParamsDict = {
			"x": cp.PyMeasurement("x", obj[cathead == "x_geo"], 0.0, ""),
			"y": cp.PyMeasurement("y", obj[cathead == "y_geo"], 0.0, ""),
			"z": cp.PyMeasurement("z", obj[cathead == "z_geo"], 0.0, ""),
			"x_min": cp.PyMeasurement("x_min", obj[cathead == "x_min"], 0.0, ""),
			"x_max": cp.PyMeasurement("x_max", obj[cathead == "x_max"], 0.0, ""),
			"y_min": cp.PyMeasurement("y_min", obj[cathead == "y_min"], 0.0, ""),
			"y_max": cp.PyMeasurement("y_max", obj[cathead == "y_max"], 0.0, ""),
			"z_min": cp.PyMeasurement("z_min", obj[cathead == "z_min"], 0.0, ""),
			"z_max": cp.PyMeasurement("z_max", obj[cathead == "z_max"], 0.0, ""),
			}
		newSource.setParameters(newParamsDict)
		initcatalog.insert(newSource)
	
	moduleParametrizer = cp.PyModuleParametrisation()
	moduleParametrizer.setFlags(Parameters["parameters"]["optimiseMask"], Parameters["parameters"]["fitBusyFunction"])
	
	cube = cube.astype("<f4", copy=False)
	mask = mask.astype("<i2", copy=False)
	
	moduleParametrizer.run(cube, mask, initcatalog)
	results = moduleParametrizer.getCatalog()
	
	# Append the results to the objects array or reset
	replParam = ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max", "id", "x", "y", "z", "n_pix"]
	origParam = ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max", "id", "x", "y", "z", "n_pix"]
	d = results.getSources()
	
	# Select data set with maximum number of parameters
	parsListLen = [len(d[list(d.keys())[i]].getParameters()) for i in range(0, len(d))]
	index = parsListLen.index(max(parsListLen))
	
	# Add parameter names from parameterisation
	pars = d[list(d.keys())[index]].getParameters()
	cathead = list(cathead)
	newunits = {
		"id": "-",
		"x": "pix",
		"y": "pix",
		"z": "pix",
		"err_x": "pix",
		"err_y": "pix",
		"err_z": "pix",
		"x_min": "pix",
		"x_max": "pix",
		"y_min": "pix",
		"y_max": "pix",
		"z_min": "chan",
		"z_max": "chan",
		"w50": "chan",
		"w20": "chan",
		"err_w50": "chan",
		"err_w20": "chan",
		"wm50": "chan",
		"f_wm50": dunits,
		"ell_maj": "pix",
		"ell_min": "pix",
		"ell_pa": "deg",
		"ell3s_maj": "pix",
		"ell3s_min": "pix",
		"ell3s_pa": "deg",
		"kin_pa": "deg",
		"f_int": dunits,
		"bf_flag": "-",
		"bf_chi2": "-",
		"bf_z": "chan",
		"bf_a": dunits,
		"bf_b1": "chan**(-1)",
		"bf_b2": "chan**(-1)",
		"bf_c": "chan**(-2)",
		"bf_xe": "chan",
		"bf_xp": "chan",
		"bf_w": "chan",
		"bf_w50": "chan",
		"bf_w20": "chan",
		"bf_f_peak": dunits,
		"bf_f_int": dunits,
		"rms": dunits,
		"f_peak": dunits,
		"snr_int": "-"
		}
	catformt = list(catformt)
	catparunits = list(catparunits)
	
	for i in sorted(pars):
		if i not in replParam:
			cathead.append(i)
			catformt.append("%12.4f")
			catparunits.append(newunits[i])
	
	# Extend the parameter array
	tmpObjects = np.empty((objects.shape[0], len(cathead)))
	tmpObjects[:, :] = np.nan
	tmpObjects[:, 0:objects.shape[1]] = objects
	objects = tmpObjects
	
	for i in d:
		source_dict = d[i].getParameters()
		
		# Check source index
		index = int(source_dict["id"].getValue())
		for j in sorted(source_dict):
			if j in replParam:
				objects[index - 1][cathead.index(origParam[replParam.index(j)])] = source_dict[j].getValue()
			else:
				objects[index - 1][cathead.index(j)] = source_dict[j].getValue()
	
	objects = np.array(objects)
	cathead = np.array(cathead)
	catparunits = np.array(catparunits)
	catformt = np.array(catformt)
	
	err.message("Parameterisation complete.")
	
	return cube, mask, objects, cathead, catformt, catparunits
