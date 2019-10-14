#!/usr/bin/python

import numpy as np
import scipy.ndimage as nd
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
	dilateChanMax = Parameters["parameters"]["dilateChanMax"]
	
	# Stops dilating when (flux_new - flux_old) / flux_new < dilateThreshold
	sourceIDs = np.unique(mask)
	# remove first element which should be zero
	if sourceIDs[0] == 0:
		sourceIDs = np.delete(sourceIDs,0)
	
	for i in range(0, len(sourceIDs)):
		obj = objects[i]
		xmin = max(0, obj[list(cathead).index("x_min")] - dilatePixMax)
		xmax = min(cube.shape[2] - 1, obj[list(cathead).index("x_max")] + dilatePixMax)
		ymin = max(0, obj[list(cathead).index("y_min")] - dilatePixMax)
		ymax = min(cube.shape[1] - 1, obj[list(cathead).index("y_max")] + dilatePixMax)
		zmin = max(0, obj[list(cathead).index("z_min")] - dilateChanMax)
		zmax = min(cube.shape[0] - 1, obj[list(cathead).index("z_max")] + dilateChanMax)
		
		[zmin, zmax, ymin, ymax, xmin, xmax] = map(int, [zmin, zmax, ymin, ymax, xmin, xmax])
		
		objcube = cube[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1].copy()
		objmask = mask[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1].copy()
		allmask = mask[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1].copy()
		otherobjs = (allmask > 0) * (allmask != sourceIDs[i])
		
		if (otherobjs).sum():
			# Ensure that objects with different source IDs within dilatePixMax, dilateChanMax are not
			# included in the flux growth calculation
			err.warning("Object {0:d} has possible overlapping objects within {1:d} pix, {2:d} chan.".format(sourceIDs[i], dilatePixMax, dilateChanMax))
			objcube[(allmask > 0) * (allmask != sourceIDs[i])] = 0
		
		fluxes = []
		
		# Loop through Z dilation kernels until the flux converges or the maximum allowed Z dilation is reached
		for dilchan in range(dilateChanMax + 1):
			dd = dilchan * 2 + 1
			dilstruct = np.ones((dd,1,1))
			fluxes.append(objcube[nd.morphology.binary_dilation(objmask==sourceIDs[i], structure=dilstruct)].sum())
			if dilchan > 0 and (fluxes[-1] - fluxes[-2]) / fluxes[-1] < dilateThreshold:
				dilchan -= 1
				break
		# Pick the best Z dilation kernel for current object and update mask
		dd = dilchan * 2 + 1
		dilstruct = np.ones((dd,1,1))
		# Only grow the mask of object sourceIDs[i] even when other objects are present in objmask
		objmask[nd.morphology.binary_dilation(objmask==sourceIDs[i], structure=dilstruct).astype(int) == 1] = sourceIDs[i]

		# Loop through XY dilation kernels until the flux converges or the maximum allowed XY dilation is reached
		for dilpix in range(dilatePixMax + 1):
			dd = dilpix * 2 + 1
			dilstruct = (np.sqrt(((np.indices((dd, dd)) - dilpix)**2).sum(axis=0)) <= dilpix).astype(int)
			dilstruct.resize((1, dilstruct.shape[0], dilstruct.shape[1]))
			fluxes.append(objcube[nd.morphology.binary_dilation(objmask==sourceIDs[i], structure=dilstruct)].sum())
			if dilpix > 0 and (fluxes[-1] - fluxes[-2]) / fluxes[-1] < dilateThreshold:
				dilpix -= 1
				break
		# Pick the best XY dilation kernel for current object and update mask
		dd = dilpix * 2 + 1
		dilstruct = (np.sqrt(((np.indices((dd, dd)) - dilpix)**2).sum(axis=0)) <= dilpix).astype(int)
		dilstruct.resize((1, dilstruct.shape[0], dilstruct.shape[1]))
		# Only grow the mask of object sourceIDs[i] even when other objects are present in objmask
		objmask[nd.morphology.binary_dilation(objmask==sourceIDs[i], structure=dilstruct).astype(int) == 1] = sourceIDs[i]
		
		err.message("Mask of source {0:d} dilated by {2:d} chan and then by {1:d} pix.".format(sourceIDs[i], dilpix, dilchan))
		# Put back in objmask objects != sourceIDs[i] that may have been inside objmask before 
		# dilation or may have been temporarily replaced by the dilated object sourceIDs[i]
		if (otherobjs).sum():
			objmask[otherobjs] = allmask[otherobjs]
		mask[zmin:zmax+1, ymin:ymax+1, xmin:xmax+1] = objmask
		
		# Update n_pix, x_geo and n_chan
		n_pix = objmask[objmask == sourceIDs[i]].sum() / sourceIDs[i]
		ind = np.vstack(np.where(objmask == sourceIDs[i]))
		cgeo = (ind.sum(axis=1)).astype(float) / float(n_pix)
		x_geo, y_geo, z_geo = cgeo[2] + xmin, cgeo[1] + ymin, cgeo[0] + zmin
		zmin, zmax = min(ind[0]), max(ind[0]) + 1
		n_chan = zmax - zmin
		
		# Update n_los
		objmask[objmask != sourceIDs[i]] = 0
		maskSumA0 = objmask.sum(axis=0)
		maskSumA0[maskSumA0 > 1] = 1
		n_los = maskSumA0.sum()
		

		del objcube
		del objmask
		del allmask
		del otherobjs
	
		objects[i,list(cathead).index("x_min")]  = max(0, obj[list(cathead).index("x_min")] - dilpix)
		objects[i,list(cathead).index("x_max")]  = min(cube.shape[2] - 1, obj[list(cathead).index("x_max")] + dilpix)
		objects[i,list(cathead).index("y_min")]  = max(0, obj[list(cathead).index("y_min")] - dilpix)
		objects[i,list(cathead).index("y_max")]  = min(cube.shape[1] - 1, obj[list(cathead).index("y_max")] + dilpix)
		objects[i,list(cathead).index("z_min")]  = max(0, obj[list(cathead).index("z_min")] - dilchan)
		objects[i,list(cathead).index("z_max")]  = min(cube.shape[0] - 1, obj[list(cathead).index("z_max")] + dilchan)
		objects[i,list(cathead).index("n_pix")]  = n_pix
		objects[i,list(cathead).index("n_chan")] = n_chan
		objects[i,list(cathead).index("n_los")]  = n_los
		objects[i,list(cathead).index("x_geo")]  = x_geo
		objects[i,list(cathead).index("y_geo")]  = y_geo
		objects[i,list(cathead).index("z_geo")]  = z_geo
		
		
        
	return mask, objects



# ==========================================
# FUNCTION: Call C++ parameterisation module
# ==========================================

def parametrise(cube, mask, objects, cathead, catformt, catparunits, Parameters, dunits):
	cathead = np.array(cathead)
	objects = np.array(objects)
	initcatalog = cp.PySourceCatalog()
	
	for obj in objects:
		# Check flags
		source_flag = create_source_flags(cube, mask, cathead, obj[cathead == "id"], obj[cathead == "x_min"], obj[cathead == "x_max"], obj[cathead == "y_min"], obj[cathead == "y_max"], obj[cathead == "z_min"], obj[cathead == "z_max"])
		
		newSource = cp.PySource()
		newSource.ID = obj[cathead == "id"]
		newParamsDict = {
			"flag": cp.PyMeasurement("flag", source_flag, 0.0, ""),
			"x": cp.PyMeasurement("x", obj[cathead == "x_geo"], 0.0, ""),
			"y": cp.PyMeasurement("y", obj[cathead == "y_geo"], 0.0, ""),
			"z": cp.PyMeasurement("z", obj[cathead == "z_geo"], 0.0, ""),
			"x_min": cp.PyMeasurement("x_min", obj[cathead == "x_min"], 0.0, ""),
			"x_max": cp.PyMeasurement("x_max", obj[cathead == "x_max"], 0.0, ""),
			"y_min": cp.PyMeasurement("y_min", obj[cathead == "y_min"], 0.0, ""),
			"y_max": cp.PyMeasurement("y_max", obj[cathead == "y_max"], 0.0, ""),
			"z_min": cp.PyMeasurement("z_min", obj[cathead == "z_min"], 0.0, ""),
			"z_max": cp.PyMeasurement("z_max", obj[cathead == "z_max"], 0.0, "")
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
		"flag": "-",
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
			catformt.append("%18.6e")
			catparunits.append(newunits[i])
	
	# Extend the parameter array
	tmpObjects = np.empty((objects.shape[0], len(cathead)))
	tmpObjects[:, :] = np.nan
	tmpObjects[:, 0:objects.shape[1]] = objects
	objects = tmpObjects
	
	for i in d:
		source_dict = d[i].getParameters()
		
		# Check source index
		ID = int(source_dict["id"].getValue())
		IDind = cathead.index("id")
		
		ind = np.where(objects[:,IDind]==ID)[0][0]

		for j in sorted(source_dict):
			if j in replParam:
				objects[ind][cathead.index(origParam[replParam.index(j)])] = source_dict[j].getValue()
			else:
				objects[ind][cathead.index(j)] = source_dict[j].getValue()
	
	objects = np.array(objects)
	cathead = np.array(cathead)
	catparunits = np.array(catparunits)
	catformt = np.array(catformt)
	
	# if mask optimization is enabled, some parameters from the linker have to be updated
	if Parameters["parameters"]["optimiseMask"]:
		for i in range(objects.shape[0]):
			# bounding box coordinates
			coord = []
			for c in ["x_min", "x_max", "y_min", "y_max", "z_min", "z_max"]:
				coord.append(int(objects[i, cathead == c]))
			# cut out object submask
			submask = mask[coord[4]:coord[5] + 1, coord[2]:coord[3] + 1, coord[0]:coord[1] + 1]
			objID = objects[i, cathead == "id"]
			submask[submask!=objID] = 0
			
			# Update n_pix, x_geo and n_chan
			n_pix = submask.sum() / objID
			ind = np.vstack(np.where(submask == objID))
			cgeo = (ind.sum(axis=1)).astype(float) / float(n_pix)
			x_geo, y_geo, z_geo = cgeo[2] + coord[0], cgeo[1] + coord[2], cgeo[0] + coord[4]
			zmin, zmax = min(ind[0]), max(ind[0]) + 1
			n_chan = zmax - zmin
			
			# Update n_los
			submaskSumA0 = submask.sum(axis=0)
			submaskSumA0[submaskSumA0 > 1] = 1
			n_los = submaskSumA0.sum()
			
			objects[i, cathead == "n_pix"]  = n_pix
			objects[i, cathead == "n_chan"] = n_chan
			objects[i, cathead == "n_los"]  = n_los
			objects[i, cathead == "x_geo"]  = x_geo
			objects[i ,cathead == "y_geo"]  = y_geo
			objects[i, cathead == "z_geo"]  = z_geo
		del submask
	
	err.message("Parameterisation complete.")
	
	return cube, mask, objects, cathead, catformt, catparunits




# ===============================
# FUNCTION: Generate source flags
# ===============================

def create_source_flags(cube, mask, cathead, obj_id, x_min, x_max, y_min, y_max, z_min, z_max):
	flag = 0
	[z_min, z_max, y_min, y_max, x_min, x_max] = map(int, [z_min, z_max, y_min, y_max, x_min, x_max])
	
	# Check spatial boundaries
	if x_min == 0 or x_max == cube.shape[2] or y_min == 0 or y_max == cube.shape[1]:
		flag |= 1
	
	# Check spectral boundaries
	if z_min == 0 or z_max == cube.shape[0]:
		flag |= 2
	
	# Check for NaN and other sources
	obj_cube = cube[z_min:z_max+1, y_min:y_max+1, x_min:x_max+1].copy()
	obj_mask = mask[z_min:z_max+1, y_min:y_max+1, x_min:x_max+1].copy()
	dil = nd.morphology.binary_dilation(obj_mask==obj_id, structure=np.ones([3, 3, 3]))
	dil_values = obj_cube[dil]
	dil_masked = obj_mask[dil]
	if np.isnan(dil_values).any():
		flag |= 4
	if(dil_masked[dil_masked > 0] != obj_id).any():
		flag |= 8
	
	return flag



# ==================================================
# FUNCTION: Derive linker parameters from input mask
# ==================================================

def parameters_from_mask(dict_Header, mask):
	# If the user only provided an indexed mask, we need to create the objects, catParNames, catParFormt, catParUnits and dunits arrays
	
	# Set catalogue header
	if "bunit" in dict_Header:
		dunits = dict_Header["bunit"]
	else:
		dunits = "-"
	
	# Set parameter info arrays
	catParNames = ("id", "x_geo", "y_geo", "z_geo", "x", "y", "z", "x_min", "x_max", "y_min", "y_max", "z_min", "z_max", "n_pix", "n_chan", "n_los")
	catParUnits = ("-", "pix", "pix", "chan", "pix", "pix", "chan", "pix", "pix", "pix", "pix", "chan", "chan", "-", "chan", "-")
	catParFormt = ("%10i", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%7i", "%7i", "%7i", "%7i", "%7i", "%7i", "%8i", "%7i", "%7i")
	catParNames = np.array(catParNames)
	catParUnits = np.array(catParUnits)
	catParFormt = np.array(catParFormt)
	
	# Create object array from index mask
	indices = np.unique(mask).astype(int)
	
	# Remove first element if it is zero
	if indices[0] == 0:
		indices = np.delete(indices, 0)
	
	# Make object array
	objects = np.empty((len(indices), len(catParNames)))
	objects[:,:] = np.nan
	
	# Fill with object parameters needed for parameterisation
	for i in range(len(indices)):
				# Add object IDs
				objects[i, catParNames=="id"] = indices[i]
				# Add boundaries and geometric centre
				xmin, xmax, ymin, ymax, zmin, zmax, xgeo, ygeo, zgeo, n_pix, n_los, n_chan = src_params_from_mask(mask, indices[i])
				objects[i, catParNames=="x_min"] = xmin
				objects[i, catParNames=="x_max"] = xmax
				objects[i, catParNames=="y_min"] = ymin
				objects[i, catParNames=="y_max"] = ymax
				objects[i, catParNames=="z_min"] = zmin
				objects[i, catParNames=="z_max"] = zmax
				objects[i, catParNames=="x_geo"] = xgeo
				objects[i, catParNames=="y_geo"] = ygeo
				objects[i, catParNames=="z_geo"] = zgeo
				objects[i, catParNames=="n_pix"] = n_pix
				objects[i, catParNames=="n_chan"] = n_chan
				objects[i, catParNames=="n_los"] = n_los
	
	return catParNames, catParUnits, catParFormt, objects, dunits



# ===================================================
# FUNCTION: Extract basic source parameters from mask
# ===================================================

def src_params_from_mask(mask, ID):
	# Get coordinates of pixels, that contain the selected ID
	indices = np.vstack(np.where(mask==ID))
	
	# Min and max (x,y,z) values represent the bounding boxes
	xmin, xmax = min(indices[2]), max(indices[2]) + 1
	ymin, ymax = min(indices[1]), max(indices[1]) + 1
	zmin, zmax = min(indices[0]), max(indices[0]) + 1
	
	# Total number of pixes in the source
	n_pix = indices.shape[1]
	
	# Geometric center
	cgeo = (indices.sum(axis=1)).astype(float) / float(n_pix)
	xgeo, ygeo, zgeo = cgeo[2], cgeo[1], cgeo[0]
	
	# Number of lines of sight
	collapsedMap = mask[zmin:zmax, ymin:ymax, xmin:xmax].sum(axis=0)
	collapsedMap[collapsedMap > 0] = 1
	n_los = collapsedMap.sum()
	
	# Number of channels
	n_chan = zmax - zmin
	
	return xmin, xmax, ymin, ymax, zmin, zmax, xgeo, ygeo, zgeo, n_pix, n_los, n_chan
