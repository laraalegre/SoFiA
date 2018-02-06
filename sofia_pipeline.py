#! /usr/bin/env python

# import default python libraries
import numpy as np
import sys, os
import string
from time import time

from scipy import __version__ as scipy_version
from astropy import __version__ as astropy_version

# import source finding modules
sys.path.insert(0, os.environ['SOFIA_MODULE_PATH'])
from sofia import functions
from sofia import readoptions
from sofia import import_data
from sofia import sigma_cube
from sofia import pyfind
from sofia import wavelet_finder
from sofia import addrel
from sofia import threshold_filter
from sofia import smooth_cube
from sofia import write_filtered_cube
from sofia import writemask
from sofia import writemoment2
from sofia import write_catalog
from sofia import linker
from sofia import cubelets
from sofia import parametrisation
from sofia import wcs_coordinates
from sofia import flag_cube
from sofia import CNHI
from sofia import version



# --------------------------------------
# ---- FUNCTION TO CHECK OVERWRITES ----
# --------------------------------------

def checkOverwrite(outputFile, isFile=True, isDir=False):
	if isFile and os.path.exists(outputFile):
		sys.stderr.write("ERROR: SoFiA tried to create the output file %s, but that file already exists.\n" % outputFile)
		sys.stderr.write("       You can do one of the following:\n")
		sys.stderr.write("       1) enable overwrite in GUI or parameter file\n")
		sys.stderr.write("       2) rename existing file\n")
		sys.stderr.write("       3) change base name and/or output directory in GUI or parameter file\n")
		sys.exit(1)
	elif isDir and os.path.exists(outputFile) and os.listdir(outputFile) != []:
		sys.stderr.write("ERROR: SoFiA tried to create the output directory %s, but that directory already exists and is not empty.\n" % outputFile)
		sys.stderr.write("       You can do one of the following:\n")
		sys.stderr.write("       1) enable overwrite in GUI or parameter file\n")
		sys.stderr.write("       2) rename existing directory\n")
		sys.stderr.write("       3) change base name and/or output directory in GUI or parameter file\n")
		sys.exit(1)
	return



# --------------------------------------------
# ---- FUNCTION TO PRINT PROGRESS MESSAGE ----
# --------------------------------------------

def printProgressMessage(message):
	print ("\n--- %.3f seconds since start" % (time() - t0))
	print ("--- %s: %s" % (version.getVersion(full=True), message))
	sys.stdout.flush()
	return



# --------------------------------------
# ---- FUNCTION TO PRINT TIME STAMP ----
# --------------------------------------

def printProgressTime():
	print ("\n--- %.3f seconds since start" % (time() - t0))
	sys.stdout.flush()
	return



# -----------------------------------------------
# ---- Check if parameter file name provided ----
# -----------------------------------------------

if len(sys.argv) != 2:
	sys.stderr.write("\n\033[1;4mUsage:\033[24m sofia_pipeline.py \033[3m<filename>\033[0m\n\nThe filename of a valid SoFiA parameter file must be specified. Please\nadd the full path if the file is not located in the current directory.\n\n")
	sys.exit(1)



# -----------------------------------------------
# ---- Print some initial status information ----
# -----------------------------------------------

print ("--------------------------")
print ("Running the SoFiA pipeline")
print ("--------------------------")
print ("Using  SoFiA   " + version.getVersion())
print ("       Python  " + str(sys.version_info[0]) + "." + str(sys.version_info[1]) + "." + str(sys.version_info[2]))
print ("       NumPy   " + np.__version__)
print ("       SciPy   " + scipy_version)
print ("       Astropy " + astropy_version)
print ("--------------------------")



# --------------------------------
# ---- START OF SoFiA PIPELINE ---
# --------------------------------

t0 = time()



# ---------------------------------
# ---- READ DEFAULT PARAMETERS ----
# ---------------------------------

printProgressMessage("Reading default parameters")
default_file = os.getenv('SOFIA_PIPELINE_PATH').replace('sofia_pipeline.py', 'SoFiA_default_input.txt')
Parameters = readoptions.readPipelineOptions(default_file)



# ------------------------------
# ---- READ USER PARAMETERS ----
# ------------------------------

printProgressMessage("Reading user parameters")

# This reads in a file with parameters and creates a dictionary:
parameter_file = sys.argv[1]
print ('Parameters extracted from: ' + parameter_file)
print
sys.stdout.flush()
User_Parameters = readoptions.readPipelineOptions(parameter_file)
if not User_Parameters:
	sys.stderr.write("ERROR: No valid parameter settings found in parameter file.\n\n")
	sys.exit(1)

# Overwrite default parameters with user parameters (if exist):
for task in iter(User_Parameters):
	if task in Parameters:
		for key in iter(User_Parameters[task]):
			if key in Parameters[task]:
				Parameters[task][key] = User_Parameters[task][key]

# Define the base and output directory name used for output files (defaults to input file 
# name if writeCat.basename is found to be invalid):
outputBase = Parameters['writeCat']['basename']
outputDir  = Parameters['writeCat']['outputDir']

if outputDir and not os.path.isdir(outputDir):
	sys.stderr.write("ERROR: The specified output directory does not exist:\n")
	sys.stderr.write("       %s\n" % outputDir)
	sys.exit(1)

if not outputBase or outputBase.isspace() or "/" in outputBase or "\\" in outputBase or outputBase == "." or outputBase == "..":
	outroot = Parameters['import']['inFile'].split('/')[-1]
	if (outroot.lower()).endswith(".fits") and len(outroot) > 5:
		outroot = outroot[0:-5]
else:
	outroot = outputBase

if not outputDir or not os.path.isdir(outputDir) or outputDir.isspace():
	outroot = Parameters['import']['inFile'][0:len(Parameters['import']['inFile']) - len(Parameters['import']['inFile'].split('/')[-1])] + outroot
else:
	if outputDir[-1] != '/': outputDir += '/'
	outroot = outputDir + outroot



# -------------------------------------------
# ---- CHECK FOR FILES TO BE OVERWRITTEN ----
# -------------------------------------------


outputFilteredCube  = '%s_filtered.fits' % outroot
outputSkellamPDF    = '%s_skel.pdf' % outroot
outputScatterPDF    = '%s_scat.pdf' % outroot
outputContoursPDF   = '%s_cont.pdf' % outroot
outputMaskCube      = '%s_mask.fits' % outroot
outputMom0Image     = '%s_mom0.fits' % outroot
outputNrchImage     = '%s_nrch.fits' % outroot
outputMom1Image     = '%s_mom1.fits' % outroot
outputCubeletsDir   = '%s/objects/' % outroot
outputCatXml        = '%s_cat.xml' % outroot
outputCatAscii      = '%s_cat.ascii' % outroot
outputCatSQL        = '%s_cat.sql' % outroot
outputCatAsciiDebug = '%s_cat.debug.ascii' % outroot

if not Parameters['writeCat']['overwrite']:
	# Output filtered cube
	if Parameters['steps']['doWriteFilteredCube'] and (Parameters['steps']['doSmooth'] or Parameters['steps']['doScaleNoise'] or Parameters['steps']['doWavelet']):
		checkOverwrite(outputFilteredCube)

	# Reliability plots
	if Parameters['steps']['doReliability'] and Parameters['steps']['doMerge'] and Parameters['reliability']['makePlot']:
		checkOverwrite(outputSkellamPDF)
		checkOverwrite(outputScatterPDF)
		checkOverwrite(outputContoursPDF)

	# Output mask
	if Parameters['steps']['doWriteMask']:
		checkOverwrite(outputMaskCube)
		
	# Moment images
	if Parameters['steps']['doMom0']:
		checkOverwrite(outputMom0Image)
		checkOverwrite(outputNrchImage)
	if Parameters['steps']['doMom1']:
		checkOverwrite(outputMom1Image)

	# Cubelet directory
	if Parameters['steps']['doCubelets'] and Parameters['steps']['doMerge']:
		checkOverwrite(outputCubeletsDir, isFile=False, isDir=True)

	# Output catalogues
	if Parameters['steps']['doWriteCat'] and Parameters['steps']['doMerge'] and Parameters['writeCat']['writeXML']:
		checkOverwrite(outputCatXml)
	if Parameters['steps']['doWriteCat'] and Parameters['steps']['doMerge'] and Parameters['writeCat']['writeASCII']:
		checkOverwrite(outputCatAscii)



# --------------------------------------------------------------
# ---- DEFINE LINKER'S OUTPUT AND CHECK RELIBILITY SETTINGS ----
# --------------------------------------------------------------

if Parameters['steps']['doMerge']:
	# Define parameters returned by the linker module
	catParNames = ('id', 'x_geo', 'y_geo', 'z_geo', 'x', 'y', 'z', 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max', 'n_pix', 'snr_min', 'snr_max', 'snr_sum', 'x_p', 'y_p', 'z_p', 'x_n', 'y_n', 'z_n', 'snr_sum_p', 'snr_sum_n', 'snr_mean', 'snr_std', 'snr_rms', 'w20', 'w50', 'w20_cfd', 'w50_cfd', 'n_x', 'n_y', 'n_chan', 'n_los')
	catParUnits = ('-', 'pix', 'pix', 'chan', 'pix', 'pix', 'chan', 'pix', 'pix', 'pix', 'pix', 'chan', 'chan', '-', '-', '-', '-', 'pix', 'pix', 'chan', 'pix', 'pix', 'chan', '-', '-', '-', '-', '-', 'chan', 'chan', 'chan', 'chan', 'pix', 'pix', 'chan', '-')
	catParFormt = ('%10i', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%7i', '%7i', '%7i', '%7i', '%7i', '%7i', '%8i', '%12.3e', '%12.3e', '%12.3e', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%12.3e', '%12.3e', '%12.3e', '%12.3e', '%12.3e', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%7i', '%7i', '%7i', '%7i')
	
	if Parameters['steps']['doReliability']:
		# Check that the parameters to be used for the reliability calculation are included in catParNames
		for pp in Parameters['reliability']['parSpace']:
			if pp not in catParNames:
				sys.stderr.write("ERROR: You asked to calculate the sources' reliability in the parameter space:\n")
				sys.stderr.write("       " + str(Parameters['reliability']['parSpace']) + "\n")
				sys.stderr.write("       Unfortunately, the parameter %s is not recognised by SoFiA.\n" % (pp))
				sys.stderr.write("       The allowed parameter names are:\n")
				for allowedPar in catParNames: sys.stderr.write("         %s\n" % allowedPar)
				sys.stderr.write("       Please use parameter names from the above list and try again.\n")
				sys.exit(1)



# ---------------------
# ---- IMPORT DATA ----
# ---------------------

printProgressMessage("Reading data cube(s)")
kwargs = Parameters['import'].copy()
kwargs.update(Parameters['flag'])
kwargs.update({'doFlag':Parameters['steps']['doFlag']})
np_Cube, dict_Header, mask, subcube = import_data.read_data(Parameters['steps']['doSubcube'], **kwargs)


# -------------------------
# ---- PRECONDITIONING ----
# -------------------------

if Parameters['steps']['doFlag'] or Parameters['steps']['doSmooth'] or Parameters['steps']['doScaleNoise'] or Parameters['steps']['doWavelet']:
	printProgressMessage("Running input filters")

## ---- FLAGGING ----
#if Parameters['steps']['doFlag']:
	#np_Cube = flag_cube.flag(np_Cube,**Parameters['flag'])

# ---- SMOOTHING ----
if Parameters['steps']['doSmooth']:
	np_Cube = smooth_cube.smooth(np_Cube, **Parameters['smooth'])

# ---- SIGMA CUBE ----
if Parameters['steps']['doScaleNoise']:
	np_Cube = sigma_cube.sigma_scale(np_Cube, **Parameters['scaleNoise'])

# --- WAVELET ---
if Parameters['steps']['doWavelet']:
	print ('Running wavelet filter')
	# WARNING: There is a lot of time and memory overhead from transposing the cube forth and back!
	# WARNING: This will need to be addressed in the future.
	np_Cube = np.transpose(np_Cube, axes=[2, 1, 0])
	np_Cube = wavelet_finder.denoise_2d1d(np_Cube, **Parameters['wavelet'])
	np_Cube = np.transpose(np_Cube, axes=[2, 1, 0])
	np_Cube = np_Cube.copy()


# --- WRITE FILTERED CUBE ---
if Parameters['steps']['doWriteFilteredCube'] and (Parameters['steps']['doSmooth'] or Parameters['steps']['doScaleNoise'] or Parameters['steps']['doWavelet']):
	print ('SoFiA: Writing filtered cube')
	write_filtered_cube.writeFilteredCube(np_Cube, dict_Header, Parameters, outputFilteredCube, Parameters['writeCat']['compress'])

if Parameters['steps']['doFlag'] or Parameters['steps']['doSmooth'] or Parameters['steps']['doScaleNoise'] or Parameters['steps']['doWavelet']:
	print ("Filtering complete")



# -----------------
# ---- FILTERS ----
# -----------------

printProgressMessage("Running source finder")

# Apply the different filters that each create a mask.

# --- PYFIND ---
if Parameters['steps']['doSCfind']:
	print ('Running S+C filter')
	mask |= pyfind.SCfinder_mem(np_Cube, dict_Header, t0, **Parameters['SCfind'])

# --- CNHI ---	
if Parameters['steps']['doCNHI']:
	print ('Running CNHI filter')
	mask = mask + CNHI.find_sources(np_Cube, mask, **Parameters['CNHI'])
 
# --- THRESHOLD ---	
if Parameters['steps']['doThreshold']:
	print ('Running threshold filter')
	threshold_filter.filter(mask, np_Cube, dict_Header, **Parameters['threshold'])

print ('Source finding complete.')

# Check whether positivity flag is set; if so, remove negative pixels from mask:
if Parameters['merge']['positivity']:
	sys.stderr.write("------------------------------------------------------------------------------\n")
	sys.stderr.write("WARNING: Enabling mask.positivity is dangerous and will render some of SoFiA's\n")
	sys.stderr.write("         most  powerful  algorithms useless,  including mask  optimisation and\n")
	sys.stderr.write("         reliability calculation.  Only use this option if you are fully aware\n")
	sys.stderr.write("         of its risks and consequences!\n")
	sys.stderr.write("------------------------------------------------------------------------------\n")
	mask = np.bitwise_and(np.greater(mask, 0), np.greater(np_Cube, 0))

# Check whether any voxel is detected
NRdet = (mask > 0).sum()
if not NRdet:
	sys.stderr.write("WARNING: No voxels detected and included in the mask yet! Exiting pipeline.\n")
	sys.exit()
else:
	print (str(NRdet) +  ' of ' + str(np.array(mask.shape).prod()) + ' pixels detected (' + str(100.0 * float(NRdet) / float(np.array(mask.shape).prod())) + "%)")



# -----------------
# ---- MERGING ----
# -----------------

if Parameters['steps']['doMerge'] and NRdet:
	printProgressMessage("Merging detections")
	objects = []
	objects, mask = linker.link_objects(np_Cube, objects, mask, Parameters['merge']['radiusX'], Parameters['merge']['radiusY'], Parameters['merge']['radiusZ'], Parameters['merge']['minSizeX'], Parameters['merge']['minSizeY'], Parameters['merge']['minSizeZ'])
	if not objects:
		sys.stderr.write("WARNING: No objects remain after merging. Exiting pipeline.\n")
		sys.exit()
	objects = np.array(objects)
	print ('Merging complete')
	NRdet = len(objects)
	NRdetNeg = (np.array(objects)[:,16] < 0).sum()
	print (str(NRdet) + ' sources detected: ' + str(NRdet - NRdetNeg) + ' positive, ' + str(NRdetNeg) + ' negative.')
	# Set catalogue header
	if 'bunit' in dict_Header: dunits = dict_Header['bunit']
	else: dunits = '-'



# -------------------------------------
# ---- OUTPUT FOR DEBUGGING (MASK) ----
# -------------------------------------

if Parameters['steps']['doDebug'] and NRdet:
	printProgressMessage("Writing all-source mask cube for debugging")
	writemask.writeMask(mask, dict_Header, Parameters, '%s_mask.debug_all.fits' % outroot, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])


# ----------------------------------------------------
# ---- ESTIMATE RELIABILITY FROM NEGATIVE SOURCES ----
# ----------------------------------------------------

if Parameters['steps']['doReliability'] and Parameters['steps']['doMerge'] and NRdet and not NRdetNeg:
	printProgressTime()
	sys.stderr.write("------------------------------------------------------------------------------\n")
	sys.stderr.write("ERROR: You asked SoFiA to calculate  the reliability  of the detected sources.\n")
	sys.stderr.write("       Unfortunately,  this calculation  cannot be done  because there  are no\n")
	sys.stderr.write("       negative sources in the catalogue of detections. This may occur because\n")
	sys.stderr.write("       of your source-finding and/or filtering settings.\n")
	sys.stderr.write("       Negative sources  are strictly necessary  to calculate the reliability.\n")
	sys.stderr.write("       You can do one of the following:\n")
	sys.stderr.write("       (1) Switch off the reliability calculation.\n")
	sys.stderr.write("       (2) Modify the source-finding and/or filtering settings in order to de-\n")
	sys.stderr.write("           tect negative sources.\n")
	sys.stderr.write("------------------------------------------------------------------------------\n")
	sys.exit(1)

elif Parameters['steps']['doReliability'] and Parameters['steps']['doMerge'] and NRdet and NRdetNeg:
	# ---- MEASURE GLOBAL SIGMA AND NORMALISE PARAMETERS----
	printProgressMessage("Measuring cube noise to divide flux parameters by global rms")
	maxNrVox = 1e+6 # maximum nr of voxels over which to calculate the global RMS. Sampling below is set accordingly.
	sampleRms = max(1, int((float(np.array(np_Cube.shape).prod()) / maxNrVox)**(1.0 / min(3, len(np_Cube.shape)))))
	globalrms = functions.GetRMS(np_Cube, rmsMode='negative', zoomx=1, zoomy=1, zoomz=1, verbose=True, sample=sampleRms)
	printProgressTime()

	# normalise flux parameters to global rms
	# (this is done also if weights were applied, in case they are prop. to 1/sigma but not exactly = 1/sigma)
	objects = np.array(objects)
	objects[:,catParNames.index('snr_min')] /= globalrms
	objects[:,catParNames.index('snr_max')] /= globalrms
	objects[:,catParNames.index('snr_sum')] /= globalrms
	objects = [list(item) for item in list(objects)]

	# ---- CALCULATE RELIABILITY ----
	printProgressMessage("Determining reliability")
	objects, reliable = addrel.EstimateRel(np.array(objects), outroot, catParNames, **Parameters['reliability'])
	print ('The following sources have been detected: ' + str(reliable))
	catParNames = tuple(list(catParNames) + ['n_pos',  'n_neg',  'rel'])
	catParUnits = tuple(list(catParUnits) + ['-',      '-',      '-'])
	catParFormt = tuple(list(catParFormt) + ['%12.3e', '%12.3e', '%12.6f'])

elif Parameters['steps']['doMerge'] and NRdet:
	printProgressTime()
	reliable = list(np.array(objects)[np.array(objects)[:,16] > 0,0].astype(int)) # select all positive sources
	print ('The following sources have been detected: ' + str(reliable))

else:
	printProgressTime()
	reliable = [1,] # if not merging, all detected voxels have ID = 1 and here they are set to be reliable



# ------------------------------------------
# ---- OUTPUT FOR DEBUGGING (CATALOGUE) ----
# ------------------------------------------

if Parameters['steps']['doDebug']:
	printProgressMessage("Writing all-source debugging catalogue including all parameters relevant for the reliability calculation")
	write_catalog.write_catalog_from_array('ASCII', objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'], outputCatAsciiDebug, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'], Parameters['parameters']['getUncertainties'])



# ------------------------------------------------------
# ---- REMOVE UNNECESSARY PARAMETERS FROM CATALOGUE ----
# ------------------------------------------------------

if Parameters['steps']['doMerge'] and NRdet:
	objects, catParNames, catParUnits, catParFormt = np.array(objects), list(catParNames), list(catParUnits), list(catParFormt)
	removecols = ['snr_min', 'snr_max', 'snr_sum', 'x_p', 'y_p', 'z_p', 'x_n', 'y_n', 'z_n', 'snr_sum_p', 'snr_sum_n', 'snr_mean', 'snr_std', 'snr_rms', 'w20', 'w50', 'w20_cfd', 'w50_cfd', 'n_pos', 'n_neg', 'n_x', 'n_y']
	for remcol in removecols:
		if remcol in catParNames:
			index = catParNames.index(remcol)
			del(catParNames[index])
			del(catParUnits[index])
			del(catParFormt[index])
			objects = np.delete(objects, [index], axis=1)
	objects, catParNames, catParUnits, catParFormt = [list(item) for item in list(objects)], tuple(catParNames), tuple(catParUnits), tuple(catParFormt)



# --------------------------------------------------
# ---- REMOVE NON RELIABLE AND NEGATIVE SOURCES ----
# --------------------------------------------------

if Parameters['steps']['doMerge'] and NRdet:
	printProgressMessage("Removing unreliable sources")

	# make sure that reliable is sorted
	relList = list(reliable)
	relList.sort()
	reliable = np.array(relList)

	# remove non reliable sources in the objects array
	relObjects = []
	for rr in reliable:
		relObjects.append([len(relObjects) + 1] + list(objects[rr - 1]))
	relObjects = np.array(relObjects)
	objects = relObjects

	tmpCatParNames = list(catParNames);
	#tmpCatParNames[0] = 'id_old'   WARNING: This is fatal, because it implicitly assumes that the first entry is 'id'!!!
	#                                        New items must always be INSERTED at a particular position without ever 
	#                                        OVERWRITING existing items in the list!
	# Better:
	tmpCatParNames.insert(1, "id_old");
	catParNames = tuple(tmpCatParNames);

	tmpCatParFormt = list(catParFormt);
	tmpCatParFormt.insert(1, "%10i");
	catParFormt= tuple(tmpCatParFormt);

	tmpCatParUnits = list(catParUnits);
	tmpCatParUnits.insert(1, "-");
	catParUnits= tuple(tmpCatParUnits);

	# in the mask file
	mask *= -1
	index = 1
	catParNames = np.array(catParNames)
	for rr in reliable:
		objrr = objects[objects[:,1] == rr][0]
		Xmin  = int(objrr[catParNames == 'x_min'])
		Ymin  = int(objrr[catParNames == 'y_min'])
		Zmin  = int(objrr[catParNames == 'z_min'])
		Xmax  = int(objrr[catParNames == 'x_max'])
		Ymax  = int(objrr[catParNames == 'y_max'])
		Zmax  = int(objrr[catParNames == 'z_max'])
		mask[Zmin:Zmax+1, Ymin:Ymax+1, Xmin:Xmax+1][mask[Zmin:Zmax+1, Ymin:Ymax+1, Xmin:Xmax+1] == -rr] = index
		index += 1
	mask[mask < 0] = 0
	catParNames = tuple(catParNames)


	newRel = []
	for i in range(0, len(relObjects)):
		newRel.append(i + 1)
	reliable = np.array(newRel)
	NRdet = objects.shape[0]



# -------------------------------------------------------------------------------
# ---- RELOAD ORIGINAL DATA CUBE FOR PARAMETERISATION IF IT HAS BEEN CHANGED ----
# -------------------------------------------------------------------------------

if Parameters['steps']['doSmooth'] or Parameters['steps']['doScaleNoise'] or Parameters['import']['weightsFile'] or Parameters['import']['weightsFunction']:
	Parameters['import']['weightsFile'] = ''
	Parameters['import']['maskFile'] = ''
	Parameters['import']['weightsFunction'] = ''
	del np_Cube, dict_Header
	kwargs = Parameters['import'].copy()
	kwargs.update({'cubeOnly':True})
	np_Cube, dict_Header = import_data.read_data(Parameters['steps']['doSubcube'], **kwargs)



# ----------------------------------------
# ---- OUTPUT FOR DEBUGGING (MOMENTS) ----
# ----------------------------------------

if Parameters['steps']['doDebug'] and NRdet:
	printProgressMessage("Writing pre-optimisation mask and moment maps for debugging")
	debug = 1
	#writemask.writeMask(mask, dict_Header, Parameters, '%s_mask.debug_rel.fits'%outroot,Parameters['writeCat']['compress'])
	#mom0_Image = writemoment2.writeMoment0(np_Cube, mask, outroot, debug, dict_Header,Parameters['writeCat']['compress'])
	#writemoment2.writeMoment1(np_Cube, mask, outroot, debug, dict_Header, mom0_Image,Parameters['writeCat']['compress'])



# ----------------------
# ---- PARAMETERISE ----
# ----------------------

if Parameters['steps']['doParameterise'] and Parameters['steps']['doMerge'] and NRdet:
	printProgressMessage("Parameterising sources")
	
	# Print warning message about statistical uncertainties
	if Parameters['parameters']['getUncertainties']:
		sys.stderr.write("------------------------------------------------------------\n")
		sys.stderr.write("WARNING:    You have requested statistical uncertainties for\n")
		sys.stderr.write("         several source parameters. Please be aware that the\n")
		sys.stderr.write("         calculation of statistical uncertainties depends on\n")
		sys.stderr.write("         a number of assumptions that may not be met. Hence,\n")
		sys.stderr.write("         the resulting numbers  may not be representative of\n")
		sys.stderr.write("         the true uncertainties of those parameters, in par-\n")
		sys.stderr.write("         ticular in the presence of systematic errors.\n")
		sys.stderr.write("------------------------------------------------------------\n")
	
	if Parameters['parameters']['dilateMask']: mask = parametrisation.dilate(np_Cube, mask, objects, catParNames, Parameters)
	np_Cube, mask, objects, catParNames, catParFormt, catParUnits = parametrisation.parametrise(np_Cube, mask, objects, catParNames, catParFormt, catParUnits, Parameters, dunits)
	catParNames = tuple(catParNames)
	catParUnits = tuple(catParUnits)
	catParFormt = tuple(catParFormt)

	print ('Parameterisation complete')



# --------------------
# ---- WRITE MASK ----
# --------------------

if Parameters['steps']['doWriteMask'] and NRdet:
	printProgressMessage("Writing mask cube")
	writemask.writeMask(mask, dict_Header, Parameters, outputMaskCube, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])


	
# ------------------------
# ---- STORE CUBELETS ----
# ------------------------

if Parameters['steps']['doCubelets'] and Parameters['steps']['doMerge'] and NRdet:
	printProgressMessage("Writing cubelets")
	objects = np.array(objects)
	cathead = np.array(catParNames)
	cubelets.writeSubcube(np_Cube, dict_Header, mask, objects, cathead, outroot, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])


# ----------------------------
# ---- MAKE MOM0 and MOM1 ----
# ----------------------------

if (Parameters['steps']['doMom0'] or Parameters['steps']['doMom1']) and NRdet:
	printProgressMessage("Writing moment maps")
	debug = 0
	writemoment2.writeMoments(np_Cube, mask, outroot, debug, dict_Header,Parameters['writeCat']['compress'], Parameters['steps']['doMom0'], Parameters['steps']['doMom1'], Parameters['writeCat']['overwrite'])

	## read the original data cube again if it is needed for further steps
	#print "\n--- SoFiA: Reloading original data cube ---"
	#np_Cube, dict_Header, mask, subcube = import_data.read_data(Parameters['steps']['doSubcube'], **Parameters['import'])


# ----------------------------------------------------
# ---- CORRECT COORDINATES IF WORKING ON SUBCUBES ----
# ----------------------------------------------------

if len(subcube) and Parameters['steps']['doMerge'] and NRdet:
	printProgressMessage("Correcting parameters for sub-cube offset")
	# list of parameters to correct for X, Y and Z offset
	corrX=['x_geo', 'x', 'x_min', 'x_max']
	corrY=['y_geo', 'y', 'y_min', 'y_max']
	corrZ=['z_geo', 'z', 'z_min', 'z_max', 'bf_z']

	if subcube[0]:
		for pp in corrX:
			if pp in catParNames: objects[:,list(catParNames).index(pp)] += subcube[0]
	if subcube[2]:
		for pp in corrY:
			if pp in catParNames: objects[:,list(catParNames).index(pp)] += subcube[2]
	if subcube[4]:
		for pp in corrZ:
			if pp in catParNames: objects[:,list(catParNames).index(pp)] += subcube[4]



# ---------------------------------------------------
# ---- APPEND PARAMETER VALUES IN PHYSICAL UNITS ----
# ---------------------------------------------------

if Parameters['steps']['doMerge'] and NRdet and Parameters['steps']['doWriteCat']:
	printProgressMessage("Adding WCS position to catalogue")
	objects, catParNames, catParFormt, catParUnits = wcs_coordinates.add_wcs_coordinates(objects, catParNames, catParFormt, catParUnits, Parameters)



# --------------------
# ---- STORE DATA ----
# --------------------

if Parameters['steps']['doWriteCat'] and Parameters['steps']['doMerge'] and NRdet:
	printProgressMessage("Writing output catalogue")

	if 'rms' in catParNames:
		catParFormt=list(catParFormt)
		catParFormt[list(catParNames).index('rms')] = '%12.4e'
		catParFormt=tuple(catParFormt)

	if Parameters['writeCat']['writeXML'] and Parameters['steps']['doMerge'] and NRdet:
		write_catalog.write_catalog_from_array('XML', objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'], outputCatXml, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'], Parameters['parameters']['getUncertainties'])

	if Parameters['writeCat']['writeASCII'] and Parameters['steps']['doMerge'] and NRdet:
		write_catalog.write_catalog_from_array('ASCII', objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'], outputCatAscii, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'], Parameters['parameters']['getUncertainties'])

	if Parameters['writeCat']['writeSQL'] and Parameters['steps']['doMerge'] and NRdet:
		write_catalog.write_catalog_from_array('SQL', objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'], outputCatSQL, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'], Parameters['parameters']['getUncertainties'])



printProgressMessage("Pipeline finished")
