#! /usr/bin/env python

# import default python libraries
import numpy as np
import sys, os
import string
from time import time

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
from sofia import linker
from sofia import store_xml
from sofia import store_ascii
from sofia import cubelets
from sofia import parametrisation
from sofia import wcs_coordinates
from sofia import flag_cube
from sofia import CNHI

t0=time()

# ---------------------------------
# ---- READ DEFAULT PARAMETERS ----
# ---------------------------------

print "\n--- SoFiA: Reading default parameters ---"
sys.stdout.flush()

# This reads in the default parameters:
default_file = '%s/SoFiA_default_input.txt'%(os.path.dirname(os.path.realpath(__file__)))
#default_file = 'SoFiA_default_input.txt'
Parameters = readoptions.readPipelineOptions(default_file)



# ------------------------------
# ---- READ USER PARAMETERS ----
# ------------------------------

print "\n--- %.3f seconds since start"%(time()-t0)
print "\n--- SoFiA: Reading user parameters ---"
sys.stdout.flush()

# This reads in a file with parameters and creates a dictionary:
parameter_file = sys.argv[1]
print 'Parameters extracted from: ', parameter_file
print
User_Parameters = readoptions.readPipelineOptions(parameter_file)

# Overwrite default parameters with user parameters (if exist):
for task in User_Parameters.iterkeys():
        if(task in Parameters):
                for key in User_Parameters[task].iterkeys():
                        if(key in Parameters[task]):
                                Parameters[task][key] = User_Parameters[task][key]

# Define the base and output directory name used for output files (defaults to input file 
# name if writeCat.basename is found to be invalid):
outputBase = Parameters['writeCat']['basename']
outputDir  = Parameters['writeCat']['outputDir']

if((not outputBase) or outputBase.isspace() or ("/" in outputBase) or ("\\" in outputBase) or (outputBase == ".") or (outputBase == "..")):
    outroot = Parameters['import']['inFile'].split('/')[-1]
    if((outroot.lower()).endswith(".fits") and len(outroot) > 5):
        outroot = outroot[0:-5]
else:
    outroot = outputBase

if((not outputDir) or (not os.path.isdir(outputDir)) or (outputDir.isspace())):
    outroot = Parameters['import']['inFile'][0:len(Parameters['import']['inFile'])-len(Parameters['import']['inFile'].split('/')[-1])]+outroot
else:
    if outputDir[-1] != '/': outputDir += '/'
    outroot = outputDir + outroot


# ---------------------
# ---- IMPORT DATA ----
# ---------------------

print "\n--- %.3f seconds since start"%(time()-t0)
print "\n--- SoFiA: Reading data cube(s) ---"
sys.stdout.flush()

np_Cube, dict_Header, mask, subcube = import_data.read_data(Parameters['steps']['doSubcube'],**Parameters['import'])


# -------------------------
# ---- PRECONDITIONING ----
# -------------------------

print "\n--- %.3f seconds since start"%(time()-t0)
print "\n--- SoFiA: Running input filters ---"
sys.stdout.flush()

# ---- FLAGGING ----
if Parameters['steps']['doFlag']:
	np_Cube = flag_cube.flag(np_Cube,**Parameters['flag'])

# ---- SMOOTHING ----
if Parameters['steps']['doSmooth']:
	np_Cube = smooth_cube.smooth(np_Cube, **Parameters['smooth'])	

# ---- RFI ----


# ---- SIGMA CUBE ----
if Parameters['steps']['doScaleNoise']:
	np_Cube = sigma_cube.sigma_scale(np_Cube, **Parameters['scaleNoise'])

globalrms=functions.GetRMS(np_Cube,rmsMode='negative',zoomx=1,zoomy=1,zoomz=1,verbose=True)

# --- WAVELET ---
if Parameters['steps']['doWavelet']:
        print 'Running wavelet filter'
        # WARNING: There is a lot of time and memory overhead from transposing the cube forth and back!
        # WARNING: This will need to be addressed in the future.
        np_Cube = np.transpose(np_Cube, axes=[2, 1, 0])
        np_Cube = wavelet_finder.denoise_2d1d(np_Cube, **Parameters['wavelet'])
        np_Cube = np.transpose(np_Cube, axes=[2, 1, 0])

# --- WRITE FILTERED CUBE ---
if Parameters['steps']['doWriteFilteredCube'] and (Parameters['steps']['doSmooth'] or Parameters['steps']['doScaleNoise'] or Parameters['steps']['doWavelet']):
        print "SoFiA: Writing filtered cube"
        write_filtered_cube.writeFilteredCube(np_Cube, dict_Header, Parameters, '%s_filtered.fits' % outroot, Parameters['writeCat']['compress'])

print "Filtering complete"
print 

# -----------------
# ---- FILTERS ----
# -----------------

print "\n--- %.3f seconds since start"%(time()-t0)
print "\n--- SoFiA: Running source finder ---"
sys.stdout.flush()

# apply the different filters that each create a mask.
# create an empty mask, the size of the cube:


# --- PYFIND ---
if Parameters['steps']['doSCfind']:
	print 'Running S+C filter'
	pyfind_mask = pyfind.SCfinder_mem(np_Cube, dict_Header, t0, **Parameters['SCfind'])
	mask = mask + pyfind_mask

# --- CNHI ---	
if Parameters['steps']['doCNHI']:
	print 'Running CNHI filter'
	mask = mask + CNHI.find_sources(np_Cube, mask, **Parameters['CNHI'])
 
# --- THRESHOLD ---	
if Parameters['steps']['doThreshold']:
	print 'Running threshold filter'
	#mask+=threshold_filter.filter(np_Cube, dict_Header, **Parameters['threshold'])
	threshold_filter.filter(mask, np_Cube, dict_Header, **Parameters['threshold'])

print 'Source finding complete.'
print

# Check whether any voxel is detected
NRdet = (mask > 0).sum()
if not NRdet:
	sys.stderr.write("WARNING: No voxels detected! EXITING pipeline.\n")
	print 
	sys.exit()



# -----------------
# ---- MERGING ----
# -----------------

if Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Merging detections ---"
	sys.stdout.flush()
	objects, mask = linker.link_objects(np_Cube, mask, **Parameters['merge'])
	if not objects:
		sys.stderr.write("WARNING: No objects remain after merging. Exiting pipeline.\n")
		sys.exit()
	print 'Merging complete'
	print
	NRdet = len(objects)
	# set catalog header	
	if 'bunit' in dict_Header: dunits=dict_Header['bunit']
	else: dunits='-'
	catParNames = ('id','x_geo','y_geo','z_geo','x','y','z','x_min','x_max','y_min','y_max','z_min','z_max','n_pix','snr_min','snr_max','snr_sum')
	catParUnits = ('-','pix','pix','chan','pix','pix','chan','pix','pix','pix','pix','chan','chan','-','-','-','-')
	catParFormt = ('%10i', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%10.3f', '%7i', '%7i', '%7i', '%7i', '%7i', '%7i', '%8i', '%12.3e', '%12.3e', '%12.3e')
	# normalise flux parameters to global rms (this is done also if weights were applied, in case they are prop. to 1/sigma but not exactly = 1/sigma)
	print 'Dividing flux parameters by global cube rms'
	objects=np.array(objects)
	objects[:,catParNames.index('snr_min')]/=globalrms
	objects[:,catParNames.index('snr_max')]/=globalrms
	objects[:,catParNames.index('snr_sum')]/=globalrms
	objects=[list(jj) for jj in list(objects)]


# -------------------------------------
# ---- OUTPUT FOR DEBUGGING (MASK) ----
# -------------------------------------

if Parameters['steps']['doDebug'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing all-source mask cube for debugging ---"
	sys.stdout.flush()
	writemask.writeMask(mask, dict_Header, Parameters, '%s_mask.debug_all.fits'%outroot,Parameters['writeCat']['compress'])



# ----------------------------------------------------
# ---- ESTIMATE RELIABILITY FROM NEGATIVE SOURCES ----
# ----------------------------------------------------

if Parameters['steps']['doReliability'] and Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Determining reliability ---"
	sys.stdout.flush()
	objects,reliable = addrel.EstimateRel(np.array(objects), outroot, catParNames, **Parameters['reliability'])
	print 'The following sources have been detected:', reliable
	print
	catParNames = tuple(list(catParNames) + ['n_pos',  'n_neg',  'rel'])
	catParUnits = tuple(list(catParUnits) + ['-','-','-'])
	catParFormt = tuple(list(catParFormt) + ['%12.3e', '%12.3e', '%12.6f'])
elif Parameters['steps']['doMerge'] and NRdet:
	reliable = list(np.array(objects)[np.array(objects)[:,16] > 0,0].astype(int)) # select all positive sources
	print 'The following sources have been detected:', reliable
else: reliable=[1,] # if not merging, all detected voxels have ID = 1 and here they are set to be reliable



# ------------------------------------------
# ---- OUTPUT FOR DEBUGGING (CATALOGUE) ----
# ------------------------------------------

if Parameters['steps']['doDebug']:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing all-source catalogue for debugging ---"
	#sys.stdout.flush()
	store_ascii.make_ascii_from_array(objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'], outroot+'_cat.debug.ascii',Parameters['writeCat']['compress'])



# --------------------------------------------------
# ---- REMOVE NON RELIABLE AND NEGATIVE SOURCES ----
# --------------------------------------------------

if Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Removing unreliable sources ---"
	sys.stdout.flush()

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
	for rr in reliable:
        	mask[mask == -rr] = index
        	index += 1
	mask[mask < 0] = 0


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
	np_Cube, dict_Header, bla, blabla = import_data.read_data(Parameters['steps']['doSubcube'],**Parameters['import'])



# ----------------------------------------
# ---- OUTPUT FOR DEBUGGING (MOMENTS) ----
# ----------------------------------------

if Parameters['steps']['doDebug'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing pre-optimisation mask and moment maps for debugging ---"
	sys.stdout.flush()
	debug=1
	writemask.writeMask(mask, dict_Header, Parameters, '%s_mask.debug_rel.fits'%outroot,Parameters['writeCat']['compress'])
	mom0_Image = writemoment2.writeMoment0(np_Cube, mask, outroot, debug, dict_Header,Parameters['writeCat']['compress'])
	writemoment2.writeMoment1(np_Cube, mask, outroot, debug, dict_Header, mom0_Image,Parameters['writeCat']['compress'])



# ----------------------
# ---- PARAMETERISE ----
# ----------------------

#results=PyCatalog.PySourceCatalog()
#results.readDuchampFile("duchamp-results.txt")

if Parameters['steps']['doParameterise'] and Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Parameterising sources ---"
	sys.stdout.flush()
#	np_Cube, dict_Header, mask, objects, catParNames, catParFormt = parametrisation.parametrise(np_Cube, dict_Header, mask, objects, catParNames, catParFormt, Parameters)
	if Parameters['parameters']['dilateMask']: mask = parametrisation.dilate(np_Cube,mask,objects,catParNames,Parameters)
	np_Cube, mask, objects, catParNames, catParFormt, catParUnits = parametrisation.parametrise(np_Cube, mask, objects, catParNames, catParFormt, catParUnits, Parameters, dunits)
	catParNames=tuple(catParNames)
	catParUnits=tuple(catParUnits)
	catParFormt=tuple(catParFormt)

        
	print 'Parameterisation complete'
	print



# ----------------------------------------------------
# ---- CORRECT COORDINATES IF WORKING ON SUBCUBES ----
# ----------------------------------------------------

if len(subcube) and Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Correcting parameters for sub-cube offset ---"
	sys.stdout.flush()
	# list of parameters to correct for X, Y and Z offset
	corrX=['x_geo','x','x_min','x_max']
	corrY=['y_geo','y','y_min','y_max']
	corrZ=['z_geo','z','z_min','z_max','bf_z']

	if subcube[0]:
		for pp in corrX:
			if pp in catParNames: objects[:,list(catParNames).index(pp)]+=subcube[0]
	if subcube[2]:
		for pp in corrY:
			if pp in catParNames: objects[:,list(catParNames).index(pp)]+=subcube[2]
	if subcube[4]:
		for pp in corrZ:
			if pp in catParNames: objects[:,list(catParNames).index(pp)]+=subcube[4]



# --------------------
# ---- WRITE MASK ----
# --------------------

if Parameters['steps']['doWriteMask'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing mask cube ---"
	sys.stdout.flush()
	writemask.writeMask(mask, dict_Header, Parameters, '%s_mask.fits'%outroot, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])



# ----------------------------
# ---- MAKE MOM0 and MOM1 ----
# ----------------------------

if Parameters['steps']['doMom0'] or Parameters['steps']['doMom1']:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing moment maps ---"
	sys.stdout.flush()
	debug = 0
	writemoment2.writeMoments(np_Cube, mask, outroot, debug, dict_Header,Parameters['writeCat']['compress'], Parameters['steps']['doMom0'], Parameters['steps']['doMom1'], Parameters['writeCat']['overwrite'])


# --------------------
# ---- MAKE MOM1  ----
# --------------------
#
#if Parameters['steps']['doMom1'] and NRdet:
#	print "\n--- %.3f seconds since start"%(time()-t0)
#	print "\n--- SoFiA: Writing moment-1 map ---"
#	sys.stdout.flush()
#	debug = 0
#	writemoment2.writeMoment1(np_Cube, mask, outroot, debug, dict_Header, mom0_Image,Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])



# ------------------------
# ---- STORE CUBELETS ----
# ------------------------

if Parameters['steps']['doCubelets'] and Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing cubelets ---"
	sys.stdout.flush()
	objects = np.array(objects)
	cathead = np.array(catParNames)
	cubelets.writeSubcube(np_Cube, dict_Header, mask, objects, cathead, outroot, Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])



# ---------------------------------------------------
# ---- APPEND PARAMETER VALUES IN PHYSICAL UNITS ----
# ---------------------------------------------------

if Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Adding WCS position to catalogue ---"
	sys.stdout.flush()
	objects, catParNames, catParFormt, catParUnits = wcs_coordinates.add_wcs_coordinates(objects,catParNames,catParFormt,catParUnits,Parameters)



# --------------------
# ---- STORE DATA ----
# --------------------

if Parameters['steps']['doWriteCat'] and Parameters['steps']['doMerge'] and NRdet:
	print "\n--- %.3f seconds since start"%(time()-t0)
	print "\n--- SoFiA: Writing output catalogue ---"
	sys.stdout.flush()
	if Parameters['writeCat']['writeXML'] and Parameters['steps']['doMerge'] and NRdet:
		store_xml.make_xml_from_array(objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'],outroot + '_cat.xml',Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])
		#store_xml.make_xml(results, outroot + '_cat.xml', Parameters['writeCat']['overwrite'])
	if Parameters['writeCat']['writeASCII'] and Parameters['steps']['doMerge'] and NRdet:
		store_ascii.make_ascii_from_array(objects, catParNames, catParUnits, catParFormt, Parameters['writeCat']['parameters'], outroot+'_cat.ascii',Parameters['writeCat']['compress'], Parameters['writeCat']['overwrite'])
		#store_ascii.make_ascii(results, Parameters['writeCat']['parameters'], outroot + '_cat.ascii', Parameters['writeCat']['overwrite'])



print "\n--- %.3f seconds since start"%(time()-t0)
print "\n--- SoFiA: Pipeline finished ---"
sys.stdout.flush()
