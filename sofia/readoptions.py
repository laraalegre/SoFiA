#! /usr/bin/env python
import re
import sys
import traceback
import ast


# -----------------------------------------
# FUNCTION: Convert string to Boolean value
# -----------------------------------------

def str2bool(s):
	return s in ["True", "true", "TRUE", "Yes", "yes", "YES"]


# -------------------------------------------
# FUNCTION: Read settings from parameter file
# -------------------------------------------

def readPipelineOptions(filename = "pipeline.options"):
	# Try to open parameter file
	try:
		f = open(filename, "r")
	except IOError as e:
		sys.stderr.write("FATAL ERROR: Failed to read parameter file: " + str(filename) + "\n")
		sys.stderr.write(str(e) + "\n")
		sys.exit(1);
	
	# Extract lines from parameter file
	lines = f.readlines()
	f.close()
	
	# Remove leading/trailing whitespace, empty lines and comments
	lines = [line.strip() for line in lines]
	lines = [line for line in lines if len(line) > 0 and line[0] != "#"]
	
	# Some additional setup
	datatypes = allowedDataTypes()
	tasks = {}
	
	# Loop over all lines:
	for line in lines:
		# Extract parameter name and value
		try:
			parameter, value = tuple(line.split("=", 1))
			parameter = parameter.strip()
			value = value.split("#")[0].strip()
			module, parname = tuple(parameter.split(".", 1))
		except:
			sys.stderr.write("FATAL ERROR: Failed to read parameter: " + str(line) + "\n")
			sys.stderr.write("             Expected format: module.parameter = value\n")
			sys.exit(1)
		
		# Ensure that module and parameter names are not empty
		if len(module) < 1 or len(parname) < 1:
			sys.stderr.write("FATAL ERROR: Failed to read parameter: " + str(line) + "\n")
			sys.stderr.write("             Expected format: module.parameter = value\n")
			sys.exit(1)
		
		subtasks = tasks
		if module not in subtasks: subtasks[module] = {}
		subtasks = subtasks[module]
		
		if parname in subtasks:
			sys.stderr.write("WARNING: Multiple definitions of parameter " + str(parameter) + " encountered.\n")
			sys.stderr.write('         Ignoring all additional definitions.\n')
			continue
		
		if parameter in datatypes:
			try:
				if datatypes[parameter]   == "bool":  subtasks[parname] = str2bool(value)
				elif datatypes[parameter] == "float": subtasks[parname] = float(value)
				elif datatypes[parameter] == "int":   subtasks[parname] = int(value)
				elif datatypes[parameter] == "array": subtasks[parname] = ast.literal_eval(value)
				else: subtasks[parname] = str(value)
			except:
				sys.stderr.write("FATAL ERROR: Failed to parse parameter value:\n")
				sys.stderr.write("             " + str(line) + "\n")
				sys.stderr.write("             Expected data type: " + str(datatypes[parameter]) + "\n")
				sys.exit(1)
		else:
			sys.stderr.write("WARNING: Ignoring unknown parameter: " + str(parameter) + " = " + str(value) + "\n")
			continue
	
	return tasks


# ----------------------------------------------------------
# FUNCTION: Return list of allowed parameter names and types
# ----------------------------------------------------------

def allowedDataTypes():
	# Define data types of individual parameters. This is used to
	# convert the input values into the correct data type.
	# Ensure that all new parameters get added to this list!
	# Parameters not listed here will be ignored by the parser!
	return {"steps.doSubcube": "bool", \
	        "steps.doFlag": "bool", \
	        "steps.doSmooth": "bool", \
	        "steps.doScaleNoise": "bool", \
	        "steps.doSCfind": "bool", \
	        "steps.doThreshold": "bool", \
	        "steps.doWavelet": "bool", \
	        "steps.doCNHI": "bool", \
	        "steps.doMerge": "bool", \
	        "steps.doReliability": "bool", \
	        "steps.doParameterise": "bool", \
	        "steps.doWriteFilteredCube": "bool", \
	        "steps.doWriteNoiseCube": "bool", \
	        "steps.doWriteMask": "bool", \
	        "steps.doWriteCat": "bool", \
	        "steps.doMom0": "bool", \
	        "steps.doMom1": "bool", \
	        "steps.doCubelets": "bool", \
	        "steps.doDebug": "bool", \
	        "steps.doOptical": "bool", \
	        "import.inFile": "string", \
	        "import.weightsFile": "string", \
	        "import.maskFile": "string", \
	        "import.weightsFunction": "string", \
	        "import.subcubeMode": "string", \
	        "import.subcube": "array", \
	        "flag.regions": "array",\
	        "flag.file": "string", \
	        "optical.sourceCatalogue": "string", \
	        "optical.spatSize": "float", \
	        "optical.specSize": "float", \
	        "optical.storeMultiCat": "bool", \
	        "smooth.kernel": "string", \
	        "smooth.edgeMode": "string", \
	        "smooth.kernelX": "float", \
	        "smooth.kernelY": "float", \
	        "smooth.kernelZ": "float", \
	        "scaleNoise.method": "string", \
	        "scaleNoise.statistic": "string", \
	        "scaleNoise.fluxRange": "string", \
	        "scaleNoise.edgeX": "int", \
	        "scaleNoise.edgeY": "int", \
	        "scaleNoise.edgeZ": "int", \
	        "scaleNoise.scaleX": "bool", \
	        "scaleNoise.scaleY": "bool", \
	        "scaleNoise.scaleZ": "bool", \
	        "scaleNoise.windowSpatial": "int", \
	        "scaleNoise.windowSpectral": "int", \
	        "scaleNoise.gridSpatial": "int", \
	        "scaleNoise.gridSpectral": "int", \
	        "scaleNoise.interpolation": "string", \
	        "SCfind.threshold": "float", \
	        "SCfind.sizeFilter": "float", \
	        "SCfind.maskScaleXY": "float", \
	        "SCfind.maskScaleZ": "float", \
	        "SCfind.edgeMode": "string", \
	        "SCfind.rmsMode": "string", \
	        "SCfind.fluxRange": "string", \
	        "SCfind.kernels": "array", \
	        "SCfind.kernelUnit": "string", \
	        "SCfind.verbose": "bool", \
	        "CNHI.pReq": "float", \
	        "CNHI.qReq": "float", \
	        "CNHI.minScale": "int", \
	        "CNHI.maxScale": "int", \
	        "CNHI.medianTest": "bool", \
	        "CNHI.verbose": "int", \
	        "wavelet.threshold": "float", \
	        "wavelet.scaleXY": "int", \
	        "wavelet.scaleZ": "int", \
	        "wavelet.positivity": "bool", \
	        "wavelet.iterations": "int", \
	        "threshold.threshold": "float", \
	        "threshold.clipMethod": "string", \
	        "threshold.rmsMode": "string", \
	        "threshold.fluxRange": "string", \
	        "threshold.verbose": "bool", \
	        "merge.radiusX": "int", \
	        "merge.radiusY": "int", \
	        "merge.radiusZ": "int", \
	        "merge.minSizeX": "int", \
	        "merge.minSizeY": "int", \
	        "merge.minSizeZ": "int", \
	        "merge.positivity": "bool", \
	        "reliability.parSpace": "array", \
	        "reliability.logPars": "array", \
	        "reliability.autoKernel": "bool", \
	        "reliability.scaleKernel": "float", \
	        "reliability.usecov": "bool", \
	        "reliability.negPerBin": "float", \
	        "reliability.skellamTol": "float", \
	        "reliability.kernel": "array", \
	        "reliability.fMin": "float", \
	        "reliability.threshold": "float", \
	        "reliability.makePlot": "bool", \
	        "parameters.getUncertainties": "bool", \
	        "parameters.fitBusyFunction": "bool", \
	        "parameters.optimiseMask": "bool", \
	        "parameters.dilateMask": "bool", \
	        "parameters.dilateThreshold":"float",\
	        "parameters.dilatePixMax":"int",\
	        "parameters.dilateChan":"int",\
	        "writeCat.overwrite":"bool", \
	        "writeCat.compress":"bool", \
	        "writeCat.outputDir": "string", \
	        "writeCat.basename": "string", \
	        "writeCat.writeASCII": "bool", \
	        "writeCat.writeXML": "bool", \
	        "writeCat.writeSQL": "bool", \
	        "writeCat.parameters": "array" }
