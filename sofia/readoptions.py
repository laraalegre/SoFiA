#! /usr/bin/env python
from ast import literal_eval
from sofia import error as err
from sofia import __version__ as sofia_version


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
		err.error("Failed to read parameter file: " + str(filename) + "\n" + str(e), fatal=True)
	
	# Extract lines from parameter file
	lines = f.readlines()
	f.close()
	
	# Remove leading/trailing whitespace and empty lines
	lines = [line.strip() for line in lines]
	
	# Check for version number
	for line in lines:
		if "# Creator: SoFiA" in line:
			par_file_version = line[17:]
			if par_file_version != sofia_version:
				err.warning(
					"The parameter file was created with a different version of SoFiA\n"
					"(" + str(par_file_version) + ") than the one you are currently using (" + str(sofia_version) + ").\n"
					"Some settings defined in the parameter file may not be recognised\n"
					"by SoFiA, which could lead to unexpected results.", frame=True)
	
	# Remove comments
	lines = [line for line in lines if len(line) > 0 and line[0] != "#"]
	
	# Some additional setup
	datatypes = allowedDataTypes()
	tasks = {}
	pedantic = True     # Controls whether to exit on unrecognised or redefined parameters
	pedantic_counter = 0
	
	# Loop over all lines:
	for line in lines:
		# Extract parameter name and value
		try:
			parameter, value = tuple(line.split("=", 1))
			parameter = parameter.strip()
			value = value.split("#")[0].strip()
			module, parname = tuple(parameter.split(".", 1))
		except:
			err.error("Failed to read parameter: " + str(line) + "\nExpected format: module.parameter = value", fatal=True)
		
		# Ensure that module and parameter names are not empty
		if len(module) < 1 or len(parname) < 1:
			err.error("Failed to read parameter: " + str(line) + "\nExpected format: module.parameter = value", fatal=True)
		
		subtasks = tasks
		if module not in subtasks: subtasks[module] = {}
		subtasks = subtasks[module]
		
		if parname in subtasks:
			err.warning("Multiple definitions of parameter " + str(parameter) + " encountered.\nIgnoring all additional definitions.")
			pedantic_counter += 1
			continue
		
		if parameter in datatypes:
			try:
				if datatypes[parameter]   == "bool":  subtasks[parname] = str2bool(value)
				elif datatypes[parameter] == "float": subtasks[parname] = float(value)
				elif datatypes[parameter] == "int":   subtasks[parname] = int(value)
				elif datatypes[parameter] == "array": subtasks[parname] = literal_eval(value)
				else: subtasks[parname] = str(value)
			except:
				err.error("Failed to parse parameter value:\n" + str(line) + "\nExpected data type: " + str(datatypes[parameter]), fatal=True)
			if parameter == "pipeline.pedantic": pedantic = subtasks[parname]  # Update 'pedantic' setting
		else:
			err.warning("Ignoring unknown parameter: " + str(parameter) + " = " + str(value))
			pedantic_counter += 1
			continue
	
	if pedantic and pedantic_counter:
		err.error("Multiply-defined or unrecognised parameter(s) encountered.\nPlease check your parameter file or set pipeline.pedantic = false\nto ignore unknown or already defined parameters.", fatal=True)
	
	return tasks


# ----------------------------------------------------------
# FUNCTION: Return list of allowed parameter names and types
# ----------------------------------------------------------

def allowedDataTypes():
	# Define data types of individual parameters. This is used to
	# convert the input values into the correct data type.
	# Ensure that all new parameters get added to this list!
	# Parameters not listed here will be ignored by the parser!
	return {"pipeline.pedantic": "bool", \
	        "steps.doSubcube": "bool", \
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
