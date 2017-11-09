#! /usr/bin/env python
import re
import sys
import traceback
import ast


def str2bool(s):
	if s in ['True', 'true', 'TRUE', 'Yes', 'yes', 'YES']:
		return True
	return False


def readPipelineOptions(filename = "pipeline.options"):
	try:
		f = open(filename, 'r')
	except IOError as e:
		sys.stderr.write("ERROR: Failed to read parameter file.\n")
		sys.stderr.write(e + '\n')
		sys.exit(1);
	
	lines = f.readlines()
	f.close()
	
	# Remove leading/trailing whitespace and empty lines:
	lines = [line.strip() for line in lines]
	lines = [line for line in lines if len(line) > 0]
	lines2 = []
	
	# The following piece of code is meant to allow line breaks in parameter settings 
	# indicated by a double backwards slash (\\):
	for l in lines:
		if len(lines2) > 0:
			if lines2[-1][-1] == "\\":
				lines2[-1] = lines2[-1][:-1] + l
				continue
		lines2.append(l)
	lines = lines2
	tasks = {}
	
	# Define data types of individual parameters. This is used to 
	# convert the input values into the correct data type.
	# Ensure that all new parameters get added to this list!
	datatypes = {"steps.doSubcube": "bool", \
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
	             "flag.regions": "array", \
	             "optical.sourceCatalogue": "string", \
	             "optical.spatSize": "float", \
	             "optical.specSize": "float", \
	             "optical.storeMultiCat": "bool", \
	             "smooth.kernel": "string", \
	             "smooth.edgeMode": "string", \
	             "smooth.kernelX": "float", \
	             "smooth.kernelY": "float", \
	             "smooth.kernelZ": "float", \
	             "scaleNoise.statistic": "string", \
	             "scaleNoise.edgeX": "int", \
	             "scaleNoise.edgeY": "int", \
	             "scaleNoise.edgeZ": "int", \
	             "scaleNoise.scaleX": "bool", \
	             "scaleNoise.scaleY": "bool", \
	             "scaleNoise.scaleZ": "bool", \
	             "SCfind.threshold": "float", \
	             "SCfind.sizeFilter": "float", \
	             "SCfind.maskScaleXY": "float", \
	             "SCfind.maskScaleZ": "float", \
	             "SCfind.edgeMode": "string", \
	             "SCfind.rmsMode": "string", \
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
	             "writeCat.parameters": "array"
	}
	
	# Loop through all lines:
	for linenr, line in enumerate(lines):
		if line[0] == "#":
			continue
		
		parameter, value = tuple(line.split("=", 1))
		if len(parameter) < 1:
			sys.stderr.write("FATAL ERROR: Parameter name missing in line %i of parameter file %s:\n%s\n"%(linenr+1,filename,line))
			sys.exit(1)
		
		subtasks = tasks
		
		while True:
			module = parameter.split(".")[0]
			module = module.strip()
			
			if len(module) < 1:
				sys.stderr.write("FATAL ERROR: (Sub)key name too short in line %i of parameter file %s:\n%s\n"%(linenr+1,filename,line))
				sys.exit(1)
			
			parameter = str(".").join(parameter.split(".")[1:])
			parameter = parameter.strip()
			
			if not subtasks.has_key(module):
				subtasks[module] = {}
			subtasks = subtasks[module]
			
			if parameter.count(".") == 0:
				if subtasks.has_key(parameter):
					sys.stderr.write("WARNING: Parameter already present in line %i of parameter file %s:\n%s\n"%(linenr+1,filename,line))
					sys.stderr.write('         Will ignore repeated parameter.\n')
					break
				try:
					value = value.split('#')[0].strip()
					searchKey = module + "." + parameter
					
					if searchKey in datatypes:
						if datatypes[searchKey]   == "bool":  subtasks[parameter] = str2bool(value)
						elif datatypes[searchKey] == "float": subtasks[parameter] = float(value)
						elif datatypes[searchKey] == "int":   subtasks[parameter] = int(value)
						elif datatypes[searchKey] == "array": subtasks[parameter] = ast.literal_eval(value)
						else: subtasks[parameter] = str(value)
					else:
						sys.stderr.write("WARNING: No data type defined for parameter %s.\n"%(searchKey))
						sys.stderr.write("         Guessing type based on value.\n")
						subtasks[parameter] = ast.literal_eval(value)
				except:
					sys.stderr.write('FATAL ERROR: Could not parse option in line %i of user parameter file %s:\n'%(linenr+1,filename))
					sys.stderr.write('             %s\n'%(line))
					sys.stderr.write('             Expecting data type: %s\n'%(datatypes[searchKey]))
					sys.exit(1)
				break
	
	return tasks
