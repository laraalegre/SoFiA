#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np


def make_ascii_from_array(objects, catHeader, catUnits, catFormat, parList, outName, flagCompress, flagOverwrite, flagUncertainties):
	if flagCompress: outName += ".gz";
	print "Writing ASCII catalogue: ", outName;
	
	# Do we need to write all parameters?
	if parList == ["*"] or not parList: parList = list(catHeader);
	
	# Remove undefined parameters
	parList = [item for item in parList if item in catHeader];
	
	# Remove statistical uncertainties if not requested
	if not flagUncertainties:
		for item in ["err_x", "err_y", "err_z", "err_w20", "err_w50"]:
			while item in parList: parList.remove(item);
	
	# Add SoFiA version number to header
	version = "[unknown]";
	fileVersionPath = os.environ["SOFIA_PIPELINE_PATH"];
	fileVersionPath = fileVersionPath.replace("sofia_pipeline.py", "VERSION");
	try:
		with open(fileVersionPath) as fileVersion:
			for line in fileVersion:
			   if line: version = line.strip();
	except:
		sys.stderr.write("WARNING: Failed to read SoFiA version number.\n");
	header = "SoFiA catalogue (version %s)\n" % version;
	
	# Determine header sizes based on variable-length formatting
	lenCathead = [];
	for j in catFormat: lenCathead.append(int(j.split("%")[1].split("e")[0].split("f")[0].split("i")[0].split("d")[0].split(".")[0].split("s")[0]) + 1);
	
	# Create header
	headerName = "";
	headerUnit = "";
	headerCol  = "";
	outFormat  = "";
	colCount   =  0;
	
	for par in parList:
		index = list(catHeader).index(par);
		headerName += catHeader[index].rjust(lenCathead[index]);
		headerUnit += catUnits[index].rjust(lenCathead[index]);
		headerCol  += ("(%i)" % (colCount + 1)).rjust(lenCathead[index]);
		outFormat  += catFormat[list(catHeader).index(par)] + " ";
		colCount += 1;
	header += headerName[3:] + '\n' + headerUnit[3:] + '\n' + headerCol[3:];
	
	# Create catalogue
	outObjects = [];
	for obj in objects:
		outObjects.append([]);
		for par in parList: outObjects[-1].append(obj[list(catHeader).index(par)]);
	
	# Write catalogue with check of overwrite flag
	if not flagOverwrite and os.path.exists(outName):
		sys.stderr.write("ERROR: Output file exists: " + outName + ".\n");
	else:
		np.savetxt(outName, np.array(outObjects, dtype=object), fmt=outFormat, header=header);
