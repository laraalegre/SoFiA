#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np


def make_ascii_from_array(objects, cathead, catunits, catfmt, store_pars, outname, compress, flagOverwrite, flagUncertainties):
	if compress: outname += ".gz";
	print "Writing ASCII catalogue: ", outname;
	
	# Do we need to write all parameters?
	parList = list(cathead) if (store_pars == ["*"]) else store_pars;
	
	# Do we want statistical uncertainties?
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
	for j in catfmt: lenCathead.append(int(j.split("%")[1].split("e")[0].split("f")[0].split("i")[0].split("d")[0].split(".")[0].split("s")[0]) + 1);
	
	# Create header
	headerName = "";
	headerUnit = "";
	headerCol  = "";
	outFormat  = "";
	colCount   =  0;
	
	for par in parList:
		if par in cathead:
			index = list(cathead).index(par);
			headerName += cathead[index].rjust(lenCathead[index]);
			headerUnit += catunits[index].rjust(lenCathead[index]);
			headerCol  += ("(%i)" % (colCount + 1)).rjust(lenCathead[index]);
			outFormat  += catfmt[list(cathead).index(par)] + " ";
			colCount += 1;
		else:
			sys.stderr.write("WARNING: Skipping undefined parameter \'" + str(par) + "\'.\n");
	header += headerName[3:] + '\n' + headerUnit[3:] + '\n' + headerCol[3:];
	
	# Create catalogue
	outObjects = [];
	for obj in objects:
		outObjects.append([]);
		for par in parList:
			if par in cathead: outObjects[-1].append(obj[list(cathead).index(par)]);
	
	# Write catalogue with check of overwrite flag
	if not flagOverwrite and os.path.exists(outname):
		sys.stderr.write("ERROR: Output file exists: " + outname + ".\n");
	else:
		np.savetxt(outname, np.array(outObjects, dtype=object), fmt=outFormat, header=header);
