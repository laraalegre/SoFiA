#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
from numpy import *
import numpy as np


def make_sql_from_array(objects, cathead, catunits, catfmt, store_pars, outname, compress, flagOverwrite):
	print "Store the results to SQL file: ", outname;
	
	
	# Recover SoFiA version number:
	version = "[unknown]";
	fileVersionPath = os.environ["SOFIA_PIPELINE_PATH"];
	fileVersionPath = fileVersionPath.replace("sofia_pipeline.py", "VERSION");
	try:
		with open(fileVersionPath) as fileVersion:
			for line in fileVersion:
			   if line: version = line.strip();
	except:
		sys.stderr.write("WARNING: Failed to read SoFiA version number.\n");
	
	
	# No idea why this would be necessary:
	objects = np.array(objects);
	
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(outname):
		sys.stderr.write("ERROR: Output file exists: " + outname + ".\n");
	else:
		title = "-- SoFiA catalogue (version %s)\n\n"%version;
		
		# Open output file:
		try:
			fp = open(outname, "w");
		except:
			sys.stderr.write("ERROR: Failed to write to SQL output file: " + outname + ".\n");
			return;
		
		# Write some header information:
		fp.write(title);
		fp.write("SET SQL_MODE = \"NO_AUTO_VALUE_ON_ZERO\";\n\n");
		
		# Write table structure:
		counter = 0;
		fp.write("CREATE TABLE IF NOT EXISTS `SoFiA-Catalogue` (");
		for i in range(0, len(cathead)):
			if store_pars == ["*"] or cathead[i] in store_pars:
				if counter > 0: fp.write(", ");
				fp.write("`" + cathead[i] + "` ");
				if   "i" in catfmt[i] or "d" in catfmt[i]: fp.write("INT");
				elif "f" in catfmt[i] or "e" in catfmt[i]: fp.write("DOUBLE");
				else: fp.write("VARCHAR(256)");
				fp.write(" NOT NULL");
				counter += 1;
		if store_pars == ["*"] or "id" in store_pars: fp.write(", PRIMARY KEY (`id`), KEY (`id`)");
		fp.write(") DEFAULT CHARSET=utf8 COMMENT=\'SoFiA source catalogue\';\n\n");
		
		# Write data:
		counter = 0;
		fp.write("INSERT INTO `SoFiA-Catalogue` (");
		for i in range(0, len(cathead)):
			if store_pars == ["*"] or cathead[i] in store_pars:
				if counter > 0: fp.write(", ");
				fp.write("`" + cathead[i] + "`");
				counter += 1;
		fp.write(") VALUES\n");
		
		for i in range(0, len(objects)):
			counter = 0;
			fp.write("(");
			for j in range(0, len(objects[i])):
				if store_pars == ["*"] or cathead[j] in store_pars:
					if counter > 0: fp.write(", ");
					if "f" in catfmt[j] or "e" in catfmt[j]: fp.write(str(float(objects[i][j])));
					elif "i" in catfmt[j] or "d" in catfmt[j]: fp.write(str(int(objects[i][j])));
					else: fp.write("\'" + str(objects[i][j]) + "\'");
					counter += 1;
			if(i < len(objects) - 1): fp.write("),\n");
			else: fp.write(");\n");
		
		# Close output file:
		fp.close();
	
	return;
