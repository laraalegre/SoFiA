#! /usr/bin/env python
import os
import sys

# ---------------------------------------
# Function to return SoFiA version number
# ---------------------------------------

def getVersion(full = False):
	version = "0.0";
	fileVersionPath = os.environ["SOFIA_PIPELINE_PATH"];
	fileVersionPath = fileVersionPath.replace("sofia_pipeline.py", "VERSION");
	
	try:
		with open(fileVersionPath) as fileVersion:
			for line in fileVersion:
			   if line: version = line.strip();
	except:
		sys.stderr.write("WARNING: Failed to read SoFiA version number.\n");
	
	if(full): return "SoFiA " + str(version)
	return str(version)
