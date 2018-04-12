#! /usr/bin/env python
import os
from sofia import error as err

# ---------------------------------------
# Function to return SoFiA version number
# ---------------------------------------

def getVersion(full = False):
	version = "0.0.0"
	fileVersionPath = os.environ["SOFIA_PIPELINE_PATH"]
	fileVersionPath = fileVersionPath.replace("sofia_pipeline.py", "VERSION")
	
	try:
		with open(fileVersionPath) as fileVersion:
			version = (fileVersion.readline()).strip()
	except:
		err.warning("Failed to read SoFiA version number.")
	
	if full: return "SoFiA " + str(version)
	return str(version)
