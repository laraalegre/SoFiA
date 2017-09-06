#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys;
import os;
import numpy as np;
from xml.etree.ElementTree import Element, SubElement, tostring;
from xml.dom import minidom;
from gzip import open as gzopen;


# --------------------------------
# Function to auto-indent XML code
# --------------------------------

def prettify(elem):
	# Indent is set to "" here to save disk space; default would be "\t".
	rough_string = tostring(elem, "utf-8");
	reparsed = minidom.parseString(rough_string);
	return reparsed.toprettyxml(indent="");


# -----------------------------------
# Function to create SQL header entry
# -----------------------------------

def sqlHeaderItem(item):
	return "`" + item + "`";


# ---------------------------------
# Function to create SQL data entry
# ---------------------------------

def sqlDataItem(item, dataFormat):
	if "f" in dataFormat or "e" in dataFormat: return str(float(item));
	if "i" in dataFormat or "d" in dataFormat: return str(int(item));
	return "\'" + str(item) + "\'";


# -----------------------------------------
# Function to create SQL data format string
# -----------------------------------------

def sqlFormat(item):
	if "f" in item or "e" in item: return " double NOT NULL";
	if "i" in item or "d" in item: return " int NOT NULL";
	return " varchar(256) NOT NULL";


# ----------------------------------
# Function to write source catalogue
# ----------------------------------

def write_catalog_from_array(mode, objects, catHeader, catUnits, catFormat, parList, outName, flagCompress, flagOverwrite, flagUncertainties):
	# Check output format and compression
	availableModes = ["ASCII", "XML", "SQL"];
	if mode not in availableModes:
		sys.stderr.write("WARNING: Unknown catalogue format: %s. Defaulting to ASCII.\n" % mode);
		mode = "ASCII";
	modeIndex = availableModes.index(mode);
	
	if flagCompress: outName += ".gz";
	print "Writing " + availableModes[modeIndex] + " catalogue: " + outName;
	
	# Exit if file exists and overwrite flag is set to false
	if not flagOverwrite and os.path.exists(outName):
		sys.stderr.write("ERROR: Output file exists: " + outName + ".\n");
		return;
	
	
	# Do we need to write all parameters?
	if parList == ["*"] or not parList: parList = list(catHeader);
	
	# Remove undefined parameters
	parList = [item for item in parList if item in catHeader];
	
	# Remove statistical uncertainties if not requested
	if not flagUncertainties:
		for item in ["err_x", "err_y", "err_z", "err_w20", "err_w50"]:
			while item in parList: parList.remove(item);
	
	# Check whether there is anything left
	if not len(parList):
		sys.stderr.write("ERROR: No valid output parameters selected. No output catalogue written.\n");
		return;
	
	
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
	
	
	# Create and write catalogue in requested format
	# -------------------------------------------------------------------------
	if mode == "XML":
		# Define basic XML header information
		votable          = Element("VOTABLE");
		resource         = SubElement(votable, "RESOURCE", name="SoFiA catalogue (version %s)" % version);
		description      = SubElement(resource, "DESCRIPTION");
		description.text = "Source catalogue from the Source Finding Application (SoFiA) version %s" % version;
		coosys           = SubElement(resource, "COOSYS", ID="J2000");
		table            = SubElement(resource, "TABLE", ID="sofia_cat", name="sofia_cat");
		
		# Load list of parameters and unified content descriptors (UCDs)
		ucdList = {};
		fileUcdPath = os.environ["SOFIA_PIPELINE_PATH"];
		fileUcdPath = fileUcdPath.replace("sofia_pipeline.py", "SoFiA_source_parameters.dat");
		
		try:
			with open(fileUcdPath) as fileUcd:
				for line in fileUcd:
					(key, value) = line.split();
					ucdList[key] = value;
		except:
			sys.stderr.write("WARNING: Failed to read UCD file.\n");
		
		# Create parameter fields
		for par in parList:
			ucdEntity = ucdList[par] if par in ucdList else "";
			index = list(catHeader).index(par);
			if catFormat[index] == "%30s":
				field = SubElement(table, "FIELD", name=par, ucd=ucdEntity, datatype="char", arraysize="30", unit=catUnits[index]);
			else:
				field = SubElement(table, "FIELD", name=par, ucd=ucdEntity, datatype="float", unit=catUnits[index]);
		
		# Create data table entries
		data = SubElement(table, "DATA");
		tabledata = SubElement(data, "TABLEDATA");

		for obj in objects:
			tr = SubElement(tabledata, "TR");
			for par in parList:
				td = SubElement(tr, "TD");
				index = list(catHeader).index(par);
				td.text = (catFormat[index] % obj[index]).strip();
		
		# Write XML catalogue:
		try:
			f1 = gzopen(outName, "wb") if flagCompress else open(outName, "w");
		except:
			sys.stderr.write("ERROR: Failed to write to XML output file: " + outName + ".\n");
			return;
		f1.write(prettify(votable));
		f1.close;
	
	# -----------------------------------------------------------------End-XML-
	
	elif mode == "SQL":
		# Record if there is an ID column in the catalogue
		# (if no ID is present, we will later create one for use as primary key)
		noID = "id" not in parList;
		
		# Write some header information:
		content = "-- SoFiA catalogue (version %s)\n\nSET SQL_MODE = \"NO_AUTO_VALUE_ON_ZERO\";\n\n" % version;
		
		# Construct and write table structure:
		flagProgress = False;
		content += "CREATE TABLE IF NOT EXISTS `SoFiA-Catalogue` (\n";
		if noID: content += "  `id` INT NOT NULL,\n";
		for par in parList:
			index = list(catHeader).index(par);
			if flagProgress: content += ",\n";
			content += "  " + sqlHeaderItem(par) + sqlFormat(catFormat[index]);
			flagProgress = True;
		content += ",\n  PRIMARY KEY (`id`),\n  KEY (`id`)\n) DEFAULT CHARSET=utf8 COMMENT=\'SoFiA source catalogue\';\n\n";
		
		# Insert data:
		flagProgress = False;
		content += "INSERT INTO `SoFiA-Catalogue` (";
		if noID: content += "`id`, ";
		for par in parList:
			if flagProgress: content += ", ";
			content += sqlHeaderItem(par);
			flagProgress = True;
		content += ") VALUES\n";
		
		source_count = 0;
		for obj in objects:
			flagProgress = False;
			source_count += 1;
			content += "(";
			if noID: content += str(source_count) + ", ";
			
			for par in parList:
				index = list(catHeader).index(par);
				if flagProgress: content += ", ";
				content += sqlDataItem(obj[index], catFormat[index]);
				flagProgress = True;
			
			if(source_count < len(objects)): content += "),\n";
			else: content += ");\n";
		
		# Write catalogue
		try:
			fp = gzopen(outName, "wb") if flagCompress else open(outName, "w");
		except:
			sys.stderr.write("ERROR: Failed to write to SQL output file: " + outName + ".\n");
			return;
		fp.write(content);
		fp.close();
	
	# -----------------------------------------------------------------End-SQL-
	
	else: # mode == "ASCII" by default
		# Determine header sizes based on variable-length formatting
		lenCathead = [];
		for j in catFormat: lenCathead.append(int(j.split("%")[1].split("e")[0].split("f")[0].split("i")[0].split("d")[0].split(".")[0].split("s")[0]) + 1);
		
		# Create header
		headerName = "";
		headerUnit = "";
		headerCol  = "";
		outFormat  = "";
		colCount   =  0;
		header     = "SoFiA catalogue (version %s)\n" % version;
		
		for par in parList:
			index = list(catHeader).index(par);
			headerName += catHeader[index].rjust(lenCathead[index]);
			headerUnit += catUnits[index].rjust(lenCathead[index]);
			headerCol  += ("(%i)" % (colCount + 1)).rjust(lenCathead[index]);
			outFormat  += catFormat[index] + " ";
			colCount += 1;
		header += headerName[3:] + '\n' + headerUnit[3:] + '\n' + headerCol[3:];
		
		# Create catalogue
		outObjects = [];
		for obj in objects:
			outObjects.append([]);
			for par in parList: outObjects[-1].append(obj[list(catHeader).index(par)]);
		
		# Write ASCII catalogue
		try:
			np.savetxt(outName, np.array(outObjects, dtype=object), fmt=outFormat, header=header);
		
		except:
			sys.stderr.write("ERROR: Failed to write to ASCII catalogue: " + outName + ".\n");
			return;
	
	# ---------------------------------------------------------------End-ASCII-
	
	return;
