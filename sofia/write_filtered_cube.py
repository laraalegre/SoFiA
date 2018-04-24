#! /usr/bin/env python
from astropy.io import fits
from numpy import nanmin, nanmax
from sofia import __version_full__ as sofia_version_full
from sofia import __astropy_arg_overwrite__

def removeOptions(dictionary):
	modDictionary = dictionary
	for key in modDictionary["steps"]:
		if modDictionary["steps"][key] == "False":
			modDictionary.pop(key.split("do")[1].lower(), None)
	return modDictionary


def recursion(dictionary, optionsList, optionsDepth, counter=0):
	if type(dictionary) == type({}):
		for k in dictionary:
			optionsList.append(str(k))
			optionsDepth.append(counter)
			recursion(dictionary[k], optionsList, optionsDepth, counter=counter+1)
	else:
		optionsList[len(optionsList) - 1] += "=" + str(dictionary)
		counter = 0


def writeFilteredCube(cube, header, dictionary, filename, compress):
	optionsList = []
	optionsDepth = []
	dictionary = removeOptions(dictionary)
	recursion(dictionary, optionsList, optionsDepth)
	headerList = []
	
	for i in range(0, len(optionsList)):
		if len(optionsList[i].split("=")) > 1:
			tmpString = optionsList[i]
			depthNumber = optionsDepth[i]
			j = i - 1
			while depthNumber > 0:
				if optionsDepth[i] > optionsDepth[j]:
					tmpString = optionsList[j] + "." + tmpString
					depthNumber = optionsDepth[j]
				j -= 1
			headerList.append(tmpString)
	for option in headerList: header.add_history(option)
	
	hdu = fits.PrimaryHDU(data = cube, header = header)
	hdu.header["DATAMIN"] = nanmin(cube)
	hdu.header["DATAMAX"] = nanmax(cube)
	hdu.header["ORIGIN"] = sofia_version_full
	
	if compress: filename += ".gz"
	hdu.writeto(filename,output_verify="warn", **__astropy_arg_overwrite__)
