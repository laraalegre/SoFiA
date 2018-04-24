#! /usr/bin/env python
import os
from astropy.io import fits
from sofia import error as err
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
	return


def writeMask(cube, header, dictionary, filename, compress, flagOverwrite):
	header.add_history("SoFiA source finding")
	optionsList = []
	optionsDepth = []
	dictionary = removeOptions(dictionary)
	recursion(dictionary,optionsList,optionsDepth)
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
				#end if
				j -= 1
			#end while
			headerList.append(tmpString)
		#end if
	#end for
	for option in headerList: header.add_history(option)
	if cube.max() < 32767: cube=cube.astype("int16")
	
	# add axes required to make the shape of the mask cube equal to the shape of the input datacube
	while header["naxis"] > len(cube.shape): cube.resize(tuple([1,] + list(cube.shape)))
	
	hdu = fits.PrimaryHDU(data=cube, header=header)
	hdu.header["BUNIT"] = "source_ID"
	hdu.header["DATAMIN"] = cube.min()
	hdu.header["DATAMAX"] = cube.max()
	hdu.header["ORIGIN"] = sofia_version_full
	
	name = filename
	if compress: name += ".gz"
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		err.error("Output file exists: " + name + ".", fatal=False)
		return
	else:
		hdu.writeto(name, output_verify="warn", **__astropy_arg_overwrite__)
	
	return
