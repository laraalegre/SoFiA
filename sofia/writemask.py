#! /usr/bin/env python
from astropy.io import fits
import os
import sys

def removeOptions(dictionary):
	modDictionary = dictionary
	for key in modDictionary['steps']:
		if modDictionary['steps'][key] == 'False':
			modDictionary.pop(key.split('do')[1].lower(), None)
	return modDictionary


def recursion(dictionary, optionsList, optionsDepth, counter=0):
	if type(dictionary) == type({}):
		for k in dictionary:
			optionsList.append(str(k))
			optionsDepth.append(counter)
			recursion(dictionary[k], optionsList, optionsDepth, counter=counter+1)
	else:
		optionsList[len(optionsList) - 1] += '=' + str(dictionary)
		counter = 0


def writeMask(cube, header, dictionary, filename, compress, flagOverwrite):
	header.add_history('SoFiA source finding')
	optionsList = []
	optionsDepth = []
	dictionary = removeOptions(dictionary)
	recursion(dictionary,optionsList,optionsDepth)
	headerList = []
	for i in range(0, len(optionsList)):
		if len(optionsList[i].split('=')) > 1:
			tmpString = optionsList[i]
			depthNumber = optionsDepth[i]
			j = i - 1
			while depthNumber > 0:
				if optionsDepth[i] > optionsDepth[j]:
					tmpString = optionsList[j] + '.' + tmpString
					depthNumber = optionsDepth[j]
				j -= 1
			headerList.append(tmpString)
	for option in headerList:
		header.add_history(option)
	if cube.max() < 32767:
		cube=cube.astype('int16')

	# add axes required to make the shape of the mask cube equal to the shape of the input datacube
        while header['naxis']>len(cube.shape): cube.resize(tuple([1,]+list(cube.shape)))

	hdu = fits.PrimaryHDU(data=cube, header=header)
	hdu.header['BUNIT'] = 'source_ID'
	hdu.header['DATAMIN'] = cube.min()
	hdu.header['DATAMAX'] = cube.max()
	#hdulist = fits.HDUList([hdu])
	#name = os.path.splitext(filename)[0] + '_mask.fits'
	name = filename
	if compress:
		name += '.gz'

	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		hdu.writeto(name, output_verify='warn', clobber=True)
