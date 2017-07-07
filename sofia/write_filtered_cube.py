#! /usr/bin/env python
from astropy.io import fits
import os

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


def writeFilteredCube(cube, header, dictionary, filename, compress):
    header.add_history('SoFiA input filtering')
    optionsList = []
    optionsDepth = []
    dictionary = removeOptions(dictionary)
    recursion(dictionary, optionsList, optionsDepth)
    headerList = []
    for i in range(0, len(optionsList)):
    #{
        if len(optionsList[i].split('=')) > 1:
        #{
            tmpString = optionsList[i]
            depthNumber = optionsDepth[i]
            j = i - 1
            while depthNumber > 0:
            #{
                if optionsDepth[i] > optionsDepth[j]:
                #{
                    tmpString = optionsList[j] + '.' + tmpString
                    depthNumber = optionsDepth[j]
                #}
                j -= 1
            #}
            headerList.append(tmpString)
        #}
    #}
    for option in headerList:
    #{
        header.add_history(option)
    #}
    #if cube.max() < 65535: cube = cube.astype('int16')
    hdu = fits.PrimaryHDU(data = cube, header = header)
    hdu.header['datamin'] = cube.min()
    hdu.header['datamax'] = cube.max()
    #hdulist = fits.HDUList([hdu])
    #name = os.path.splitext(filename)[0] + '_mask.fits'
    
    if compress:
    #{
        filename += '.gz'
    #}
    hdu.writeto(filename,output_verify='warn', clobber=True)
