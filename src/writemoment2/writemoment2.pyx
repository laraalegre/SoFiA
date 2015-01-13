#!/usr/bin/env python
import astropy.io.fits as pyfits
import os
import numpy as np
cimport numpy as np
from libc.math cimport isnan

def removeOptions(dictionary):
  modDictionary = dictionary
  for key in modDictionary['steps']:
    if modDictionary['steps'][key] == 'False':
      modDictionary.pop(key.split('do')[1].lower(), None)
  return modDictionary


def recursion(dictionary,optionsList,optionsDepth,counter=0):
  if type(dictionary)==type({}):
    for k in dictionary:
      optionsList.append(str(k))
      optionsDepth.append(counter)
      recursion(dictionary[k],optionsList,optionsDepth,counter=counter+1)
  else:
    optionsList[len(optionsList)-1]+='='+str(dictionary)
    counter = 0

def writeMoment0(datacube,maskcube,filename,debug,header,compress):
  print 'Writing moment-0' # in units of header['bunit']*km/s
  if 'CELLSCAL' in header:
    if header['CELLSCAL'] == '1/F':
      print 'WARNING: CELLSCAL keyword with value 1/F found'
  datacube = np.array(datacube, dtype=np.single)
  maskcube = np.array(maskcube, dtype=np.single)
  m0 = mom0(datacube,maskcube)  
  op=0
  if 'vopt' in header['ctype3'].lower() or 'vrad' in header['ctype3'].lower() or 'velo' in header['ctype3'].lower() or 'felo' in header['ctype3'].lower():
    if not 'cunit3' in header: dkms=abs(header['cdelt3'])/1e+3 # assuming m/s
    elif header['cunit3'].lower()=='km/s': dkms=abs(header['cdelt3'])
    bunitExt = '.km/s'
  elif 'freq' in header['ctype3'].lower():
    #if not 'cunit3' in header or header['cunit3'].lower()=='hz': dkms=abs(header['cdelt3'])/1.42040575177e+9*2.99792458e+5 # assuming Hz
    #elif header['cunit3'].lower()=='khz': dkms=abs(header['cdelt3'])/1.42040575177e+6*2.99792458e+5
    if not 'cunit3' in header or ('cunit3' in header and header['cunit3'].lower()=='hz'):
      bunitExt = '.Hz'
    else: 
      bunitExt = '.'+header['cunit3']
    dkms=1. # no scaling, avoids crashing
  hdu = pyfits.PrimaryHDU(data=m0*dkms,header=header)
  hdu.header['bunit']+=bunitExt
  hdu.header['datamin']=(m0*dkms).min()
  hdu.header['datamax']=(m0*dkms).max()
  del(hdu.header['crpix3'])
  del(hdu.header['crval3'])
  del(hdu.header['cdelt3'])
  del(hdu.header['ctype3'])
  if debug: hdu.writeto('%s_mom0.debug.fits'%filename,output_verify='warn',clobber=True)
  else: 
    name = '%s_mom0.fits'%filename
    if compress:
      name += '.gz'
    hdu.writeto(name,output_verify='warn',clobber=True)
  return m0


def writeMoment1(datacube,maskcube,filename,debug,header,m0,compress):
  print 'Writing moment-1'
  datacube = np.array(datacube, dtype=np.single)
  maskcube = np.array(maskcube, dtype=np.single)
  m1 = mom1(datacube,maskcube,m0,header['crpix3'],header['crval3'],header['cdelt3'])
  # convert it to km/s (using radio velocity definition to go from Hz to km/s)
  if 'vopt' in header['ctype3'].lower() or 'vrad' in header['ctype3'].lower() or 'velo' in header['ctype3'].lower() or 'felo' in header['ctype3'].lower():
    if not 'cunit3' in header:
      m1/=1e+3 # assuming m/s
      bunitExt = 'km/s'
    elif header['cunit3'].lower()=='km/s': 
      bunitExt = 'km/s'
    else:
      bunitExt = header['cunit3']
  elif 'freq' in header['ctype3'].lower():
    #if not 'cunit3' in header or header['cunit3'].lower()=='hz': m1*=2.99792458e+5/1.42040575177e+9 # assuming Hz
    #elif header['cunit3'].lower()=='khz': m1*=2.99792458e+5/1.42040575177e+6
    if not 'cunit3' in header or ('cunit3' in header and header['cunit3'].lower()=='hz'):
      bunitExt = 'Hz'
    else: 
      bunitExt = header['cunit3']
    dkms=1. # no scaling, avoids crashing
  # calculate moment 1
  hdu = pyfits.PrimaryHDU(data=m1,header=header)
  hdu.header['bunit']=bunitExt
  hdu.header['datamin']=np.nanmin(m1)
  hdu.header['datamax']=np.nanmax(m1)
  del(hdu.header['crpix3'])
  del(hdu.header['crval3'])
  del(hdu.header['cdelt3'])
  del(hdu.header['ctype3'])
  if debug: hdu.writeto('%s_mom1.debug.fits'%filename,output_verify='warn',clobber=True)
  else: 
    name = '%s_mom1.fits'%filename
    if compress:
      name += '.gz'
    hdu.writeto(name,output_verify='warn',clobber=True)


def mom0(cube1, cube2):
  
  cdef:
    int i, j, k
    double sum
    double [:,:] mom0 = np.zeros((cube1.shape[1],cube1.shape[2]))
    float[:,:,:] cube = cube1
    float[:,:,:] mask = cube2
    
  for j in range(cube.shape[1]):
    for k in range(cube.shape[2]):
      sum = 0
      for i in range(cube.shape[0]):
        if not isnan(cube[i,j,k]) and mask[i,j,k]!=0:
          sum += cube[i,j,k]
      mom0[j,k] = sum

  return np.array(mom0)

def mom1(cube1, cube2, cube3, int cpx, float cval, float cdelt):

  cdef:
    int i, j, k
    double sum
    double [:,:] mom1 = np.zeros((cube1.shape[1],cube1.shape[2]))
    float[:,:,:] cube = cube1
    float[:,:,:] mask = cube2
    double[:,:] mom0 = cube3
  
  for j in range(cube.shape[1]):
    for k in range(cube.shape[2]):
      sum = 0
      for i in range(cube.shape[0]):
        if not isnan(cube[i,j,k]) and mask[i,j,k]!=0:
          sum += ((i+1-cpx)*cdelt+cval)*cube[i,j,k]
      if mom0[j,k] != 0 and not isnan(mom0[j,k]):
        mom1[j,k] = sum/mom0[j,k]
      else:
        mom1[j,k] = np.nan

  return np.array(mom1)
