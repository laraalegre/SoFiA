#! /usr/bin/env python

def GaussianNoise(F,N0,s0):
    return N0*np.exp(-F**2/2/s0**2)

def SizeFilter(mskt,sfx,sfy,sfz,sbx,sby,sbz,zt,sizeFilter,edgeMode='constant',verbose=0):
	mskt=nd.filters.gaussian_filter(mskt,[0,mt.sqrt(sfy**2+sby**2)/2.355,mt.sqrt(sfx**2+sbx**2)/2.355],mode=edgeMode)
	if zt=='b': mskt=nd.filters.uniform_filter1d(mskt,max(sbz,sfz),axis=0,mode=edgeMode)
	elif zt=='g': mskt=nd.filters.gaussian_filter1d(mskt,max(sbz,sfz/2.355),axis=0,mode=edgeMode)
	mskt[mskt< sizeFilter]=0
	mskt[mskt>=sizeFilter]=1
	return mskt

#def MaskedCube(incube,msk,replace_value):
#	maskedcube=np.copy(incube)
#	maskedcube[msk]=np.sign(incube[msk])*np.minimum(abs(incube[msk]),replace_value)
#	# this only decreases the absolute value of voxels already in the mask, or leaves it unchanged
#	# if already lower than replace_value; the sign is unchanged
#	return maskedcube

def SortKernels(kernels):
	# Sorting kernels
	uniquesky=[]
	velsmooth=[]
	velfshape=[]
	for jj in np.array(kernels):
		if list(jj[:2].astype(float)) not in uniquesky: uniquesky.append(list(jj[:2].astype(float)))
	uniquesky=[kk[1] for kk in sorted([(float(jj[0]),jj) for jj in uniquesky])]

	for jj in uniquesky:
		velsmooth.append([])
		velfshape.append([])
		for ii in np.array(kernels):
			if list(ii[:2].astype(float))==jj:
				velsmooth[-1].append(int(ii[2]))
				velfshape[-1].append(ii[3])
	return uniquesky,velsmooth,velfshape

#def SCfinder(cube,header,t0,kernels=[[0,0,0,'b'],],threshold=3.5,sizeFilter=0,maskScaleXY=2.,maskScaleZ=2.,kernelUnit='pixel',edgeMode='constant',rmsMode='negative',verbose=0):
#	# Create binary mask array
#	msk=np.zeros(cube.shape,'bool')
#	found_nan=np.isnan(cube).sum()
#
#	# Measure noise in original cube
#	rms=GetRMS(cube,rmsMode=rmsMode,zoomx=1,zoomy=1,zoomz=1,verbose=verbose)
#
#	# Sort kernels
#	uniquesky,velsmooth,velfshape=SortKernels(kernels)
#
#	# Loop over all xy kernels
#	for jj in range(len(uniquesky)):
#		[kx,ky]=uniquesky[jj]
#		if kernelUnit=='world' or kernelUnit=='w':
#			kx=abs(float(kx)/header['cdelt1']/3600)
#			ky=abs(float(ky)/header['cdelt2']/3600)
#		if verbose: 
#			print "\n--- %.3f seconds since start"%(time()-t0)
#			print '    Filter %2.0f %2.0f %2.0f %s ...'%(kx,ky,0,'-')
#		sys.stdout.flush()
#
#		mskxy=np.zeros(cube.shape,'bool')
#
#		# Get previous xy kernel
#		if jj: [kxold,kyold]=uniquesky[jj-1]
#		else: [kxold,kyold]=[0,0]
#
#		# Gaussian angular smoothing of *clipped* cube if needed
#		if kx+ky:
#			# smooth starting from the latest xy cube (clipped to the detection threshold times maskScaleXY)
#			cubexy=nd.filters.gaussian_filter(np.clip(cubexy,-maskScaleXY*threshold*rmsxy,maskScaleXY*threshold*rmsxy),[0,mt.sqrt(ky**2-kyold**2)/2.355,mt.sqrt(kx**2-kxold**2)/2.355],mode=edgeMode)
#			if found_nan: cubexy[np.isnan(cube)]=np.nan
#			rmsxy=GetRMS(cubexy,rmsMode=rmsMode,zoomx=1,zoomy=1,zoomz=1,verbose=verbose)
#			if found_nan: cubexy=np.nan_to_num(cubexy)
#		elif found_nan: cubexy,rmsxy=np.nan_to_num(cube),rms
#		else: cubexy,rmsxy=np.copy(cube),rms
#
#		# Loop over all z kernels
#		for ii in range(len(velsmooth[jj])):
#			kz=velsmooth[jj][ii]
#			kt=velfshape[jj][ii]
#			if kernelUnit=='world' or kernelUnit=='w': kz=abs(float(kz)/header['cdelt3']*1000)
#
#			# Velocity smoothing of *clipped* cube if needed
#			if kz:
#				if verbose: 
#					print "\n--- %.3f seconds since start"%(time()-t0)
#					print '    Filter %2.0f %2.0f %2.0f %s ...'%(kx,ky,kz,kt)
#				sys.stdout.flush()
#				if not ii:
#					rmsxyz=rmsxy/mt.sqrt(kz)
#					print '!!!'
#				# smooth starting from the latest xy cube
#				# before smoothing voxels already detected at the current angular resolution are brought down to +/-maskScaleZ*rmsxyz*threshold
#				if kt=='b': cubexyz=nd.filters.uniform_filter1d(MaskedCube(cubexy,mskxy,maskScaleZ*rmsxyz*threshold),kz,axis=0,mode=edgeMode)
#				elif kt=='g': cubexyz=nd.filters.gaussian_filter1d(MaskedCube(cubexy,mskxy,maskScaleZ*rmsxyz*threshold),kz/2.355,axis=0,mode=edgeMode)
#				if found_nan: cubexyz[np.isnan(cube)]=np.nan
#				rmsxyz=GetRMS(cubexyz,rmsMode=rmsMode,zoomx=1,zoomy=1,zoomz=1,verbose=verbose)
#				if found_nan: cubexyz=np.nan_to_num(cubexyz)
#			else: cubexyz,rmsxyz=cubexy,rmsxy
#
#			# Add detected voxels to mask
#			if sizeFilter:
#				if verbose: print '      Adding detected voxels to xy mask with size filtering'
#				sys.stdout.flush()
#				# Get beam FWHM in pixels
#				if 'BMAJ' in header.keys(): bmaj=header['bmaj']/abs(header['cdelt2']) # assumed to be along y axis
#				else: bmaj=0
#				if 'BMIN' in header.keys(): bmin=header['bmin']/abs(header['cdelt1']) # assumed to be along x axis
#				else: bmin=bmaj
#				mskxy=mskxy+SizeFilter(((cubexyz>=threshold*rmsxyz)+(cubexyz<=-threshold*rmsxyz)).astype('float32'),kx,ky,kz,bmin,bmaj,3,kt,sizeFilter,edgeMode=edgeMode,verbose=0).astype('bool')
#			else:
#				if verbose: print '      Adding detected voxels directly to xy mask ...'
#				sys.stdout.flush()
#				mskxy=mskxy+(cubexyz>=threshold*rmsxyz)+(cubexyz<=-threshold*rmsxyz)
#
#		if verbose: print '    Adding xy mask to final mask ...'
#		sys.stdout.flush()
#		msk=msk+mskxy
#
#	msk[np.isnan(cube)]=False
#	return msk

def SCfinder_mem(cube,header,t0,kernels=[[0,0,0,'b'],],threshold=3.5,sizeFilter=0,maskScaleXY=2.,maskScaleZ=2.,kernelUnit='pixel',edgeMode='constant',rmsMode='negative',verbose=0):
    # Create binary mask array
    msk=np.zeros(cube.shape,'bool')
    found_nan=np.isnan(cube).sum()
    # Set dn x dn x dn box boundaries where to measure noise (value of dn set in next line)
    dn=100
    n0=max(0,int(float(cube.shape[0]-dn)/2))
    n1=max(0,int(float(cube.shape[1]-dn)/2))
    n2=max(0,int(float(cube.shape[2]-dn)/2))
    # Measure noise in original (sub-) cube
    rms=GetRMS(cube[n0:cube.shape[0]-n0,n1:cube.shape[1]-n1,n2:cube.shape[2]-n2],rmsMode=rmsMode,zoomx=1,zoomy=1,zoomz=1,verbose=verbose)
    #rms_sample=GetRMS(cube,rmsMode=rmsMode,zoomx=10,zoomy=10,zoomz=10,verbose=verbose)
    # Loop over all kernels
    for jj in kernels:
    	[kx,ky,kz,kt]=jj
        if verbose:
        	print "\n--- %.3f seconds since start"%(time()-t0)
	    	print '    Filter %s %s %s %s ...'%(kx,ky,kz,kt)
        if kernelUnit=='world' or kernelUnit=='w':
        	if verbose: print '    Converting filter size to pixels ...'
        	kx=abs(float(kx)/header['cdelt1'])
        	ky=abs(float(ky)/header['cdelt2'])
        	kz=abs(float(kz)/header['cdelt3'])
        if kt=='b':
        	if kz!=int(mt.ceil(kz)) and verbose: print '    WARNING: Rounding width of boxcar z kernel to next integer'
        	kz=int(mt.ceil(kz))
        sys.stdout.flush()
        smoothedcube=cube*1.
        if found_nan: smoothedcube=np.nan_to_num(smoothedcube)
        smoothedcube[(smoothedcube>0)*msk]=+maskScaleXY*rms
        smoothedcube[(smoothedcube<0)*msk]=-maskScaleXY*rms
        if kx+ky: smoothedcube=nd.filters.gaussian_filter(smoothedcube,[0,ky/2.355,kx/2.355],mode=edgeMode)
        if kz:
            if kt=='b': smoothedcube=nd.filters.uniform_filter1d(smoothedcube,kz,axis=0,mode=edgeMode)
            elif kt=='g': smoothedcube=nd.filters.gaussian_filter1d(smoothedcube,kz/2.355,axis=0,mode=edgeMode)
        if found_nan: smoothedcube[np.isnan(cube)]=np.nan
        #smoothedrms=GetRMS(smoothedcube,rmsMode=rmsMode,zoomx=10,zoomy=10,zoomz=10,verbose=verbose)/rms_sample*rms
        smoothedrms=GetRMS(smoothedcube[n0:cube.shape[0]-n0,n1:cube.shape[1]-n1,n2:cube.shape[2]-n2],rmsMode=rmsMode,zoomx=1,zoomy=1,zoomz=1,verbose=verbose)
        if found_nan: smoothedcube=np.nan_to_num(smoothedcube)
        msk=msk+(smoothedcube>=threshold*smoothedrms)+(smoothedcube<=-threshold*smoothedrms)
        del(smoothedcube)
    return msk

import astropy.io.fits as pf
import numpy as np
import math as mt
import scipy
from scipy import optimize
import scipy.ndimage as nd
from sys import argv
import string
from os import path
import sys
from functions import *
from time import time
