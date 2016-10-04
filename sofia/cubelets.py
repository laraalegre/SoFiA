#! /usr/bin/env python
import numpy as np
import astropy.io.fits as pyfits
import os
import sys
from sofia import writemoment
import math
from scipy.ndimage import map_coordinates


def regridMaskedChannels(datacube,maskcube,header):
  
  maskcubeFlt = maskcube.astype('float')
  maskcubeFlt[maskcubeFlt>1] = 1
  
  print 'Regridding...'
  sys.stdout.flush()
  from scipy import interpolate
  z=(np.arange(1.,header['naxis3']+1)-header['crpix3'])*header['cdelt3']+header['crval3']
  if header['ctype3']=='VELO-HEL':
    pixscale=(1-header['crval3']/2.99792458e+8)/(1-z/2.99792458e+8)
  else:
    sys.stderr.write("WARNING: Cannot convert axis3 coordinates to frequency. Will ignore the effect of CELLSCAL = 1/F.\n")
    pixscale=np.ones((header['naxis3']))
  x0,y0=header['crpix1']-1,header['crpix2']-1
  xs=np.arange(datacube.shape[2],dtype=float)-x0
  ys=np.arange(datacube.shape[1],dtype=float)-y0
  for zz in range(datacube.shape[0]):
    regrid_channel=interpolate.RectBivariateSpline(ys*pixscale[zz],xs*pixscale[zz],datacube[zz])
    datacube[zz]=regrid_channel(ys,xs)
    regrid_channel_mask=interpolate.RectBivariateSpline(ys*pixscale[zz],xs*pixscale[zz],maskcubeFlt[zz])
    maskcubeFlt[zz] = regrid_channel_mask(ys,xs)
  
  maskMin = maskcubeFlt.min()
  datacube[abs(maskcubeFlt)<=abs(maskMin)] = np.nan
  del maskcubeFlt
  
  return datacube

def writeSubcube(cube, header, mask, objects, cathead, outroot, compress, flagOverwrite):
    
    # strip path variable to get the file name and the directory separately
    splitroot = outroot.split('/')
    cubename  = splitroot[-1]
    if len(splitroot)>1:
      outputDir = '/'.join(splitroot[:-1])+'/objects/'
    else:
      outputDir = './objects/'
      
    # check if output directory exists and create it if not
    if os.path.exists(outputDir)==False:
      os.system('mkdir '+outputDir)

    
    # read all important information (central pixels & values, increments) from the header
    dX    = header['CDELT1']
    dY    = header['CDELT2']
    dZ    = header['CDELT3']
    cValX = header['CRVAL1']
    cValY = header['CRVAL2']
    cValZ = header['CRVAL3']
    cPixX = header['CRPIX1']-1
    cPixY = header['CRPIX2']-1
    cPixZ = header['CRPIX3']-1
    #specTypeX = header['CTYPE3']
    #specTypeY = header['CTYPE3']
    #specUnitY = header['BUNIT']
    #specUnitX = header['CUNIT3']
    cubeDim = cube.shape
    
    for obj in objects:
      
	# centers and bounding boxes 
	#obj = np.array(objects[rr - 1])
	Xc = obj[cathead=='x'][0]
	Yc = obj[cathead=='y'][0]
	Zc = obj[cathead=='z'][0]
	Xmin = obj[cathead=='x_min'][0]
	Ymin = obj[cathead=='y_min'][0]
	Zmin = obj[cathead=='z_min'][0]
	Xmax = obj[cathead=='x_max'][0]
	Ymax = obj[cathead=='y_max'][0]
	Zmax = obj[cathead=='z_max'][0]
	
	# if center of mass estimation is wrong replace by geometric center
	if Xc < 0 or Xc > cubeDim[2]-1:  Xc = obj[cathead=='x_geo'][0]
	if Yc < 0 or Yc > cubeDim[1]-1:  Yc = obj[cathead=='y_geo'][0]
	if Zc < 0 or Zc > cubeDim[0]-1:  Zc = obj[cathead=='z_geo'][0]
    
        cPixXNew = int(Xc)
        cPixYNew = int(Yc)
        cPixZNew = int(Zc)
	
	# largest distance of source limits from the center
	maxX = 2*max(abs(cPixXNew-Xmin),abs(cPixXNew-Xmax))
	maxY = 2*max(abs(cPixYNew-Ymin),abs(cPixYNew-Ymax))
	maxZ = 2*max(abs(cPixZNew-Zmin),abs(cPixZNew-Zmax))
	
	# calculate the new bounding box for the mass centered cube
	XminNew = cPixXNew - maxX
	if XminNew < 0: XminNew = 0
	YminNew = cPixYNew - maxY
	if YminNew < 0: YminNew = 0
	ZminNew = cPixZNew - maxZ
	if ZminNew < 0: ZminNew = 0
	XmaxNew = cPixXNew + maxX
	if XmaxNew > cubeDim[2]-1: XmaxNew = cubeDim[2]-1
	YmaxNew = cPixYNew + maxY
	if YmaxNew > cubeDim[1]-1: YmaxNew = cubeDim[1]-1
	ZmaxNew = cPixZNew + maxZ
	if ZmaxNew > cubeDim[0]-1: ZmaxNew = cubeDim[0]-1
	
	# calculate the center with respect to the cutout cube
	cPixXCut = cPixX - XminNew
	cPixYCut = cPixY - YminNew
	cPixZCut = cPixZ - ZminNew
	
	# update header keywords:
	header['CRPIX1'] = cPixXCut+1
	header['CRPIX2'] = cPixYCut+1
	header['CRPIX3'] = cPixZCut+1
	
	# extract the cubelet
	subcube = cube[ZminNew:ZmaxNew+1,YminNew:YmaxNew+1,XminNew:XmaxNew+1]

	# update header keywords:
	header['NAXIS1'] = subcube.shape[2]
	header['NAXIS2'] = subcube.shape[1]
	header['NAXIS3'] = subcube.shape[0]

	# write the cubelet
	hdu = pyfits.PrimaryHDU(data=subcube,header=header)
	hdulist = pyfits.HDUList([hdu])
	name = outputDir + cubename + '_' + str(int(obj[0])) + '.fits'
	if compress:
		name += '.gz'
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		hdulist.writeto(name, output_verify = 'warn', clobber = True)
	
	hdulist.close()
	
	# make PV diagram
	if 'kin_pa' in cathead:
		pv_sampling=10
		pv_r=np.arange(-max(subcube.shape[1:]),max(subcube.shape[1:])-1+1./pv_sampling,1./pv_sampling)
		pv_y=Yc-YminNew+pv_r*math.cos(float(obj[cathead=='kin_pa'][0])/180*math.pi)
		pv_x=Xc-XminNew-pv_r*math.sin(float(obj[cathead=='kin_pa'][0])/180*math.pi)
		pv_x,pv_y=pv_x[(pv_x>=0)*(pv_x<=subcube.shape[2]-1)],pv_y[(pv_x>=0)*(pv_x<=subcube.shape[2]-1)]
		pv_x,pv_y=pv_x[(pv_y>=0)*(pv_y<=subcube.shape[1]-1)],pv_y[(pv_y>=0)*(pv_y<=subcube.shape[1]-1)]
		pv_x.resize((1,pv_x.shape[0]))
		pv_y.resize((pv_x.shape))
		pv_coords=np.concatenate((pv_y,pv_x),axis=0)
		pv_array=[]
		for jj in range(subcube.shape[0]):
			plane=map_coordinates(subcube[jj],pv_coords)
			plane=[plane[ii::pv_sampling] for ii in range(pv_sampling)]
			plane=np.array([ii[:plane[-1].shape[0]] for ii in plane])
			pv_array.append(plane.mean(axis=0))
		pv_array=np.array(pv_array)
		hdu = pyfits.PrimaryHDU(data=pv_array,header=header)
		hdulist = pyfits.HDUList([hdu])
		hdulist[0].header['CTYPE1']='PV--DIST'
		hdulist[0].header['CDELT1']=hdulist[0].header['CDELT2']
		hdulist[0].header['CRVAL1']=0
		hdulist[0].header['CRPIX1']=pv_array.shape[1]/2
		hdulist[0].header['CTYPE2']=hdulist[0].header['CTYPE3']
		hdulist[0].header['CDELT2']=hdulist[0].header['CDELT3']
		hdulist[0].header['CRVAL2']=hdulist[0].header['CRVAL3']
		hdulist[0].header['CRPIX2']=hdulist[0].header['CRPIX3']
		del hdulist[0].header['CTYPE3']
		del hdulist[0].header['CDELT3']
		del hdulist[0].header['CRVAL3']
		del hdulist[0].header['CRPIX3']
		name = outputDir+cubename+'_'+str(int(obj[0]))+'_pv.fits'
		
		# Check for overwrite flag:
		if not flagOverwrite and os.path.exists(name):
			sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
		else:
			hdulist.writeto(name,output_verify='warn',clobber=True)
		hdulist.close()

	# remove all other sources from the mask
	submask = mask[ZminNew:ZmaxNew+1,YminNew:YmaxNew+1,XminNew:XmaxNew+1].astype('int')
	submask[submask!=obj[0]] = 0
	submask[submask==obj[0]] = 1
	
	# write mask
	hdu = pyfits.PrimaryHDU(data=submask.astype('int16'),header=header)
	hdu.header['bunit']='source_ID'
	hdu.header['datamin']=submask.min()
	hdu.header['datamax']=submask.max()
	hdulist = pyfits.HDUList([hdu])
	name = outputDir+cubename+'_'+str(int(obj[0]))+'_mask.fits'
	if compress:
	  name += '.gz'
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		hdulist.writeto(name,output_verify='warn',clobber=True)
	hdulist.close()
	


	# moment 0
	if 'vopt' in header['ctype3'].lower() or 'vrad' in header['ctype3'].lower() or 'velo' in header['ctype3'].lower() or 'felo' in header['ctype3'].lower():
		if not 'cunit3' in header: dkms=abs(header['cdelt3'])/1e+3 # assuming m/s
		elif header['cunit3'].lower()=='m/s': dkms=abs(header['cdelt3'])/1e+3
		elif header['cunit3'].lower()=='km/s': dkms=abs(header['cdelt3'])
	elif 'freq' in header['ctype3'].lower():
		if not 'cunit3' in header or header['cunit3'].lower()=='hz': dkms=abs(header['cdelt3'])/1.42040575177e+9*2.99792458e+5 # assuming Hz
		elif header['cunit3'].lower()=='khz': dkms=abs(header['cdelt3'])/1.42040575177e+6*2.99792458e+5
		else: dkms=1. # no scaling, avoids crashing
	subcubeCopy = subcube.copy()
	subcubeCopy[submask==0] = 0
	if 'cellscal' in header:
		if header['cellscal'] == '1/F': subcubeCopy = regridMaskedChannels(subcubeCopy,submask,header)
	m0=np.nan_to_num(subcubeCopy).sum(axis=0)
	m0*=dkms
	hdu = pyfits.PrimaryHDU(data=m0,header=header)
	hdu.header['bunit']+='.km/s'
	hdu.header['datamin']=np.nanmin(m0)
	hdu.header['datamax']=np.nanmax(m0)
	del(hdu.header['crpix3'])
	del(hdu.header['crval3'])
	del(hdu.header['cdelt3'])
	del(hdu.header['ctype3'])
	name = outputDir+cubename+'_'+str(int(obj[0]))+'_mom0.fits'
	if compress:
	  name += '.gz'
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		hdu.writeto(name,output_verify='warn',clobber=True)
	
	
	# moment 1
	m1=(np.arange(subcubeCopy.shape[0]).reshape((subcubeCopy.shape[0],1,1))*np.ones(subcubeCopy.shape)-header['crpix3']+1)*header['cdelt3']+header['crval3']
	if 'vopt' in header['ctype3'].lower() or 'vrad' in header['ctype3'].lower() or 'velo' in header['ctype3'].lower() or 'felo' in header['ctype3'].lower():
		if not 'cunit3' in header: m1/=1e+3 # assuming m/s
		elif header['cunit3'].lower()=='km/s': pass
	elif 'freq' in header['ctype3'].lower():
		if not 'cunit3' in header or header['cunit3'].lower()=='hz': m1*=2.99792458e+5/1.42040575177e+9 # assuming Hz
		elif header['cunit3'].lower()=='khz': m1*=2.99792458e+5/1.42040575177e+6
	m0[m0==0]=np.nan
	m0/=dkms
	m1=np.divide(np.array(np.nan_to_num(m1*subcubeCopy).sum(axis=0)),m0)
	hdu = pyfits.PrimaryHDU(data=m1,header=header)
	hdu.header['bunit']='km/s'
	hdu.header['datamin']=np.nanmin(m1)
	hdu.header['datamax']=np.nanmax(m1)
	del(hdu.header['crpix3'])
	del(hdu.header['crval3'])
	del(hdu.header['cdelt3'])
	del(hdu.header['ctype3'])
	name = outputDir+cubename+'_'+str(int(obj[0]))+'_mom1.fits'
	if compress:
	  name += '.gz'
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		hdu.writeto(name,output_verify='warn',clobber=True)
	
	
	# moment 2
	m2=(np.arange(subcubeCopy.shape[0]).reshape((subcubeCopy.shape[0],1,1))*np.ones(subcubeCopy.shape)-header['crpix3']+1)*header['cdelt3']+header['crval3']
	if 'vopt' in header['ctype3'].lower() or 'vrad' in header['ctype3'].lower() or 'velo' in header['ctype3'].lower() or 'felo' in header['ctype3'].lower():
		if not 'cunit3' in header: m2/=1e+3 # assuming m/s
		elif header['cunit3'].lower()=='km/s': pass
	elif 'freq' in header['ctype3'].lower():
		if not 'cunit3' in header or header['cunit3'].lower()=='hz': m2*=2.99792458e+5/1.42040575177e+9 # assuming Hz
		elif header['cunit3'].lower()=='khz': m2*=2.99792458e+5/1.42040575177e+6
	m2 = m2**2
	m2=np.divide(np.array(np.nan_to_num(m2*subcubeCopy).sum(axis=0)),m0)
	m2-=m1*m1
	m2=np.sqrt(m2)
	hdu = pyfits.PrimaryHDU(data=m2,header=header)
	hdu.header['bunit']='km/s'
	hdu.header['datamin']=np.nanmin(m1)
	hdu.header['datamax']=np.nanmax(m1)
	del(hdu.header['crpix3'])
	del(hdu.header['crval3'])
	del(hdu.header['cdelt3'])
	del(hdu.header['ctype3'])
	name = outputDir+cubename+'_'+str(int(obj[0]))+'_mom2.fits'
	if compress:
		name += '.gz'
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		hdu.writeto(name,output_verify='warn',clobber=True)
	
	
	# spectra
	spec = np.nansum(subcubeCopy,axis=(1,2))
	
	name = outputDir + cubename + '_' + str(int(obj[0])) + '_spec.txt'
	if compress: name += '.gz'
	
	# Check for overwrite flag:
	if not flagOverwrite and os.path.exists(name):
		sys.stderr.write("ERROR: Output file exists: " + name + ".\n")
	else:
		if compress:
			import gzip
			f = gzip.open(name, 'wb')
		else:
			f = open(name, 'w')
		#f.write('# '+specTypeX+' ('+specUnitX+')'+'  '+specTypeY+' ('+specUnitY+')\n')
		
		for i in range(0,len(spec)):
			xspec = cValZ + (i+ZminNew-cPixZ) * dZ
			f.write('%9.0f %15.6e %15.6e\n'%(i+ZminNew,xspec,spec[i]))
		
		f.close()
