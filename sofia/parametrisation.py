#!/usr/bin/python
# -*- coding: utf-8 -*-

from sofia import cparametrizer as cp
import numpy as np
import scipy.ndimage as nd
import sys

def dilate(cube,mask,objects,cathead,Parameters):
    dilate_threshold=Parameters['parameters']['dilate_threshold']
    dilate_pix_max=Parameters['parameters']['dilate_pix_max']
    dilate_chan=Parameters['parameters']['dilate_chan']
    # stops dilating when (flux_new-flux_old)/flux_new < dilate_threshold
    for mm in range(1,mask.max()+1):
    	obj=objects[mm-1]
    	xmin=obj[list(cathead).index('Xmin')]-dilate_pix_max
    	xmax=obj[list(cathead).index('Xmax')]+dilate_pix_max
    	ymin=obj[list(cathead).index('Ymin')]-dilate_pix_max
    	ymax=obj[list(cathead).index('Ymax')]+dilate_pix_max
    	zmin=obj[list(cathead).index('Zmin')]-dilate_chan
    	zmax=obj[list(cathead).index('Zmax')]+dilate_chan
    	xmin=max(0,xmin)
    	xmax=min(xmax,cube.shape[2]-1)
    	ymin=max(0,ymin)
    	ymax=min(ymax,cube.shape[1]-1)
    	zmin=max(0,zmin)
    	zmax=min(zmax,cube.shape[0]-1)
    	objcube=cube[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].copy()
    	objmask=mask[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].copy()
    	allmask=mask[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1].copy()
    	otherobjs=(allmask>0)*(allmask!=mm)
    	if (otherobjs).sum():
    		# Ensure that objects!=mm within dilate_pix_max, dilate_chan are not included in the flux growth calculation
    		print 'WARNING: object %i has possible overlapping objects within %i pix, %i chan'%(mm,dilate_pix_max,dilate_chan)
    		objcube[(allmask>0)*(allmask!=mm)]=0
        fluxes=[]
        for dil in range(dilate_pix_max+1):
            dd=dil*2+1
            dilstruct=(np.sqrt(((np.indices((dd,dd))-dil)**2).sum(axis=0))<=dil).astype(int)
            dilstruct.resize((1,dilstruct.shape[0],dilstruct.shape[1]))
            dilstruct=dilstruct.repeat(dilate_chan*2+1,axis=0)
            fluxes.append(objcube[nd.morphology.binary_dilation(objmask==mm,structure=dilstruct)].sum())
            if dil>0 and (fluxes[-1]-fluxes[-2])/fluxes[-1]<dilate_threshold: break
        # pick the best dilation kernel for current object and update mask
        dil-=1
        print 'Mask dilation of source %i by %i pix and %i chan'%(mm,dil,dilate_chan)
        sys.stdout.flush()
        dd=dil*2+1
        dilstruct=(np.sqrt(((np.indices((dd,dd))-dil)**2).sum(axis=0))<=dil).astype(int)
        dilstruct.resize((1,dilstruct.shape[0],dilstruct.shape[1]))
        dilstruct=dilstruct.repeat(dilate_chan*2+1,axis=0)
        # Only grow the mask of object mm even when other objects are present in objmask
        objmask[nd.morphology.binary_dilation(objmask==mm,structure=dilstruct).astype(int)==1]=mm
        # Put back in objmask objects!=mm that may have been inside objmask before dilation or may have been temporarily replaced by the dilated object mm
        if (otherobjs).sum(): objmask[otherobjs]=allmask[otherobjs]
        mask[zmin:zmax+1,ymin:ymax+1,xmin:xmax+1]=objmask
        del(objcube)
        del(objmask)
        del(allmask)
        del(otherobjs)
    return mask

def parametrise(
        cube,
        mask,
        objects,
        cathead,
        catformt,
        catparunits,
        Parameters,
        dunits,
        ):

    cathead = np.array(cathead)
    objects = np.array(objects)
    initcatalog = cp.PySourceCatalog()
    for obj in objects:
        newSource = cp.PySource()
        newSource.ID = obj[cathead == 'ID']
        newParamsDict = {
            'X': cp.PyMeasurement('X', obj[cathead == 'Xg'], 0., ''),
            'Y': cp.PyMeasurement('Y', obj[cathead == 'Yg'], 0., ''),
            'Z': cp.PyMeasurement('Z', obj[cathead == 'Zg'], 0., ''),
            'BBOX_X_MIN': cp.PyMeasurement('BBOX_X_MIN', obj[cathead == 'Xmin'], 0., ''),
            'BBOX_X_MAX': cp.PyMeasurement('BBOX_X_MAX', obj[cathead == 'Xmax'], 0., ''),
            'BBOX_Y_MIN': cp.PyMeasurement('BBOX_Y_MIN', obj[cathead == 'Ymin'], 0., ''),
            'BBOX_Y_MAX': cp.PyMeasurement('BBOX_Y_MAX', obj[cathead == 'Ymax'], 0., ''),
            'BBOX_Z_MIN': cp.PyMeasurement('BBOX_Z_MIN', obj[cathead == 'Zmin'], 0., ''),
            'BBOX_Z_MAX': cp.PyMeasurement('BBOX_Z_MAX', obj[cathead == 'Zmax'], 0., ''),
            }
        newSource.setParameters(newParamsDict)
        initcatalog.insert(newSource)

    moduleParametrizer = cp.PyModuleParametrisation()
    moduleParametrizer.setFlags(
        Parameters['parameters']['optimiseMask'],
        Parameters['parameters']['fitBusyFunction']
        )

    cube = cube.astype('<f4', copy=False)
    mask = mask.astype('<i2', copy=False)

    moduleParametrizer.run(cube, mask, initcatalog)
    results = moduleParametrizer.getCatalog()

    # append the results to the objects array or reset
    replParam = [
        'BBOX_X_MAX',
        'BBOX_X_MIN',
        'BBOX_Y_MAX',
        'BBOX_Y_MIN',
        'BBOX_Z_MAX',
        'BBOX_Z_MIN',
        'ID',
        'X',
        'Y',
        'Z',
        'NRvox',
        ]
    origParam = [
        'Xmax',
        'Xmin',
        'Ymax',
        'Ymin',
        'Zmax',
        'Zmin',
        'ID',
        'Xm',
        'Ym',
        'Zm',
        'NRvox',
        ]
    d = results.getSources()

    # select data set with maximum number of parameters
    parsListLen = [len(d[d.keys()[i]].getParameters()) for i in range(0, len(d))]
    index = parsListLen.index(max(parsListLen))

    # add parameter names from parametrization
    pars = d[d.keys()[index]].getParameters()
    cathead = list(cathead)
    newunits = {
        'ID': '-',
        'X': 'pix',
        'Y': 'pix',
        'Z': 'pix',
        'BBOX_X_MAX': 'pix',
        'BBOX_X_MIN': 'pix',
        'BBOX_Y_MIN': 'pix',
        'BBOX_Y_MAX': 'pix',
        'BBOX_Z_MIN': 'chan',
        'BBOX_Z_MAX': 'chan',
        'W50': 'chan',
        'W20': 'chan',
        'Wm50': 'chan',
        'F_Wm50': dunits,
        'ELL_MAJ': 'pix',
        'ELL_MIN': 'pix',
        'ELL_PA': 'deg',
        'F_TOT': dunits,
        'BF_FLAG': '-',
        'BF_CHI2': '-',
        'BF_Z': 'chan',
        'BF_A': dunits,
        'BF_B1': 'chan**(-1)',
        'BF_B2': 'chan**(-1)',
        'BF_C': 'chan**(-2)',
        'BF_XE0': 'chan',
        'BF_XP0': 'chan',
        'BF_W': 'chan',
        'BF_W50': 'chan',
        'BF_W20': 'chan',
        'BF_F_PEAK': dunits,
        'BF_F_INT': dunits,
        'RMS_CUBE': dunits,
        'F_PEAK': dunits,
        }
    catformt = list(catformt)
    catparunits = list(catparunits)
    for i in sorted(pars):
        if i not in replParam:
            cathead.append(i)
            catformt.append('%12.4f')
            catparunits.append(newunits[i])

    # extend the parameter array
    tmpObjects = np.empty((objects.shape[0], len(cathead)))
    tmpObjects[:, :] = np.NAN
    tmpObjects[:, 0:objects.shape[1]] = objects
    objects = tmpObjects
    for i in d:
        source_dict = d[i].getParameters()

        # check the source index
        index = int(source_dict['ID'].getValue())
        for j in sorted(source_dict):
            if j in replParam:
                objects[index - 1][cathead.index(origParam[replParam.index(j)])] = \
                    source_dict[j].getValue()
            else:
                objects[index - 1][cathead.index(j)] = source_dict[j].getValue()

    objects = np.array(objects)
    cathead = np.array(cathead)
    catparunits = np.array(catparunits)
    catformt = np.array(catformt)
    print 'The parameterisation has completed'
    print

    # print objects.shape,cathead.shape
    return (cube, mask, objects, cathead, catformt, catparunits)


