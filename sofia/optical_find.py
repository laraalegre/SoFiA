#! /usr/bin/env python

# A function to do optically based source finding.
# Input are an data cube and optical calatogue
# this task defines a subcube around every optical source and does the source finding

import csv
import re
import os

def optical_find(cube,catalogue,defaultPars,spatSize,specSize,subcubeMode,outputCat,storeSingleCat):
    # cube: the indput cube
    # catalogue: a source catalogue with RA, DEC and Freq
    # defaultPars: the parameter file with default parameters
    # spatSize: the spatial size of the subcube (deg)
    # specSize: the spectral size of the subcube (unit of the cube)
    # outputCat: the name of the output catalogue

    # cube: numpy arrau
    # catalogue: string
    # spatSize: float
    # specSize: float
    # outputCat: string
    # storeSingleCat: bool

    
    print 'working on cube: ', cube
    print 'working on catalogue: ', catalogue


    # read the optical catalogue
    # For the moment I require a .csv file, with the following order:
    # id, ra (degrees), dec (degrees) and z (units of cube)    

    f = open(catalogue,'rb')
    reader = csv.DictReader(f)
    cat = list(reader)
    f.close()

         
    for i in range(len(cat)):
        # define the subregion:
        subcube  = [float(cat[i]['ra']),float(cat[i]['dec']),float(cat[i]['z']),spatSize,spatSize,specSize]
        
        
        new_file = 'temp_optical.txt'
        # replace the relevant parameters in the parameter file
        f1 = open(defaultPars,'r')
        f2 = open(new_file,'wt')
        
        for line in f1:
            pattern1 = 'import.subcube '
            subst1 = 'import.subcube = ' + str(subcube) + '\n'          
            pattern2 = 'import.subcubeMode'
            subst2 = 'import.subcubeMode =  ' + subcubeMode + '\n'           
            pattern3 = 'writeCat.basename'
            subst3 = 'writeCat.basename = ' + cat[i]['id'] + '\n'
            if pattern1 in line:
                f2.write(subst1)
            elif pattern2 in line:
                f2.write(subst2)
            elif pattern3 in line:
                f2.write(subst3)
            else:
                f2.write(line)
        
        f1.close()    
        f2.close()

        #os.system('pwd')
        #os.system('less temp_optical.txt')
        os.system('python SoFiA-master/sofia_pipeline.py temp_optical.txt')


    # merge the output tables:
    f=open(outputCat,'w')
    h = 0    
    for i in range(len(cat)):
        temp_ascii = cat[i]['id'] + '_cat.ascii'
        # check whether the file exists:
        n = 0
        if os.path.exists(temp_ascii) == True:
            f1=open(temp_ascii,'r')
            for line in f1:
                h += 1
                if h < 5:
                    f.write(line)  
        
                n += 1
                if  n > 4:
                    f.write(line)
            f1.close
            # throw away the intermediate catalogues if not needed
            if storeSingleCat == False:
                os.system('rm ' + temp_ascii)
    f.close    
    return

##########################
##########################


#catalogue = 'test_cat.csv'
#cube = 'nancube.fits'
#outputCat = 'full_cat.txt'
#spatSize = 10
#specSize = 10
#subcubeMode = 'pix'
#storeSingleCat = False
#defaultPars = 'par_optical.txt'


#optical_find(cube,catalogue,defaultPars,spatSize,specSize,subcubeMode,outputCat,storeSingleCat)
