#! /usr/bin/env python

# A function to do optically based source finding.
# Input are an data cube and optical calatogue
# this task defines a subcube around every optical source and does the source finding

import csv
import re
import os
import sys

sys.path.insert(0, os.environ['SOFIA_MODULE_PATH'])
from sofia import readoptions



# cube: the indput cube
# catalogue: a source catalogue with RA, DEC and Freq
# defaultPars: the parameter file with default parameters
# spatSize: the spatial size of the subcube (deg)
# specSize: the spectral size of the subcube (unit of the cube)

# cube: numpy arrau
# catalogue: string
# spatSize: float
# specSize: float
# subcubeMode: str
# storeSingleCat: bool


# ------------------------------------
# ---- Functions to merger tables ----
# ------------------------------------

def merge_ascii(outroot,cat,storeSingleCat):
    # merge the ascii tables:
    outfile = outroot + '_cat.ascii'
    f=open(outfile,'w')
    h = 0    
    for i in range(len(cat)):
        temp_ascii = outroot + '_' + cat[i]['id'] + '_cat.ascii'
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
                os.system('rm -rf ' + temp_ascii)           
    f.close    
    return



def merge_xml(outroot,cat,storeSingleCat):
    # merge the xml tables:
    outfile = outroot + '_cat.xml'
    f=open(outfile,'w')
    h = 0    
    for i in range(len(cat)):
        temp_xml = outroot + '_' + cat[i]['id'] + '_cat.xml'
        # check whether the file exists:
        n = 0
        if os.path.exists(temp_xml) == True:
            f1=open(temp_xml,'r')
            for line in f1:
                h += 1
                if h < 39:
                    f.write(line)  

                n += 1
                if  n > 38:
                    if line[0:2] != '</':
                    
                        f.write(line)
            f.write('</TR>\n')            
            f1.close
            # throw away the intermediate catalogues if not needed
            if storeSingleCat == False:
                os.system('rm -rf ' + temp_xml)           
            
    f.write('</TABLEDATA>\n')
    f.write('</DATA>\n')
    f.write('</TABLE>\n')
    f.write('</RESOURCE>\n')
    f.write('</VOTABLE>\n')
    f.close    
    return

######################################
######################################
######################################


# ---------------------------------
# ---- READ DEFAULT PARAMETERS ----
# ---------------------------------

print "\n--- SoFiA: Reading default parameters ---"
sys.stdout.flush()

# This reads in the default parameters:
default_file = '%s/SoFiA_default_input.txt'%(os.path.dirname(os.path.realpath(__file__)))
#default_file = 'SoFiA_default_input.txt'
Parameters = readoptions.readPipelineOptions(default_file)

# ------------------------------
# ---- READ USER PARAMETERS ----
# ------------------------------

# this is neede to know what kind of output is required

print "\n--- SoFiA: Reading user parameters ---"
sys.stdout.flush()

# This reads in a file with parameters and creates a dictionary:
parameter_file = sys.argv[1]
print 'Parameters extracted from: ', parameter_file
print
User_Parameters = readoptions.readPipelineOptions(parameter_file)

# Overwrite default parameters with user parameters (if exist):
for task in User_Parameters.iterkeys():
        if(task in Parameters):
                for key in User_Parameters[task].iterkeys():
                        if(key in Parameters[task]):
                                Parameters[task][key] = User_Parameters[task][key]

# Define the base and output directory name used for output files (defaults to input file 
# name if writeCat.basename is found to be invalid):
outputBase = Parameters['writeCat']['basename']
outputDir  = Parameters['writeCat']['outputDir']

if((not outputBase) or outputBase.isspace() or ("/" in outputBase) or ("\\" in outputBase) or (outputBase == ".") or (outputBase == "..")):
    outroot = Parameters['import']['inFile'].split('/')[-1]
    if((outroot.lower()).endswith(".fits") and len(outroot) > 5):
        outroot = outroot[0:-5]
else:
    outroot = outputBase

if((not outputDir) or (not os.path.isdir(outputDir)) or (outputDir.isspace())):
    outroot = Parameters['import']['inFile'][0:len(Parameters['import']['inFile'])-len(Parameters['import']['inFile'].split('/')[-1])]+outroot
else:
    if outputDir[-1] != '/': outputDir += '/'
    outroot = outputDir + outroot



cube = Parameters['import']['inFile']
specSize = Parameters['optical']['specSize']
spatSize = Parameters['optical']['spatSize']
catalogue = Parameters['optical']['sourceCatalogue']
storeSingleCat = Parameters['optical']['storeSingleCat']


print 'working on cube: ', cube
print 'working on catalogue: ', catalogue

# check whether the catalogue exists, and otherwise exit.
if os.path.isfile(catalogue) == False:
    sys.stderr.write("ERROR: The specified source catalogue does not exist.\n")
    sys.stderr.write("       Cannot find: " + catalogue + "\n")
    raise SystemExit(1)

# check whether the catalogue is a csv file, and otherwise exit
# not the most elegant "if" but somehow it doesn't work at once
if catalogue[-4:-1] != '.cs' and catalogue[-1] != 'v':
    sys.stderr.write("ERROR: The specified source catalogue is not a .csv file.\n")
    sys.stderr.write("       Cannot work on: " + catalogue + "\n")
    raise SystemExit(1)



# read the optical catalogue
# For the moment I require a .csv file, with the following order:
# id, ra (degrees), dec (degrees) and z (units of cube)    

f = open(catalogue,'rb')
reader = csv.DictReader(f)
cat = list(reader)
f.close()

# check whether the file has the right columns (id, ra, dec, z)
if cat[0]['id'] == False or cat[0]['ra'] == False or cat[0]['dec'] == False or cat[0]['z'] == False: 
    sys.stderr.write("ERROR: The specified source catalogue does not have the right input\n")
    sys.stderr.write("       it should contain at least for columns with values id, ra, dec and z \n")
    sys.stderr.write("       Cannot work on: " + catalogue + "\n")
    raise SystemExit(1)

for i in range(len(cat)):
    # define the subregion:
    subcube  = [float(cat[i]['ra']),float(cat[i]['dec']),float(cat[i]['z']),spatSize,spatSize,specSize]


    new_file = 'temp_optical.txt'
    # replace the relevant parameters in the parameter file
    f1 = open(sys.argv[1],'r')
    f2 = open(new_file,'wt')

    for line in f1:
        pattern1 = 'import.subcube '
        subst1 = 'import.subcube = ' + str(subcube) + '\n'          
        pattern2 = 'import.subcubeMode'
        subst2 = 'import.subcubeMode =  ' + 'wcs' + '\n'           
        pattern3 = 'writeCat.basename'
        subst3 = 'writeCat.basename = ' + Parameters['writeCat']['basename'] + '_' + cat[i]['id'] + '\n'
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

    os.system('python SoFiA-master/sofia_pipeline.py temp_optical.txt')
    os.system('rm -rf temp_optical.txt')









    
# --------------------
# ---- STORE DATA ----
# --------------------

if Parameters['steps']['doWriteCat'] and Parameters['steps']['doMerge']:
    print "\n--- SoFiA: Writing output catalogue ---"
    sys.stdout.flush()
    if Parameters['writeCat']['writeXML'] and Parameters['steps']['doMerge']:
        print 'Merge the xml file'
        merge_xml(outroot,cat,storeSingleCat)
				
    if Parameters['writeCat']['writeASCII'] and Parameters['steps']['doMerge']:
        print 'Merge the ascii file'
        merge_ascii(outroot,cat,storeSingleCat)
				



        

print "\n--- SoFiA: Pipeline finished ---"
sys.stdout.flush()


