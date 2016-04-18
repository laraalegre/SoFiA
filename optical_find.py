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

def merge_ascii(outroot,cat,storeMultiCat):
    # merge the ascii tables:
    outfile = outroot + '_cat.ascii'
    f=open(outfile,'w')
    h = 0    
    for i in range(len(cat)):
        temp_ascii = outroot + '_' + cat[i]['id'] + '_cat.ascii'
        # check whether the file exists:
        n = 0
        if os.path.isfile(temp_ascii) == True:
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
            if storeMultiCat == False:
            # security check whether this is a legitimate file (and not a link) before deleting anything
                if os.path.isfile(temp_ascii) == True and os.path.islink(temp_ascii) == False:
                    os.system('rm -f ' + temp_ascii)
                else:
                    sys.stderr.write("ERROR: The specified file does not seem to be valid \n")
                    sys.stderr.write("       Cannot identify: " + temp_ascii + " as a normal file\n")
                    raise SystemExit(1)

            
    f.close    
    return



def merge_xml(outroot,cat,storeMultiCat):
    # merge the xml tables:
    # this is a temporary (working) solution, should be replaced by proper xml functions in python
    outfile = outroot + '_cat.xml'
    f=open(outfile,'w')


    hlines = 0
    cat_nr = 0
    #h = 0    
    for i in range(len(cat)):
        temp_xml = outroot + '_' + cat[i]['id'] + '_cat.xml'
        # check whether the file exists:
        #n = 0
        if os.path.isfile(temp_xml) == True:
            cat_nr += 1
            f1=open(temp_xml,'r')
            if cat_nr == 1:
                # this is the first catalogue which is used for the meta information
                stop_header = 0
                # count the number of header lines and write the header
                for line in f1:
                    line_string = line
                                        
                    if stop_header == 0 and line_string[0:4] != '<TR>':
                        hlines += 1
                        f.write(line_string)
                    else:
                        stop_header = 1
                        if line_string[0:2] != '</' or line_string[0:5] == '</TR>':
                            f.write(line_string)
        
        
            else:
                # from the other catalogues only data is used
                l = 0
                for line in f1:
                    line_string = line                    
                    l+=1
                    if l > hlines and line_string[0:2] != '</' or line_string[0:5] == '</TR>':
                        f.write(line_string)
            
                                 
            
            f1.close
            # throw away the intermediate catalogues if not needed
            if storeMultiCat == False:
                # security check whether this is a legitimate file (and not a link) before deleting anything
                if os.path.isfile(temp_xml) == True and os.path.islink(temp_xml) == False:
                    os.system('rm -f ' + temp_xml)
                else:
                    sys.stderr.write("ERROR: The specified file does not seem to be valid \n")
                    sys.stderr.write("       Cannot identify: " + temp_xml + " as a normal file\n")
                    raise SystemExit(1)

                              
    # write the meta information at the end        
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


match_catalogue = outroot + '_match.csv'
    
cube = Parameters['import']['inFile']
specSize = Parameters['optical']['specSize']
spatSize = Parameters['optical']['spatSize']
catalogue = Parameters['optical']['sourceCatalogue']
storeMultiCat = Parameters['optical']['storeMultiCat']


print 'working on cube: ', cube
print 'working on catalogue: ', catalogue

# check whether the catalogue exists, and otherwise exit.
if os.path.isfile(catalogue) == False:
    sys.stderr.write("ERROR: The specified source catalogue does not exist.\n")
    sys.stderr.write("       Cannot find: " + catalogue + "\n")
    raise SystemExit(1)

# check whether the catalogue is a csv file, and otherwise exit
if catalogue[-4:] != '.csv':
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
if 'id' not in cat[0] or 'ra' not in cat[0] or 'dec' not in cat[0] or 'z' not in cat[0]: 
    sys.stderr.write("ERROR: The specified source catalogue does not have the right input\n")
    sys.stderr.write("       it should contain at least for columns with values id, ra, dec and z \n")
    sys.stderr.write("       Cannot work on: " + catalogue + "\n")
    raise SystemExit(1)


# open a file for the output catalogue
f_out = open(match_catalogue,'w')
f_out.write('# id, ra, dec, z, match \n')


for i in range(len(cat)):
    # define the subregion:
    subcube  = [float(cat[i]['ra']),float(cat[i]['dec']),float(cat[i]['z']),spatSize,spatSize,specSize]


    new_file = 'temp_optical.txt'
    # replace the relevant parameters in the parameter file
    f1 = open(sys.argv[1],'r')
    f2 = open(new_file,'wt')

    for line in f1:
        pattern0 = 'steps.doSubcube'
        subst0 = 'steps.doSubcube =  true \n'
        pattern1 = 'import.subcube'
        subst1 = 'import.subcube = ' + str(subcube) + '\n'          
        pattern2 = 'import.subcubeMode'
        subst2 = 'import.subcubeMode =  ' + 'world' + '\n'           
        pattern3 = 'writeCat.basename'
        if (Parameters['writeCat']['basename']): 
            subst3 = 'writeCat.basename = ' + Parameters['writeCat']['basename'] + '_' + cat[i]['id'] + '\n'
        else:
            temp_base = Parameters['import']['inFile'].split('/')[-1]
            if((temp_base.lower()).endswith(".fits") and len(temp_base) > 5):
                temp_base = temp_base[0:-5]
            subst3 = 'writeCat.basename = ' + temp_base + '_' + cat[i]['id'] + '\n'
        if pattern0 in line:
            f2.write(subst0)
        elif pattern1 in line and pattern2 not in line:
            f2.write(subst1)
        elif pattern2 in line:
            f2.write(subst2)
        elif pattern3 in line:
            f2.write(subst3)
        else:
            f2.write(line)

    f1.close()    
    f2.close()


    os.system('python $SOFIA_PIPELINE_PATH temp_optical.txt')
    #os.system('cp temp_optical.txt temp_'+ cat[i]['id'] + '.txt' )
    os.system('rm -f temp_optical.txt')

    # write the number of detections into the match_catalogue
    temp_ascii = outroot + '_' + cat[i]['id'] + '_cat.ascii'

    # check whether the file exists:
    n = 0
    nr_sources = 0
    if os.path.exists(temp_ascii) == True:
        f1=open(temp_ascii,'r')
        for line in f1:
            n += 1
            if  n > 4:
                nr_sources += 1
        
    f_out.write(cat[i]['id'] + ',' + cat[i]['ra'] + ',' + cat[i]['dec'] + ',' + cat[i]['z'] + ',' +str(int( nr_sources)) + '\n')
    f1.close
    
    
f_out.close()







    
# --------------------
# ---- STORE DATA ----

# --------------------

if Parameters['steps']['doWriteCat'] and Parameters['steps']['doMerge']:
    print "\n--- SoFiA: Writing output catalogue ---"
    sys.stdout.flush()
    if Parameters['writeCat']['writeXML'] and Parameters['steps']['doMerge']:
        print 'Merge the xml file'
        merge_xml(outroot,cat,storeMultiCat)
				
    if Parameters['writeCat']['writeASCII'] and Parameters['steps']['doMerge']:
        print 'Merge the ascii file'
        merge_ascii(outroot,cat,storeMultiCat)
				



        

print "\n--- SoFiA: Pipeline finished ---"
sys.stdout.flush()


