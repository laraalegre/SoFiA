#!/usr/bin/python
# -*- coding: utf-8 -*-

from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, tostring, XML
from xml.dom import minidom
from numpy import *
import sys
import os


# example page:
# http://pymotw.com/2/xml/etree/ElementTree/create.html

def prettify(elem):

    # Return a pretty-printed XML string for the Element.
    # Indent is set to "" here to save disk space; default would be "\t".
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent='')


def make_xml(objects, outname, flagOverwrite):
    print 'Store the results to xml file: ', outname
    top = Element('VOTABLE')
    resource = SubElement(top, 'RESOURCE', name='SoFiA catalogue')
    description = SubElement(resource, 'DESCRIPTION')
    description.text = 'Output catalogue from SoFiA'
    coocys = SubElement(resource, 'COOSYS', ID='J2000')
    table = SubElement(resource, 'TABLE', ID='sofia_cat', name='sofia_cat')
    description = SubElement(table, 'DESCRIPTION')
    description.text = 'Output catalogue from SoFiA'

    # write the parameters in fields:
    sources = objects.getSources()
    pars = sources.itervalues().next().getParameters()
    for i in pars:
        # this has to be updated with proper types and units
        field = SubElement(table, 'FIELD', name=i, datatype='float', unit='-')

    # write the data....
    data = SubElement(table, 'DATA')
    tabledata = SubElement(data, 'TABLEDATA')

    for i in sources:
        tr = SubElement(tabledata, 'TR')
        source_dict = sources[i].getParameters()
        for j in source_dict:
            td = SubElement(tr, 'TD')
            td.text = str(source_dict[j].getValue())

    # print prettify(top)
    # print

    # Check for overwrite flag:
    if not flagOverwrite and os.path.exists(outname):
        sys.stderr.write("ERROR: Output file exists: " + outname + ".\n")
    else:
        f1 = open(outname, 'w+')
        f1.write(prettify(top))
        f1.close

    return


def make_xml_from_array(objects, cathead, catunits, catfmt, store_pars, outname, compress, flagOverwrite):
    print 'Store the results to xml file: ', outname
    top = Element('VOTABLE')
    resource = SubElement(top, 'RESOURCE', name='SoFiA catalogue')
    description = SubElement(resource, 'DESCRIPTION')
    description.text = 'Output catalogue from SoFiA'
    coocys = SubElement(resource, 'COOSYS', ID='J2000')
    table = SubElement(resource, 'TABLE', ID='sofia_cat', name='sofia_cat')
    description = SubElement(table, 'DESCRIPTION')
    description.text = 'Output catalogue from SoFiA'
    
    # Load list of parameters and unified content descriptors (UCDs)
    parList = {};
    fileUcdPath = os.environ['SOFIA_PIPELINE_PATH'];
    fileUcdPath = fileUcdPath.replace("sofia_pipeline.py", "SoFiA_source_parameters.dat");
    
    try:
        with open(fileUcdPath) as fileUcd:
            for line in fileUcd:
               (key, value) = line.split();
               parList[key] = value;
    except:
        sys.stderr.write("WARNING: Failed to read UCD file.\n");

    # write the parameters in fields:
    if store_pars == ['*']:
        for i in cathead:
            if (i in parList): ucdEntity = parList[i];
            else: ucdEntity = "";
            field = SubElement(
                table, 'FIELD', name=i, ucd=ucdEntity, datatype='float', unit=catunits[cathead.index(i)]
                )
    else:
        for par in store_pars:
            if (par in cathead):
                index = list(cathead).index(par)
                if (cathead[index] in parList): ucdEntity = parList[cathead[index]];
                else: ucdEntity = "";
                field = SubElement(table, 'FIELD', name=cathead[index], ucd=ucdEntity, datatype='float', unit=catunits[index])
            else:
                sys.stderr.write("WARNING: Skipping undefined parameter \'" + str(par) + "\'.\n");

    # write the data....
    data = SubElement(table, 'DATA')
    tabledata = SubElement(data, 'TABLEDATA')

    if store_pars == ['*']:
        for obj in objects:
            tr = SubElement(tabledata, 'TR')
            for i in range(0, len(obj)):
                td = SubElement(tr, 'TD')
                td.text = (catfmt[i] % obj[i]).strip()
    else:
        for obj in objects:
            tr = SubElement(tabledata, 'TR')
            for par in store_pars:
                if (par in cathead):
                    td = SubElement(tr, 'TD')
                    index = list(cathead).index(par)
                    td.text = (catfmt[index] % obj[index]).strip()

    if compress: outname += '.gz'

    # Check for overwrite flag:
    if not flagOverwrite and os.path.exists(outname):
        sys.stderr.write("ERROR: Output file exists: " + outname + ".\n")
    else:
        if compress:
            import gzip
            f1 = gzip.open(outname, 'wb')
        else:
            f1 = open(outname, 'w+')
        f1.write(prettify(top))
        f1.close

    return
