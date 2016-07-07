SoFiA
=====

Introduction
------------

SoFiA, the Source Finding Application, is a new HI source finding pipeline 
intended to find and parametrise galaxies in HI data cubes. SoFiA can be 
launched from the command line, but it also comes with an easy-to-use 
graphical user interface that allows control parameters to be manipulated 
interactively. While the software is still under development, several 
stable versions of SoFiA have already been released and can be obtained 
from the SoFiA webpage on GitHub.

If you would like to stay informed about new stable releases of SoFiA and 
other important updates, you can sign up to the **SoFiA mailing list**. To do 
so, simply send an e-mail to `sofia-request [at] atnf.csiro.au` with the 
word `subscribe` in the e-mail body (note that the e-mail subject will be 
ignored).


Requirements
------------

The following packages and libraries will be required to install and run 
SoFiA:

    Operating system:
        Linux or Unix (e.g. Ubuntu, Mac OS, etc.)
        Terminal with bash or tcsh (other shells should work as well)
    Packages (Python):
        Python (≥ 2.4)
        numpy (≥ 1.8)
        scipy (≥ 0.7)
        astropy (≥ 0.2.5)
        matplotlib (≥ 1.1; optional, only needed for reliability plot)
    Packages (C++):
        GCC (≥ 4.6)
        GNU Scientific Library (≥ 1.15; including dev files)
        Qt (≥ 4.7; including dev files and qmake)

All of the above packages must be installed before SoFiA can be compiled and 
run. It is recommended that you install them through your operating system’s 
package manager. Please ensure that the development packages of the GNU 
Scientific Library and the Qt library are installed as well.


Installation
------------

To install SoFiA on your computer, open a terminal window and change into the 
folder where the downloaded file was saved. Then follow the steps below.

1. Unpack all files

   Download and unpack the zipped archive into a directory of your choice:

   > unzip SoFiA-[version].zip

   or

   > tar -xzvf SoFiA-[version].tar.gz

   where [version] is the downloaded version of SoFiA, e.g. 0.2. This will 
   unpack all files into a directory called SoFiA-[version].

2. Enter the installation directory

   > cd SoFiA-[version]

3. Compile and install the SoFiA pipeline and user interface

   > python setup.py build --force
   
   (add `--no-gui=True` if you do not wish to install the graphical user interface)

4. Set up environment variables

    Follow the instructions given at the end of the installation process to 
    define the required environment variables and paths in your .bashrc or 
    .cshrc file.

5. Launch SoFiA

    Open a new terminal window and type:

    > SoFiA &


Problems
--------

If you encounter problems when running the setup script, it is likely 
that you are either missing one of the required packages and libraries or 
that some of the packages are outdated. Please check that the required 
versions of all packages are installed and properly set up such that they 
can be found by Python and the GCC. On some systems it may be necessary to 
explicitly install the GNU C++ compiler (g++, should normally be part of 
the GCC) as well as the development packages (dev) of Qt and the GNU 
Scientific Library. Please also see the trouble shooting page on the SoFiA 
wiki for more information on a few commonly encountered problems:

* https://github.com/SoFiA-Admin/SoFiA/wiki/SoFiA-Troubleshooting


Documentation
-------------

SoFiA comes with its own built-in help browser that can be launched from 
the help menu in the user interface. Alternatively, you can use your web 
browser to open the index.html file located in the gui/doc/ sub-folder.


Version history
---------------

* SoFiA 1.0.0
  * Released 12/07/2016
  * Adds several new features, including improvements in reliability
    calculation, the inclusion of new source parameters aimed at exteded
    galaxies, and the introduction of basic GUI configuration options.
  * Fixes a large number of bugs across the entire pipeline that have
    led to frequent crashes in the past.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v1.0.0

* SoFiA 0.5.0
  * Released 22/09/2015
  * Adds several new features, including GUI sessions, automated kernel
    size determination in reliability calculation and Unified Content
    Descriptors for source parameters (UCDs).
  * Fixes a few bugs in handling WCS information and merging detections
    into sources.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.5.0

* SoFiA 0.4.0
  * Released 22/12/2014
  * Adds several new features, including 2D–1D wavelet filtering, the CNHI 
    finder, and optically motivated source finding.
  * Fixes a few bugs and greatly improves the graphical user interface and 
    documentation.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.4.0

* SoFiA 0.3.2
  * Released 25/08/2014
  * Adds several new features, including subcube processing and separate 
    output directory.
  * Fixes a bug that limited the number of detected sources to 10,000.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.3.2

* SoFiA 0.3.1
  * Released 14/07/2014
  * Adds several new features, including weights functions, flagging, and
    Qt 5 support.
  * Fixes several bugs, including reduction of memory footprint and removal
    of a security vulnerability.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.3.1

* SoFiA 0.3
  * Released 20/05/2014
  * Adds several new features, including WCS support, PV diagrams, and 
    handling of units.
  * Fixes a few bugs related to FITS import and handling of blanked pixels.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.3

* SoFiA 0.2.1
  * Released 25/02/2014
  * Some more bug fixes, including fixing the threshold finder.
  * Implementation of VO table output.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.2.1

* SoFiA 0.2
  * Released 10/02/2014
  * Fixes several bugs that were still present in the initial beta version.
  * See release notes for details:
    https://github.com/SoFiA-Admin/SoFiA/releases/tag/v0.2

* SoFiA 0.1 beta
  * Released 14/11/2013
  * Initial beta release of the package with basic functionality and built-in 
    help system; may still contain bugs and errors.


Copyright and licence
---------------------

SoFiA was created by the following people: Lars Flöer, Nadine Giese, Russell 
Jurek, Martin Meyer, Attila Popping, Paolo Serra, Tobias Westmeier, and 
Benjamin Winkel.

© 2016 The SoFiA Authors

This programme is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.

This programme is distributed in the hope that it will be useful, but **without 
any warranty**; without even the implied warranty of **merchantability** or **fitness 
for a particular purpose**. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this programme. If not, see http://www.gnu.org/licenses/.
