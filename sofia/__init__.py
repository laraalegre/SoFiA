#! /usr/bin/env python

# ===========================
# Define SoFiA version number
# ===========================
__version__ = "1.2.0-beta"
__version_full__ = "SoFiA " + __version__



# =============================
# Fix Astropy's 'clobber' issue
# =============================
# Usage: from sofia import __astropy_arg_overwrite__
#        hdu.writeto(filename, [other options], **__astropy_arg_overwrite__)
import inspect
from astropy.io.fits import writeto
if "clobber" in inspect.getargspec(writeto).args:
	__astropy_arg_overwrite__ = {"clobber" : True}
else:
	__astropy_arg_overwrite__ = {"overwrite" : True}
