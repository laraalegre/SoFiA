# Define SoFiA version number
__version__ = "1.2.0-beta"
__version_full__ = "SoFiA " + __version__



# Fix Astropy 'clobber' issue
import inspect
from astropy.io.fits import writeto
if "clobber" in inspect.getargspec(writeto).args:
	__astropy_arg_overwrite__ = "clobber"
else:
	__astropy_arg_overwrite__ = "overwrite"
