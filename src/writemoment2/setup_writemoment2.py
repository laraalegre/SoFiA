# usage: python setup_writemoment2.py build_ext --inplace

import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

setup(
	name = 'New writemoment module',
	ext_modules = cythonize("writemoment2.pyx"),
	include_dirs = [numpy.get_include()],
)
