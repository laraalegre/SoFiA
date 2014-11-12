from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'New writemoment module',
  ext_modules = cythonize("writemoment2.pyx"),
)
