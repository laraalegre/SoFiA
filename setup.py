#!/usr/bin/python
# -*- coding: utf-8 -*-
from distutils.core import setup, Extension
from distutils.version import StrictVersion, LooseVersion
import distutils.command.build

import sys
import os
import stat
import warnings

from sofia import __version__ as version

# compile the gui?
compile_gui = True

# Dependency checking
dependencies = [['numpy', '1.8'], ['scipy', None]]

for (pkg, minversion) in dependencies:
    try:
        m = __import__(pkg)
        if minversion is not None:
            if StrictVersion(m.__version__) < StrictVersion(minversion):
                if LooseVersion(m.__version__) < LooseVersion(minversion):
                    raise ValueError
                warnings.warn(
                    'Version', m.__version__,
                    'of package', pkg,
                    'might not be sufficient')
    except ImportError:
        print 'Package', pkg, 'not present.'
        sys.exit(1)
    except ValueError:
        print 'Package', pkg, 'has version', m.__version__
        print 'Version', minversion, 'required.'
        sys.exit(1)

class build(distutils.command.build.build):

    user_options = distutils.command.build.build.user_options + [
        ('no-gui=', None,
        'specify if you want the gui to be built'),
        ]

    def initialize_options(self, *args, **kwargs):
        self.no_gui = False
        distutils.command.build.build.initialize_options(self, *args, **kwargs)

    def run(self, *args, **kwargs):
        global nogui
        nogui = self.no_gui
        distutils.command.build.build.run(self, *args, **kwargs)



import numpy

# Define include directories
include_dirs = []
include_dirs.append('src')
include_dirs.append('src/linker')
include_dirs.append(numpy.get_include())

# Get external include directories from PATH
ext_include_dirs = []
ext_library_dirs = []
try:
    ext_include_dirs += os.environ['PATH'].split(os.pathsep)
except:
    print 'Could not parse PATH'
    sys.exit(1)

# Get external library paths from LD_LIBRARY_PATH, if it exists
try:
    ext_library_dirs += os.environ['LD_LIBRARY_PATH'].split(os.pathsep)
except KeyError:
    pass

if sys.platform == 'darwin':
    # Remove "bin" directories to avoid problems with clang
    ext_include_dirs = [ext_include_dir for ext_include_dir in ext_include_dirs if not 'bin'
                        in ext_include_dir]

# C/C++ source code files
# Object linking code
linker_src_base = 'src/linker/'
linker_src_files = [
    'linker.cpp',
    'RJJ_ObjGen_CreateObjs.cpp',
    'RJJ_ObjGen_DetectDefn.cpp',
    'RJJ_ObjGen_MemManage.cpp',
    'RJJ_ObjGen_ThreshObjs.cpp',
    'RJJ_ObjGen_Dmetric.cpp'
    ]
linker_src = [linker_src_base + f for f in linker_src_files]

# 2D1D wavelet finder code
wavelet_src_base = 'src/wavelet/'
wavelet_src_files = [
    'wavelet.c'
    ]
wavelet_src = [wavelet_src_base + f for f in wavelet_src_files]

# CNHI finder code
CNHI_src_base = 'src/CNHI/'
CNHI_src_files = [
    'CNHI.cpp'
    ]
CNHI_src = [CNHI_src_base + f for f in CNHI_src_files]

# moment output module
writemoment2_src_base = 'src/writemoment2/'
writemoment2_src_files = [
    'writemoment2.c'
    ]
writemoment2_src = [writemoment2_src_base + f for f in writemoment2_src_files]

# Interface to the parametrization code
parametrizer_src_base = 'src/parametrizer/'
parametrizer_src_files = [
    'BFfit.cpp',
    'cparametrizer.cpp',
    'DataCube.cpp',
    'helperFunctions.cpp',
    'MaskOptimization.cpp',
    'Measurement.cpp',
    'MetaData.cpp',
    'ModuleParametrisation.cpp',
    'Parametrization.cpp',
    'Source.cpp',
    'SourceCatalog.cpp',
    'Unit.cpp',
    ]
parametrizer_src = [parametrizer_src_base + f for f in parametrizer_src_files]

setup(
    name='sofia',
    version=version,
    ext_package='sofia',
    ext_modules=[
        Extension(
            'linker',
            linker_src,
            extra_compile_args=['-O3'],
            include_dirs=include_dirs),
        Extension(
            'wavelet',
            wavelet_src,
            extra_compile_args=['-O3'],
            include_dirs=include_dirs),
        Extension(
            'CNHI',
            CNHI_src,
            extra_compile_args=['-I. -O3'],
            include_dirs=include_dirs),
        Extension(
            'writemoment2',
            writemoment2_src,
            extra_compile_args=['-O3'],
            include_dirs=include_dirs),
        Extension(
            'cparametrizer',
            parametrizer_src,
            extra_compile_args=['-O3'],
            include_dirs=include_dirs + ext_include_dirs,
            library_dirs=ext_library_dirs + ext_include_dirs#,
            #libraries=['gsl', 'gslcblas']
            )
        ],
    package_dir={'sofia': 'sofia'},
    packages=['sofia'],  #,"cfitsio", "wcs"
    cmdclass={
        'build': build,
    },
    )

# make sofia_pipeline.py executable
os.chmod('sofia_pipeline.py', os.stat('sofia_pipeline.py').st_mode | stat.S_IXUSR)

# path to the sofia modules
cwd = os.getcwd()
for ll in os.listdir(cwd + '/build'):
    if 'lib.' in ll:
        sofiaModulesPath = cwd + '/build/' + ll

# set system variable
os.environ['SOFIA_PIPELINE_PATH'] = cwd + '/sofia_pipeline.py'
print cwd + '/sofia_pipeline.py'
print os.environ['SOFIA_PIPELINE_PATH']

# compile the SoFiA gui
if nogui=='True': compile_gui = False
else: compile_gui = True

if compile_gui:
    os.chdir('gui')
    os.system('qmake; make')
    if sys.platform == 'darwin':
        os.chmod(
            'SoFiA.app/Contents/MacOS/SoFiA',
            os.stat('SoFiA.app/Contents/MacOS/SoFiA').st_mode | stat.S_IXUSR
            )
        sofiaApplicationPath = cwd + '/gui/SoFiA.app/Contents/MacOS'
    else:
        os.chmod(
            'SoFiA',
            os.stat('SoFiA').st_mode | stat.S_IXUSR
            )
        sofiaApplicationPath = cwd + '/gui'
    os.chdir('../')

print '\n'
print '-------------------------------------------------------------------------'
print '\n'
print '\033[1;32mInstallation complete.\033[0m\n'
print '\033[1mPlease add the following lines to your shell configuration file:\033[0m\n'
print '\033[3;4mFor BASH (~/.bashrc):\033[0m'
print '    export SOFIA_MODULE_PATH="' + sofiaModulesPath + '"'
print '    export SOFIA_PIPELINE_PATH="' + cwd + '/sofia_pipeline.py"'
if compile_gui: print '    export PATH="$PATH:' + cwd + ':' + sofiaApplicationPath + '"\n'
else:           print '    export PATH="$PATH:' + cwd + '"\n'
print '\033[3;4mFor (T)CSH (~/.cshrc):\033[0m'
print '    setenv SOFIA_MODULE_PATH "' + sofiaModulesPath + '"'
print '    setenv SOFIA_PIPELINE_PATH "' + cwd + '/sofia_pipeline.py"'
if compile_gui:  print '    setenv PATH {$PATH}:"' + cwd + ':' + sofiaApplicationPath + '"'
else:            print '    setenv PATH {$PATH}:"' + cwd + '"'
print '\n'

## test sofia installation
# print 'If you want to test your SoFiA installation, please open a new terminal and type:'
# print '  cd ' + cwd + '/test_data; sofia_pipeline.py parameters.txt'
# print '\n'

