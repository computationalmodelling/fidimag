from distutils.core import setup
from distutils.extension import Extension
import multiprocessing
from tools import glob_files, BuildError, get_user_module_sources
from Cython.Build import cythonize
import numpy
import glob
import os


print(os.environ['CC'], os.environ['CXX'])
version = '5.0alpha'

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
INCLUDE_DIR = os.path.join(MODULE_DIR, 'local', 'include')
LIB_DIR = os.path.join(MODULE_DIR, 'local', 'lib')

RPATH = '../../local/lib'
com_link = ['-Wl,-rpath,{}'.format(LIB_DIR)]

lib_paths = [LIB_DIR, os.path.join(MODULE_DIR, 'native')]
com_inc = [numpy.get_include(), INCLUDE_DIR, os.path.join(MODULE_DIR, 'native', 'include')]
com_libs = ['fidimag']
com_args = []

if 'SUNDIALS_INC' in os.environ:
     com_inc.append(os.environ['SUNDIALS_INC'])

if 'FFTW_DIR' in os.environ:
     com_inc.append(os.environ['FFTW_INC'])



source_files = glob.glob(os.path.join('fidimag', '**', '*.pyx'), recursive=True)
print(source_files)
ext_names = ["fidimag.extensions." + s.split('/')[-1].rstrip('.pyx') for s in source_files]

print(ext_names)


ext_modules = []
for i, (module, src) in enumerate(zip(ext_names, source_files)):
     print("Compiling module {}".format(module))
     ext_modules.append(Extension(module,
                                  sources=[src],
                                  include_dirs=com_inc,
                                  libraries=com_libs,
 	                          library_dirs=lib_paths, runtime_library_dirs=lib_paths,
                                  extra_compile_args=com_args,
                                  extra_link_args=com_link,
    )
)

                        
     
nthreads = multiprocessing.cpu_count()

setup(
    name='fidimag',
    version=version,
    description='Finite difference micromagnetic code',
    packages=['fidimag',
              'fidimag.atomistic',
              'fidimag.micro',
              'fidimag.extensions',
              'fidimag.common',
    ],
    ext_modules=cythonize(ext_modules,
                          nthreads=nthreads,
                          compiler_directives={
                              'linetrace': True,
                              'language_level': '3',
                          }
    ),
)
