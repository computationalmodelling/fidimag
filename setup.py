from setuptools import setup
from setuptools import Extension
import multiprocessing
from Cython.Build import cythonize
import numpy
import os
import glob
import re
from pathlib import Path

sYellow = "\x1b[33;49m"
sBlue = "\x1b[34;49m"
sRed = "\x1b[31;49m"
sReset = "\x1b[0m"


class BuildError(Exception):
    pass


# setup.py requires relative paths:
# MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
MODULE_DIR = os.path.dirname(os.path.relpath(__file__))
INCLUDE_DIR = os.path.join(MODULE_DIR, 'local', 'include')
LIB_DIR = os.path.join(MODULE_DIR, 'local', 'lib')
LIB_DIR64 = os.path.join(MODULE_DIR, 'local', 'lib64')
# INCLUDE_DIR = 'local/include'
# LIB_DIR = 'local/lib'
# LIB_DIR64 = 'local/lib64'

# rpath is the path relative to the compiled shared object files (e.g. clib.so, etc)
# which the dynamic linker looks for the linked libraries (e.g. libsundials_*.so) in.
# We need to set it relatively in order for it to be preserved if the parent directory is moved
# hence why it is a 'relative'(r) path. Here the relative path is with respect to
# the fidimag/fidimag/extensions directory.
RPATH = '../../local/lib'
com_link = ['-Wl,-rpath,{},-rpath,{}'.format(str(LIB_DIR), str(LIB_DIR64)), '-fopenmp']

lib_paths = [LIB_DIR, LIB_DIR64]
# lib_paths = [LIB_DIR, os.path.join(MODULE_DIR, 'native')]
# com_inc = [numpy.get_include(), INCLUDE_DIR, os.path.join(MODULE_DIR, 'native', 'include')]
# com_libs = ['fidimag']
com_libs = ['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes', 'sundials_nvecserial', 'sundials_nvecopenmp', 'blas', 'lapack']
# com_args = []
com_args = ['-O3', '-Wno-cpp', '-Wno-unused-function', '-Wall', '-std=c99', '-fopenmp']
com_args_cpp = ['-O3', '-Wno-unused-function', '-Wall', '-std=c++14', '-fopenmp']

if 'SUNDIALS_INC' in os.environ:
     com_inc.append(os.environ['SUNDIALS_INC'])

if 'FFTW_INC' in os.environ:
     com_inc.append(os.environ['FFTW_INC'])

# Find .pyx files with extensions (source files)
ROOT_DIR = Path('fidimag')
source_files = [s for s in ROOT_DIR.rglob('*.pyx')]  # Paths
print(source_files)
ext_names = []
for s in source_files:
    if 'user' in str(s):
        ext_names.append("fidimag.extensions.user." + s.stem)
    else:
        ext_names.append("fidimag.extensions." + s.stem)
print(ext_names)

ext_modules = []
for i, (module, src) in enumerate(zip(ext_names, source_files)):
    print(sYellow + f"Compiling module {module}" + sReset)

    if 'fmmlib' in module:
        continue

    # src is a Path
    # srcFiles = [str(sF.resolve()) for sF in src.parent.glob('*')  # resolve -> absolute paths
    #             if sF.is_file()
    #             and sF != src.with_suffix('.c') 
    #             and str(sF).endswith(('.c', '.cpp'))
    #             ]

    srcFiles = [str(sF) for sF in src.parent.glob('*')
                if sF.is_file()
                and sF != src.with_suffix('.c')
                and sF != src.with_suffix('.cpp')
                and str(sF).endswith(('.c', '.cpp', '.pyx'))
                ]

    com_inc = [numpy.get_include(), INCLUDE_DIR]

    print(module)
    print(com_inc)
    print(com_libs)
    print(srcFiles)
    print(lib_paths)
    print(com_link)
    # print(com_args_compiler)
    for s in srcFiles:
        print(s)
    print(com_inc)

    if 'fmm' in module:
        print(sBlue + f'Using cpp for this module' + sReset)
        com_args_compiler = com_args_cpp
        lan = 'c++'
    else:
        com_args_compiler = com_args
        lan = 'c'

    ext_modules.append(Extension(module,
                                 sources=srcFiles,
                                 include_dirs=com_inc,
                                 libraries=com_libs,
                                 library_dirs=lib_paths,
                                 runtime_library_dirs=lib_paths,
                                 extra_compile_args=com_args_compiler,
                                 extra_link_args=com_link,
                                 language=lan
                                 )
    )


if 'CC' in os.environ:
    print("Using CC={}".format(os.environ['CC']))
else:
    os.environ["CC"] = "gcc"
    print("Using CC={} (set by setup.py)".format(os.environ['CC']))

USER_DIR = os.path.join("fidimag/user")


pkg_init_path = os.path.join(
    os.path.dirname(__file__), 'fidimag', '__init__.py')


def get_version():
    with open(pkg_init_path) as f:
        for line in f:
            m = re.match(r'''__version__\s*=\s*(['"])(.+)\1''', line.strip())
            if m:
                return m.group(2)
    raise Exception("Couldn't find __version__ in %s" % pkg_init_path)


nthreads = multiprocessing.cpu_count()
print(sYellow + f'Building with {nthreads} threads' + sReset)
setup(
    name='fidimag',
    version=get_version(),
    packages=['fidimag',
              'fidimag.atomistic',
              'fidimag.micro',
              'fidimag.extensions',
              'fidimag.common',
              ],
    ext_modules=cythonize(ext_modules,
                          nthreads=nthreads,
                          compiler_directives={'linetrace': True, 'language_level' : "3"}
                          ),
)
