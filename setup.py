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


ABS_MODULE_DIR = Path(__file__).parent  # This should be the abs path
MODULE_DIR = Path('.')  # This should be the rel path (setup.py requires rel paths for pip -e)
print(MODULE_DIR)
INCLUDE_DIR = ABS_MODULE_DIR / 'local/include'
# Use absolute paths for the sundials/fftw libs:
LIB_DIR = ABS_MODULE_DIR / 'local/lib'
LIB_DIR64 = ABS_MODULE_DIR / 'local/lib64'

# rpath: run-time search path for the sundials (cvode) and fftw library objects
com_link = ['-Wl,-rpath,{},-rpath,{}'.format(str(LIB_DIR), str(LIB_DIR64)), '-fopenmp']
lib_paths = [str(LIB_DIR), str(LIB_DIR64)]
com_libs = ['m', 'fftw3_omp', 'fftw3', 'sundials_core', 'sundials_cvodes', 'sundials_nvecserial', 'sundials_nvecopenmp', 'blas', 'lapack']
com_args = ['-O3', '-Wno-cpp', '-Wno-unused-function', '-Wall', '-std=c99', '-fopenmp']
com_args_cpp = ['-O3', '-Wno-unused-function', '-Wall', '-std=c++14', '-fopenmp']

# Find all .pyx files with extensions (source files) -> relative paths
ROOT_DIR = MODULE_DIR / 'fidimag'
source_files = [s for s in ROOT_DIR.rglob('*.pyx')]  # Paths

# User extensions are located in the "user" namespace within "extensions"
ext_names = []
for s in source_files:
    if 'user' in str(s):
        ext_names.append("fidimag.extensions.user." + s.stem)
    else:
        ext_names.append("fidimag.extensions." + s.stem)

com_inc = [numpy.get_include(), str(INCLUDE_DIR)]

if 'SUNDIALS_INC' in os.environ:
     com_inc.append(os.environ['SUNDIALS_INC'])

if 'FFTW_INC' in os.environ:
     com_inc.append(os.environ['FFTW_INC'])

ext_modules = []
for i, (module, src) in enumerate(zip(ext_names, source_files)):
    print(sYellow + f"Compiling module {module}" + sReset)

    if 'fmmlib' in module:
        continue

    # "python -m build ..." can use absolute paths
    # srcFiles = [str(sF.resolve()) for sF in src.parent.glob('*')  # resolve -> absolute paths
    #             if sF.is_file()
    #             and sF != src.with_suffix('.c') 
    #             and str(sF).endswith(('.c', '.cpp'))
    #             ]

    # src is a Path
    srcFiles = [str(sF) for sF in src.parent.glob('*')
                if sF.is_file()
                and sF != src.with_suffix('.c')
                and sF != src.with_suffix('.cpp')
                and str(sF).endswith(('.c', '.cpp', '.pyx'))
                ]


    if 'fmm' in module:
        print(sBlue + f'Using cpp for this module: {module}' + sReset)
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


def get_version():
    with open('fidimag/__init__.py') as f:
        for line in f:
            m = re.match(r'''__version__\s*=\s*(['"])(.+)\1''', line.strip())
            if m:
                return m.group(2)
    raise Exception("Couldn't find __version__ in %s" % pkg_init_path)


nthreads = 0  # Disabled parallel compilation due to Python 3.14 multiprocessing issues (0 = no multiprocessing)
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
