from distutils.core import setup
from distutils.extension import Extension
import multiprocessing
from tools import glob_files, BuildError, get_user_module_sources
from Cython.Build import cythonize
import numpy
import glob
import os

version = '5.0alpha'

#if 'CC' in os.environ:
#    print("Using CC={}".format(os.environ['CC']))
#else:
#    os.environ["CC"] = "gcc"
#    print("Using CC={} (set by setup.py)".format(os.environ['CC']))

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(MODULE_DIR, "fidimag")
SUNDIALS_DIR = os.path.join(SRC_DIR, "common", "sundials")
NEB_DIR = os.path.join(SRC_DIR, "common", "neb")
NEBM_DIR = os.path.join(SRC_DIR, "common", "neb_method")
ATOM_DIR = os.path.join(SRC_DIR, "atomistic", "lib")
COMMON_DIR = os.path.join(SRC_DIR, "common", "lib")
MICRO_DIR = os.path.join(SRC_DIR, "micro", "lib")
BARYAKHTAR_DIR = os.path.join(MICRO_DIR, "baryakhtar")
DEMAG_DIR = os.path.join(SRC_DIR, "common", "dipolar")
USER_DIR = os.path.join(SRC_DIR, "user")
LOCAL_DIR = os.path.join(MODULE_DIR, "local")
INCLUDE_DIR = os.path.join(LOCAL_DIR, "include")
LIB_DIR = os.path.join(LOCAL_DIR, "lib")

sources = []
sources.append(os.path.join(ATOM_DIR, 'clib.pyx'))
sources += glob_files(ATOM_DIR, excludes=["clib.cpp"])

common_sources = []
common_sources.append(os.path.join(COMMON_DIR, 'common_clib.pyx'))
common_sources += glob_files(COMMON_DIR, excludes=["common_clib.cpp"])

cvode_sources = []
cvode_sources.append(os.path.join(SUNDIALS_DIR, 'cvode.pyx'))
cvode_sources += glob_files(SUNDIALS_DIR, excludes=["cvode.cpp"])

baryakhtar_sources = []
baryakhtar_sources.append(os.path.join(BARYAKHTAR_DIR, 'baryakhtar_clib.pyx'))
baryakhtar_sources += glob_files(BARYAKHTAR_DIR,
                                  excludes=["baryakhtar_clib.cpp"])

micro_sources = []
micro_sources.append(os.path.join(MICRO_DIR, 'micro_clib.pyx'))
micro_sources += glob_files(MICRO_DIR, excludes=["micro_clib.cpp"])

# NEB Method ------------------------------------------------------------------

nebm_sources = []
nebm_sources.append(os.path.join(NEBM_DIR, "nebm_clib.pyx"))
nebm_sources += glob_files(NEBM_DIR, excludes=["nebm_clib.cpp"])

# -----------------------------------------------------------------------------

dipolar_sources = []
dipolar_sources.append(os.path.join(DEMAG_DIR, 'dipolar.pyx'))
dipolar_sources += glob_files(DEMAG_DIR, excludes=["dipolar.cpp"])


com_libs = ['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes',
            'sundials_nvecserial', 'sundials_nvecopenmp', 'blas', 'lapack']

com_args = ['-O3', '-Wno-cpp', '-Wno-unused-function']
# rpath is the path relative to the compiled shared object files (e.g. clib.so, etc)
# which the dynamic linker looks for the linked libraries (e.g. libsundials_*.so) in.
# We need to set it relatively in order for it to be preserved if the parent directory is moved
# hence why it is a 'relative'(r) path. Here the relative path is with respect to
# the fidimag/fidimag/extensions directory.
RPATH = '../../local/lib'
com_link = ['-Wl,-rpath,{}'.format(LIB_DIR)]
lib_paths = [LIB_DIR]


if 'icc' in os.environ['CC']:
    com_args.append('-openmp')
    com_link.append('-openmp')
else:
    com_args.append('-fopenmp')
    com_link.append('-fopenmp')


com_inc = [numpy.get_include(), INCLUDE_DIR]

if 'SUNDIALS_DIR' in os.environ:
    lib_paths.append(os.environ['SUNDIALS_DIR'])
    com_inc.append(os.environ['SUNDIALS_INC'])

if 'FFTW_DIR' in os.environ:
    lib_paths.append(os.environ['FFTW_DIR'])
    com_inc.append(os.environ['FFTW_INC'])

ext_modules = [
    Extension("fidimag.extensions.clib",
              sources=sources,
              include_dirs=com_inc,
              libraries=com_libs,
	          library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.common_clib",
              sources=common_sources,
              include_dirs=com_inc,
              libraries=com_libs,
	          library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.cvode",
              sources=cvode_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.micro_clib",
              sources=micro_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.nebm_clib",
              sources=nebm_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.cvode",
              sources=cvode_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.micro_clib",
              sources=micro_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.dipolar",
              sources=dipolar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
]

for folder in glob.glob(os.path.join(USER_DIR, '*/')):
    module_name = folder.split('/')[-2]
    user_sources = get_user_module_sources(folder)

    ext_modules.append(
       Extension("fidimag.extensions.user.{}".format(module_name),
          sources=user_sources,
          include_dirs=com_inc,
          libraries=com_libs,
          library_dirs=lib_paths, runtime_library_dirs=lib_paths,
          extra_compile_args=com_args,
          extra_link_args=com_link,
       ),
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
