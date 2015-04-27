"""
You shouldn't call this script directly but rather use the Makefile as in
`make build`. If you decide to use this script directly, then create the
file __init__.py in the fidimag/extensions directory.

"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import fnmatch
import os
import glob

os.environ["CC"] = "gcc"

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(MODULE_DIR, "fidimag")
CP_DIR = os.path.join(SRC_DIR, "cp")
SUNDIALS_DIR = os.path.join(CP_DIR, "sundials")
BARYAKHTAR_DIR = os.path.join(SRC_DIR, "baryakhtar", "cython")
MICRO_DIR = os.path.join(SRC_DIR, "micro", "cython")
NEB_DIR = os.path.join(CP_DIR, "neb")
LOCAL_DIR = os.path.join(MODULE_DIR, "local")
INCLUDE_DIR = os.path.join(LOCAL_DIR, "include")
LIB_DIR = os.path.join(LOCAL_DIR, "lib")


def glob_cfiles(path, excludes):
    cfiles = []
    for cfile in glob.glob(os.path.join(path, "*.c")):
        if not cfile.endswith(tuple(excludes)):
            cfiles.append(cfile)
    return cfiles

sources = []
sources.append(os.path.join(CP_DIR, 'clib.pyx'))
sources += glob_cfiles(CP_DIR, excludes=["clib.c"])

cvode_sources = []
cvode_sources.append(os.path.join(SUNDIALS_DIR, 'cvode.pyx'))
cvode_sources += glob_cfiles(SUNDIALS_DIR, excludes=["cvode.c"])

baryakhtar_sources = []
baryakhtar_sources.append(os.path.join(BARYAKHTAR_DIR, 'baryakhtar_clib.pyx'))
baryakhtar_sources += glob_cfiles(BARYAKHTAR_DIR, excludes=["baryakhtar_clib.c"])

micro_sources = []
micro_sources.append(os.path.join(MICRO_DIR, 'micro_clib.pyx'))
micro_sources += glob_cfiles(MICRO_DIR, excludes=["micro_clib.c"])

neb_sources = []
neb_sources.append(os.path.join(NEB_DIR, 'neb_clib.pyx'))
neb_sources += glob_cfiles(NEB_DIR, excludes=["neb_clib.c"])

ext_modules = [
    Extension("clib",
              sources=sources,
              include_dirs=[numpy.get_include(), INCLUDE_DIR],
              libraries=['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes', 'sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s' % LIB_DIR, '-fopenmp'],
              ),
    Extension("cvode",
              sources=cvode_sources,
              include_dirs=[numpy.get_include(), INCLUDE_DIR],
              libraries=['m', 'fftw3', 'sundials_cvodes', 'sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s' % LIB_DIR, '-fopenmp'],
              ),
    Extension("baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=[numpy.get_include(), INCLUDE_DIR],
              libraries=['m', 'fftw3', 'sundials_cvodes', 'sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s' % LIB_DIR, '-fopenmp'],
              ),
    Extension("micro_clib",
              sources=micro_sources,
              include_dirs=[numpy.get_include(), INCLUDE_DIR],
              libraries=['m', 'fftw3', 'sundials_cvodes', 'sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s' % LIB_DIR, '-fopenmp'],
              ),
    Extension("neb_clib",
              sources=neb_sources,
              include_dirs=[numpy.get_include(), INCLUDE_DIR],
              libraries=['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes', 'sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s' % LIB_DIR, '-fopenmp'],
              ),
]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)
