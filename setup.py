from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import os
import glob
import re


if 'CC' in os.environ:
    print("Using CC={}".format(os.environ['CC']))
else:
    os.environ["CC"] = "gcc"
    print("Using CC={} (set by setup.py)".format(os.environ['CC']))

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(MODULE_DIR, "fidimag")
SUNDIALS_DIR = os.path.join(SRC_DIR, "common", "sundials")
NEB_DIR = os.path.join(SRC_DIR, "common", "neb")
ATOM_DIR = os.path.join(SRC_DIR, "atomistic", "lib")
MICRO_DIR = os.path.join(SRC_DIR, "micro", "lib")
BARYAKHTAR_DIR = os.path.join(MICRO_DIR, "baryakhtar")
DEMAG_DIR = os.path.join(SRC_DIR, "common", "dipolar")

LOCAL_DIR = os.path.join(MODULE_DIR, "local")
INCLUDE_DIR = os.path.join(LOCAL_DIR, "include")
LIB_DIR = os.path.join(LOCAL_DIR, "lib")
print("LIB_DIR={}".format(LIB_DIR))


pkg_init_path = os.path.join(
    os.path.dirname(__file__), 'fidimag', '__init__.py')


def get_version():
    with open(pkg_init_path) as f:
        for line in f:
            m = re.match(r'''__version__\s*=\s*(['"])(.+)\1''', line.strip())
            if m:
                return m.group(2)
    raise Exception("Couldn't find __version__ in %s" % pkg_init_path)

version = get_version()


def glob_cfiles(path, excludes):
    cfiles = []
    for cfile in glob.glob(os.path.join(path, "*.c")):
        filename = os.path.basename(cfile)
        if not filename in tuple(excludes):
            cfiles.append(cfile)
    return cfiles

sources = []
sources.append(os.path.join(ATOM_DIR, 'clib.pyx'))
sources += glob_cfiles(ATOM_DIR, excludes=["clib.c"])

cvode_sources = []
cvode_sources.append(os.path.join(SUNDIALS_DIR, 'cvode.pyx'))
cvode_sources += glob_cfiles(SUNDIALS_DIR, excludes=["cvode.c"])

baryakhtar_sources = []
baryakhtar_sources.append(os.path.join(BARYAKHTAR_DIR, 'baryakhtar_clib.pyx'))
baryakhtar_sources += glob_cfiles(BARYAKHTAR_DIR,
                                  excludes=["baryakhtar_clib.c"])

micro_sources = []
micro_sources.append(os.path.join(MICRO_DIR, 'micro_clib.pyx'))
micro_sources += glob_cfiles(MICRO_DIR, excludes=["micro_clib.c"])

neb_sources = []
neb_sources.append(os.path.join(NEB_DIR, 'neb_clib.pyx'))
neb_sources += glob_cfiles(NEB_DIR, excludes=["neb_clib.c"])

dipolar_sources = []
dipolar_sources.append(os.path.join(DEMAG_DIR, 'dipolar.pyx'))
dipolar_sources += glob_cfiles(DEMAG_DIR, excludes=["dipolar.c"])


com_libs = ['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes',
            'sundials_nvecserial', 'sundials_nvecopenmp']

com_args = ['-std=c99']
com_link = ['-L%s' % LIB_DIR]

if 'icc' in os.environ['CC']:
    com_args.append('-openmp')
    com_link.append('-openmp')
else:
    com_args.append('-fopenmp')
    com_link.append('-fopenmp')


com_inc = [numpy.get_include(), INCLUDE_DIR]

ext_modules = [
    Extension("fidimag.extensions.clib",
              sources=sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.cvode",
              sources=cvode_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.micro_clib",
              sources=micro_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.neb_clib",
              sources=neb_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.cvode",
              sources=cvode_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.micro_clib",
              sources=micro_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.neb_clib",
              sources=neb_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.dipolar",
              sources=dipolar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
]

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
    ext_modules=cythonize(ext_modules),
)
