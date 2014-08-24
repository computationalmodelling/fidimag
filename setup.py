from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

#python setup.py build_ext --inplace

import fnmatch
import os
import glob

#print __file__
#print os.getcwd()
realpath=os.path.realpath(__file__)
pccp_path=os.path.split(realpath)[0]

cp_path=os.path.join(pccp_path,'cp')
os.chdir(pccp_path)

sources = []
sources.append(os.path.join(cp_path,'clib.pyx'))


cfiles = glob.glob(os.path.join(cp_path,'*.c'))
for cf in cfiles:
    if not cf.endswith("clib.c"):
        sources.append(cf)

print 'sources',sources


sundials_path = os.path.join(cp_path,'sundials')
cvode_sources = []
cvode_sources.append(os.path.join(sundials_path,'cvode.pyx'))
cfiles = glob.glob(os.path.join(sundials_path,'*.c'))
for cf in cfiles:
    if not cf.endswith("cvode.c"):
        cvode_sources.append(cf)


baryakhtar_path=os.path.join(pccp_path,'baryakhtar')
baryakhtar_path = os.path.join(baryakhtar_path,'cython')

baryakhtar_sources = []
baryakhtar_sources.append(os.path.join(baryakhtar_path,'baryakhtar_clib.pyx'))
cfiles = glob.glob(os.path.join(baryakhtar_path,'*.c'))
for cf in cfiles:
    if not cf.endswith("baryakhtar_clib.c"):
        baryakhtar_sources.append(cf)


libs_path = os.path.join(pccp_path,'libs')
include_path = os.path.join(libs_path,'include')
lib_path = os.path.join(libs_path,'lib')


ext_modules = [
    Extension("clib",
              sources = sources,
              include_dirs = [numpy.get_include(),include_path],
              libraries=['m','fftw3','sundials_cvodes','sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s'%lib_path,'-fopenmp'],
        ),
    Extension("cvode",
              sources = cvode_sources,
              include_dirs = [numpy.get_include(),include_path],
              libraries=['m','fftw3','sundials_cvodes','sundials_nvecserial'],
              extra_compile_args=["-fopenmp"],
              extra_link_args=['-L%s'%lib_path,'-fopenmp'],
              #extra_link_args=["-g"],
        ),
               
    Extension("baryakhtar_clib",
              sources = baryakhtar_sources,
              include_dirs = [numpy.get_include(),include_path],
              libraries=['m','fftw3','sundials_cvodes','sundials_nvecserial'],
              extra_compile_args=["-fopenmp", '-std=c99'],
              extra_link_args=['-L%s'%lib_path,'-fopenmp'],
              #extra_link_args=["-g"],
        )
    ]
    

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
