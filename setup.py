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

#print 'sources',sources
sundials_path = os.path.join(cp_path,'sundials')

sources2 = []
sources2.append(os.path.join(sundials_path,'cvode.pyx'))
for root, dirnames, filenames in os.walk(sundials_path):
    for filename in fnmatch.filter(filenames, '*.c'):
        if filename!='cvode.c':
            sources.append(os.path.join(sundials_path,filename))
        print filename

print sources2


libs_path = os.path.join(pccp_path,'libs')
include_path = os.path.join(libs_path,'include')
lib_path = os.path.join(libs_path,'lib')

print include_path
print lib_path

ext_modules = [
    Extension("clib",
              sources = sources,
              include_dirs = [numpy.get_include(),include_path],
              libraries=['m','fftw3','sundials_cvodes','sundials_nvecserial'],
              extra_compile_args=["-fopenmp"],
              extra_link_args=['-L%s'%lib_path,'-fopenmp'],
        ),
    Extension("cvode",
              sources = sources2,
              include_dirs = [numpy.get_include(),include_path],
              libraries=['m','fftw3','sundials_cvodes','sundials_nvecserial'],
              extra_compile_args=["-fopenmp"],
              extra_link_args=['-L%s'%lib_path,'-fopenmp'],
              #extra_link_args=["-g"],
        )
    ]
    

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
