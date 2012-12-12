from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

#python setup.py build_ext --inplace

import fnmatch
import os

#print __file__
#print os.getcwd()
realpath=os.path.realpath(__file__)
pccp_path=os.path.split(realpath)[0]

cp_path=os.path.join(pccp_path,'cp')
os.chdir(pccp_path)

sources = []
sources.append(os.path.join('cp','clib.pyx'))
for root, dirnames, filenames in os.walk(cp_path):
    for filename in fnmatch.filter(filenames, '*.c'):
        if filename!='clib.c':
            sources.append(os.path.join('cp',filename))

print sources
ext_modules = [
    Extension("clib",
              sources = sources,
              include_dirs = [numpy.get_include()],
              libraries=['m','fftw3'],
              extra_compile_args=["-fopenmp"],
              #extra_link_args=["-g"],
        )
    ]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
