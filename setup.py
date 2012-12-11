from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

#python setup.py build_ext --inplace

import fnmatch
import os

print __file__
print os.getcwd()
cp_path=os.path.join(os.getcwd(),'cp')

sources = []
sources.append(os.path.join(cp_path,'clib.pyx'))
for root, dirnames, filenames in os.walk(cp_path):
    for filename in fnmatch.filter(filenames, '*.c'):
        if filename!='clib.c':
            sources.append(os.path.join(root, filename))

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

