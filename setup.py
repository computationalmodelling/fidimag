from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

#python setup.py build_ext --inplace



ext_modules = [
    Extension("clib",
              sources = ['demag.c','llg.c','exch.c','anis.c','clib.pyx'],
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

