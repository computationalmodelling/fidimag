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
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
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

    # src is a Path
    srcFiles = [str(sF) for sF in src.parent.glob('*') if sF.is_file()
                and sF != src.with_suffix('.c') 
                and str(sF).endswith(('.c', '.cpp'))
                ]

    com_inc = [numpy.get_include(), INCLUDE_DIR] + [str(src.parent)]

    # print(module)
    # print(com_inc)
    # print(com_libs)
    # print(srcFiles)
    # print(lib_paths)
    # print(com_link)
    # print(com_args_compiler)
    for s in srcFiles:
        print(s)

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
                                 extra_link_args=com_link,
                                 extra_compile_args=com_args_compiler,
                                 language=lan
                                 )
    )


if 'CC' in os.environ:
    print("Using CC={}".format(os.environ['CC']))
else:
    os.environ["CC"] = "gcc"
    print("Using CC={} (set by setup.py)".format(os.environ['CC']))

# # setup.py requires relative paths:
# MODULE_DIR = os.path.dirname(os.path.relpath(__file__))
# SRC_DIR = os.path.join(MODULE_DIR, "fidimag")
# #
# SUNDIALS_DIR = os.path.join(SRC_DIR, "common", "sundials")
# NEB_DIR = os.path.join(SRC_DIR, "common", "neb")
# NEBM_DIR = os.path.join(SRC_DIR, "common", "neb_method")
# ATOM_DIR = os.path.join(SRC_DIR, "atomistic", "lib")
# COMMON_DIR = os.path.join(SRC_DIR, "common", "lib")
# MICRO_DIR = os.path.join(SRC_DIR, "micro", "lib")
# BARYAKHTAR_DIR = os.path.join(MICRO_DIR, "baryakhtar")
# DEMAG_DIR = os.path.join(SRC_DIR, "common", "dipolar")
# 
# FMMLIB_DIR = os.path.join(SRC_DIR, "atomistic", "fmmlib")
USER_DIR = os.path.join("fidimag/user")
# 
# LOCAL_DIR = os.path.join(MODULE_DIR, "local")
# INCLUDE_DIR = os.path.join(LOCAL_DIR, "include")
# LIB_DIR = os.path.join(LOCAL_DIR, "lib")
# LIB_DIR64 = os.path.join(LOCAL_DIR, "lib64")


pkg_init_path = os.path.join(
    os.path.dirname(__file__), 'fidimag', '__init__.py')


def get_version():
    with open(pkg_init_path) as f:
        for line in f:
            m = re.match(r'''__version__\s*=\s*(['"])(.+)\1''', line.strip())
            if m:
                return m.group(2)
    raise Exception("Couldn't find __version__ in %s" % pkg_init_path)


def glob_cfiles(path, excludes, extension="*.c"):
    cfiles = []
    for cfile in glob.glob(os.path.join(path, extension)):
        filename = os.path.basename(cfile)
        if filename not in tuple(excludes):
            cfiles.append(cfile)
    return cfiles
# 
# 
# sources = []
# sources.append(os.path.join(ATOM_DIR, 'clib.pyx'))
# sources += glob_cfiles(ATOM_DIR, excludes=["clib.c"])
#  
# common_sources = []
# common_sources.append(os.path.join(COMMON_DIR, 'common_clib.pyx'))
# common_sources += glob_cfiles(COMMON_DIR, excludes=["common_clib.c"])
# 
# cvode_sources = []
# cvode_sources.append(os.path.join(SUNDIALS_DIR, 'cvode.pyx'))
# cvode_sources += glob_cfiles(SUNDIALS_DIR, excludes=["cvode.c"])
# 
# baryakhtar_sources = []
# baryakhtar_sources.append(os.path.join(BARYAKHTAR_DIR, 'baryakhtar_clib.pyx'))
# baryakhtar_sources += glob_cfiles(BARYAKHTAR_DIR,
#                                   excludes=["baryakhtar_clib.c"])
# 
# micro_sources = []
# micro_sources.append(os.path.join(MICRO_DIR, 'micro_clib.pyx'))
# micro_sources += glob_cfiles(MICRO_DIR, excludes=["micro_clib.c"])
# 
# nebm_sources = []
# nebm_sources.append(os.path.join(NEBM_DIR, "nebm_clib.pyx"))
# nebm_sources += glob_cfiles(NEBM_DIR, excludes=["nebm_clib.c"])
# 
# dipolar_sources = []
# dipolar_sources.append(os.path.join(DEMAG_DIR, 'dipolar.pyx'))
# dipolar_sources += glob_cfiles(DEMAG_DIR, excludes=["dipolar.c"])
# 
# fmm_sources = []
# fmm_sources.append(os.path.join(FMMLIB_DIR, 'fmm.pyx'))
# fmm_sources += glob_cfiles(FMMLIB_DIR, excludes=["fmm.cpp"], extension="*.cpp")
# fmm_sources += glob_cfiles(FMMLIB_DIR, excludes=[], extension="*.c")
# 
# 
# 
# com_libs = ['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes',
#             'sundials_nvecserial', 'sundials_nvecopenmp', 'blas', 'lapack']
# 
# 
# com_args = ['-O3', '-Wno-cpp', '-Wno-unused-function', '-Wall']
# 
# 
# 
# # rpath is the path relative to the compiled shared object files (e.g. clib.so, etc)
# # which the dynamic linker looks for the linked libraries (e.g. libsundials_*.so) in.
# # We need to set it relatively in order for it to be preserved if the parent directory is moved
# # hence why it is a 'relative'(r) path. Here the relative path is with respect to
# # the fidimag/fidimag/extensions directory.
# RPATH = '../../local/lib'
# com_link = ['-Wl,-rpath,{},-rpath,{}'.format(LIB_DIR, LIB_DIR64)]
# lib_paths = [LIB_DIR, LIB_DIR64]
# 
# 
# com_args.append('-fopenmp')
# com_link.append('-fopenmp')
# 
# 
# com_inc = [numpy.get_include(), INCLUDE_DIR]
# 
# if 'SUNDIALS_DIR' in os.environ:
#     lib_paths.append(os.environ['SUNDIALS_DIR'])
#     com_inc.append(os.environ['SUNDIALS_INC'])
# 
# if 'FFTW_DIR' in os.environ:
#     lib_paths.append(os.environ['FFTW_DIR'])
#     com_inc.append(os.environ['FFTW_INC'])
# 
# com_args_cpp = com_args.copy()
# com_args_cpp.append('-std=c++14')
# 
# com_args.append('-std=c99')
# com_inc_cpp = com_inc.copy()
# com_link_cpp = com_link.copy()
# com_libs_cpp = com_libs.copy()
# lib_paths_cpp = lib_paths.copy()
# 
# ext_modules = [
#     Extension("fidimag.extensions.clib",
#               sources=sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.common_clib",
#               sources=common_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.cvode",
#               sources=cvode_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.baryakhtar_clib",
#               sources=baryakhtar_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.micro_clib",
#               sources=micro_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.nebm_clib",
#               sources=nebm_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.cvode",
#               sources=cvode_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.baryakhtar_clib",
#               sources=baryakhtar_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.micro_clib",
#               sources=micro_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.dipolar",
#               sources=dipolar_sources,
#               include_dirs=com_inc,
#               libraries=com_libs,
#               library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#               extra_compile_args=com_args,
#               extra_link_args=com_link,
#               ),
#     Extension("fidimag.extensions.fmm",
#               sources=fmm_sources,
#               include_dirs=com_inc_cpp,
#               libraries=com_libs_cpp,
#               library_dirs=lib_paths_cpp, runtime_library_dirs=lib_paths_cpp,
#               extra_compile_args=com_args_cpp,
#               extra_link_args=com_link_cpp,
#               language="c++",
#               )
# ]


# for folder in glob.glob(os.path.join(USER_DIR, '*/')):
#     module_name = folder.split('/')[-2]
#     print(sYellow + f'Found User Module: {module_name}' + sReset)
#     user_sources = glob.glob(folder + '/*.pyx')
#     if len(user_sources) != 1:
#         raise BuildError("User Modules are only allowed one Cython .pyx file")
# 
#     filename_string = user_sources[0].split('/')[-1][:-4]
#     if filename_string != module_name:
#         print(filename_string, module_name)
#         raise BuildError("The Cython source file in {} must match the folder name - i.e. it must be {}.pyx".format(module_name, module_name))
#     cfilename = filename_string + '.c'
#     user_sources += glob_cfiles(folder, excludes=[cfilename])
# 
# 
#     ext_modules.append(
#        Extension("fidimag.extensions.user.{}".format(module_name),
#           sources=user_sources,
#           include_dirs=com_inc,
#           libraries=com_libs,
#           library_dirs=lib_paths, runtime_library_dirs=lib_paths,
#           extra_compile_args=com_args,
#           extra_link_args=com_link,
#        ),
#     )

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
