"""
Creates working dylibs on Mac OS X by adding the full path to referenced libraries.
This fixes ImportErrors that may appear starting with SUNDIALS 2.6.

"""
import os
import glob
import re
import subprocess

THIS_FILE = os.path.dirname(os.path.abspath(__file__))
MODULE_DIR = os.path.abspath(os.path.join(THIS_FILE, os.pardir))
LIB_DIR = os.path.join(MODULE_DIR, "local", "lib")

EXTENSION_DIR = os.path.join(MODULE_DIR, "fidimag", "extensions")

sos = glob.glob(os.path.join(EXTENSION_DIR, "*.so"))

patten = re.compile(r'libsundials[\w.]+\.dylib')
def extract_library(so_file):
    cmd = ('otool', '-L', so_file)
    output = subprocess.check_output(cmd)
    for line in output.decode().split('\t'):
        m = patten.match(line)
        print(m, line)
        if m:
            lib_name =  m.group()
            full_name = os.path.join(LIB_DIR, lib_name)
            print(lib_name, full_name)
            cmd = ('install_name_tool', '-change',
                 lib_name, full_name, so_file
                 )
            os.system(' '.join(cmd))

for so in sos:
    extract_library(so)
    print("Done for %s!"%so)
