import os
import glob
import re
import subprocess

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
LIB_DIR = os.path.join(MODULE_DIR, "local", "lib")

EXTENSION_DIR = os.path.join(MODULE_DIR, "fidimag", "extensions")

sos = glob.glob(os.path.join(EXTENSION_DIR, "*.so"))

patten = re.compile(r'libsundials[\w.]+\.dylib')
def extract_library(so_file):
    cmd = ('otool', '-L', so_file)
    output = subprocess.check_output(cmd)
    for line in output.split('\t'):
        m = patten.match(line)
        if m: 
            lib_name =  m.group()
            full_name = os.path.join(LIB_DIR, lib_name)
            #print lib_name, full_name
            cmd = ('install_name_tool', '-change', 
                 lib_name, full_name, so_file
                 )
            os.system(' '.join(cmd))

for so in sos:
    extract_library(so)
    print("Done for %s!"%so)