"""
import sys
import subprocess
import os
realpath=os.path.realpath(__file__)
pccp_path=os.path.split(realpath)[0]

cmd=('python',
     os.path.join(pccp_path,'setup.py'),
     'build_ext',
     '--inplace')

try:
    FNULL = open(os.devnull, 'w')
    subprocess.check_call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError, ex:
    sys.stderr.write(ex.output)
    raise Exception("make_modules: Make failed")
"""
