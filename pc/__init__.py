import sys
import subprocess

import os
realpath=os.path.realpath(__file__)
pccp_path=os.path.split(os.path.split(realpath)[0])[0]

cmd=('python',
     os.path.join(pccp_path,'setup.py'),
     'build_ext',
     '--inplace')

try:
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError, ex:
    sys.stderr.write(ex.output)
    raise Exception("make_modules: Make failed")

from sim import Sim
from mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from demag import Demag
from dmi import DMI 
from materials import Nickel
from fileio import DataSaver, DataReader
#from show_vector import VisualSpin
