import sys


import subprocess
import os
realpath=os.path.realpath(__file__)
pccp_path=os.path.split(os.path.split(realpath)[0])[0]

cmd=('python',
     os.path.join(pccp_path,'setup.py'),
     'build_ext',
     '--inplace')

FNULL = open(os.devnull, 'w')
try:
    subprocess.check_call(cmd, stdout=FNULL, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError, ex:
    sys.stderr.write(ex.output)
    raise Exception("make_modules: Make failed")


from sim import Sim
from mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from zeeman import TimeZeeman
from demag import Demag
from dmi import DMI
from fileio import DataSaver, DataReader
from materials import UnitMaterial
#from materials import Nickel

#from show_vector import VisualSpin
