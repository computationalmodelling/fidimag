import sys
import subprocess

cmd=('python',
     'setup.py',
     'build_ext',
     '--inplace')

try:
    subprocess.check_output(cmd, stderr=subprocess.STDOUT)
except subprocess.CalledProcessError, ex:
    sys.stderr.write(ex.output)
    raise Exception("make_modules: Make failed")

from sim import Sim
from fd_mesh import FDMesh
from exchange import UniformExchange
from anisotropy import Anisotropy
from zeeman import Zeeman
from demag import Demag
from materials import Nickel
from show_vector import VisualSpin