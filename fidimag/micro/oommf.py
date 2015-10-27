"""
This code reproduces a Fidimag simulation using OOMMF

"""

import os
import logging
import subprocess
import sys
import numpy as np

import omf

from fidimag.common import CuboidMesh


def setup_logger():

    logger = logging.getLogger("oommf")

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)
    logger.setLevel(logging.INFO)

    return logger

logger = setup_logger()

# OOMMF bash variables used by Fidimag to locate OOMMF's
# installation folder
if os.environ.has_key('OOMMF_PATH'):
    OOMMF_PATH = os.environ['OOMMF_PATH']
else:
    OOMMF_PATH = '/home/ww1g11/Softwares/oommf-1.2a5/'

if os.environ.has_key('OOMMF_TKTCL_VERSION'):
    OOMMF_TKTCL_VERSION = os.environ['OOMMF_TKTCL_VERSION']
else:
    OOMMF_TKTCL_VERSION = ""

# The path where the script is being executed
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

"""
Here we create a block of text with an OOMMF simulation

The mesh parameters are taken from the Fidimag mesh:
    length_x, length_y, length_z
    cellsize_x, cellsize_y, cellsize_z

Magnetic parameters will be specified manually. This simulation has Exchange,
DMI and Demag.  OOMMF will relax the system with the RK45 method and a specific
field is saved in a OMF file (see the Schedule specification at the end)

The saturation magnetisation Ms and the initial magnetisation m0 can be
modified to get spatial variation, providing a function that is going to be
used in a proc instance at the end of the script

"""
mif_demag = """# MIF 2.1

Specify Oxs_MultiAtlas:atlas {
    atlas { Oxs_BoxAtlas:world {
        xrange {0 %(length_x)s}
        yrange {0 %(length_y)s}
        zrange {0 %(length_z)s}
        name world
    } }
}

Specify Oxs_RectangularMesh:mesh [subst {
  cellsize {%(cellsize_x)s %(cellsize_y)s %(cellsize_z)s}
  atlas :atlas
}]

Specify Oxs_UniformExchange {
  A  %(A)s
}

Specify Oxs_DMExchange6Ngbr [subst {
    default_D %(D)s
    atlas :atlas
    D {
        world world %(D)s
    }
}]

Specify Oxs_Demag {}

Specify Oxs_RungeKuttaEvolve:evolve {
  alpha 0.0
  gamma_G 0.0
}

Specify Oxs_TimeDriver [subst {
 basename %(base_name)s
 evolver :evolve
 total_iteration_limit 1
 mesh :mesh

 Ms %(Ms)s


m0 { Oxs_ScriptVectorField {
    atlas :atlas
    script_args {rawpt}
    script init_m0
    norm 1
  } }
}]

proc init_Ms {x y z} {
    %(spatial_Ms)s
}

proc init_m0 { x y z} {
    %(init_m0)s
}

Destination archive mmArchive
Schedule DataTable archive Stage 1
Schedule Oxs_%(field)s::Field archive Stage 1
"""


def gen_oommf_conf(mesh, init_m0, A=1.3e-11, Ms=8e5, D=0,
                   spatial_Ms=None, field='Demag'):
    """

    Generate an OOMMF simulation script using the base text at the beginning of
    this code. This script will generate a MIF file with the name of the field
    specified in the field variable which must be a valid OOMMF field name


    mesh        :: A Fidimag mesh which is going to be reproduced
                   in OOMMF

    init_m0     :: The initial magnetisation profile, as a function
                   in TCL language, in a string. This string will be wrapped
                   around
                            proc init_m0 { x y z} {
                                --init_m0--
                            }
                   Then we use $x, $y, $z to vary m spatially and return
                   a list with the mx, my, mz components
                   (see OOMMF manual for details)

    A, D, Ms    :: Magnetic parameters (exchange, dmi, sat magnetisation)

    spatial_Ms  :: Instead of using an uniform Ms, it can be passed
                   a function (in TCL language) to vary Ms spatialy,
                   as a STRING. This string will be wrapped around:
                           proc init_Ms { x y z} {
                               --spatial_Ms--
                           }
                   Then we can use $x, $y, $z for the spatial
                   variables to modify Ms (see OOMMF manual for details),
                   for example:
                        spatial_Ms = \"\"\"
                        if { $x * $x + $y * $y < 5e-9 * 5e-9 } {
                            return 2e4
                        } else {
                            return 0
                        }
                        \"\"\"

    field       :: A string with the name of one of the OOMMF fields
                   used in the simulation. This name is used to save the field
                   to an OHf file and the script to a MIF file.

    """
    conf_path = os.path.join(MODULE_DIR, field)
    if not os.path.exists(conf_path):
        os.makedirs(conf_path)

    # Specify the mesh with a Fidimag's mesh object
    dx = mesh.dx * mesh.unit_length
    dy = mesh.dy * mesh.unit_length
    dz = mesh.dz * mesh.unit_length

    # If we have a space variational Ms, we add the corresponding
    # option inside the TimeDriver in the code
    # The init_Ms script, which was completed with the spatial_Ms
    # variable, is added at the end of the code
    if spatial_Ms is not None:
        Ms = """{ Oxs_ScriptScalarField {
            atlas :atlas
            script_args {rawpt}
            script  init_Ms
        } }"""
    # Otherwise, just use a constant Ms for a cubic mesh
    else:
        Ms = "%0.16g" % Ms

    # Replace the corresponding options from the base script
    # at the beginning of the code
    params = {
        'length_x': "%0.16g" % (dx * mesh.nx),
        'length_y': "%0.16g" % (dy * mesh.ny),
        'length_z': "%0.16g" % (dz * mesh.nz),
        'cellsize_x': "%0.16g" % dx,
        'cellsize_y': "%0.16g" % dy,
        'cellsize_z': "%0.16g" % dz,
        'base_name': field.lower(),
        'A': "%0.16g" % A,
        'D': "%0.16g" % D,
        'Ms': Ms,
        'spatial_Ms': spatial_Ms,
        'init_m0': init_m0,
        'field': field
    }

    mif = mif_demag % params

    # Write the MIF file for OOMMF with the field name
    with open(os.path.join(conf_path, field + ".mif"), "w") as mif_file:
        mif_file.write(mif)


def run_oommf(field='Demag'):
    """
    Run the OOMMF simulation specifying a valid OOMMF field
    to be saved into an OHF file
    """

    command = ('tclsh{}'.format(OOMMF_TKTCL_VERSION),
               os.path.join(OOMMF_PATH, 'oommf.tcl'),
               'boxsi',
               '-threads',
               '1',
               field + ".mif")

    cmd = ' '.join(command)
    logger.info("About to execute '{}'".format(cmd))
    print("About to execute '{}'".format(cmd))

    save_path = os.getcwd()
    new_path = os.path.join(MODULE_DIR, field)

    os.chdir(new_path)

    os.system(cmd)

    os.chdir(save_path)


def extract_data(mesh, ovf_file):
    """
    Extract the magnetisation components from an OVF file
    """

    ovf = omf.OMF2(ovf_file)

    mx = np.zeros(mesh.nxyz)
    my = np.zeros(mesh.nxyz)
    mz = np.zeros(mesh.nxyz)

    for i in range(mesh.nx):
        for j in range(mesh.ny):
            for k in range(mesh.nz):
                id_n = mesh.index(i, j, k)
                mx[id_n] = ovf.get_mag(i, j, k, comp='x')
                my[id_n] = ovf.get_mag(i, j, k, comp='y')
                mz[id_n] = ovf.get_mag(i, j, k, comp='z')

    m = np.array([mx, my, mz])
    m.shape = (-1,)
    return m


def get_field(mesh,  field='Demag'):
    """
    Return the field components of every node, extracting the
    data from an OVF file which is named with a valid
    OOMMF field (Demag by default)
    """
    new_path = os.path.join(MODULE_DIR, field)
    file_name = '%s-Oxs_%s-Field-00-0000001.ohf' % (field.lower(), field)
    ovf_file = os.path.join(new_path, file_name)
    ovf = omf.OMF2(ovf_file)

    return ovf.get_all_mags()


# The following functions automatically extract the field
# components for 3 different OOMMF fields: Exchange, Demag and DMI

def compute_demag_field(mesh, init_m0, Ms=8e5, spatial_Ms=None, field='Demag'):

    gen_oommf_conf(
        mesh, Ms=Ms, init_m0=init_m0, spatial_Ms=spatial_Ms, field=field)

    run_oommf(field)

    m = get_field(mesh, field=field)

    new_path = os.path.join(MODULE_DIR, field)

    command = ('rm',
               '-rf',
               new_path)
    cmd = ' '.join(command)
    os.system(cmd)

    return m


def compute_exch_field(mesh, init_m0, Ms=8e5, A=1.3e-11,
                       spatial_Ms=None, field='UniformExchange'):

    gen_oommf_conf(
        mesh, Ms=Ms, init_m0=init_m0, A=A, spatial_Ms=spatial_Ms, field=field)
    run_oommf(field)

    m = get_field(mesh, field=field)

    new_path = os.path.join(MODULE_DIR, field)

    command = ('rm',
               '-rf',
               new_path)
    cmd = ' '.join(command)
    os.system(cmd)

    return m


def compute_dmi_field(mesh, init_m0, Ms=8e5, D=4e-3, field='DMExchange6Ngbr'):

    gen_oommf_conf(mesh, Ms=Ms, init_m0=init_m0, D=D, field=field)
    run_oommf(field)

    m = get_field(mesh, field=field)

    new_path = os.path.join(MODULE_DIR, field)

    command = ('rm',
               '-rf',
               new_path)
    cmd = ' '.join(command)
    os.system(cmd)

    return m

if __name__ == "__main__":

    mesh = CuboidMesh(nx=5, ny=2, nz=1, dx=1.0, dy=1.0, dz=1.0)
    m = compute_demag_field(mesh, init_m0='return "1 0 0"')
    print m
