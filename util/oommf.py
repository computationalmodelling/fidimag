import os
import sys
import os
import numpy as np
import subprocess
import omf

from pc.fd_mesh import FDMesh

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

mif_demag = """# MIF 2.1

Specify Oxs_BoxAtlas:atlas {
       xrange {0 %(length_x)s}
       yrange {0 %(length_y)s}
       zrange {0 %(length_z)s}
}

Specify Oxs_RectangularMesh:mesh [subst {
  cellsize {%(cellsize_x)s %(cellsize_y)s %(cellsize_z)s}
  atlas :atlas
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

proc init_m0 { x y z} {

    %(init_m0)s

}

Destination archive mmArchive

Schedule Oxs_Demag::Field archive Stage 1
"""

def gen_conf_demag(mesh, init_m0, Ms=8e5, base_name='demag'):

    conf_path = os.path.join(MODULE_DIR, base_name)
    if not os.path.exists(conf_path):
        os.makedirs(conf_path)
    
    dx = mesh.dx*mesh.unit_length
    dy = mesh.dy*mesh.unit_length
    dz = mesh.dz*mesh.unit_length
    
    params = {
        'length_x': "%0.16e" %(dx*mesh.nx),
        'length_y': "%0.16e" %(dy*mesh.ny),
        'length_z': "%0.16e" %(dz*mesh.nz),
        'cellsize_x': "%0.16e" % dx,
        'cellsize_y': "%0.16e" % dy,
        'cellsize_z': "%0.16e" % dz,
        'base_name': base_name,
        'Ms': "%0.16e" % Ms,
        'init_m0': init_m0,
    }
    
    mif = mif_demag%params

    with open(os.path.join(conf_path, base_name+".mif"), "w") as mif_file:
        mif_file.write(mif)

def run_oommf(base_name='demag'):

    command = ('tclsh',
           '/home/ww1g11/Softwares/oommf-1.2a5/oommf.tcl',
           'boxsi',
           '-threads',
           '1',
           base_name+".mif")
    
    cmd = ' '.join(command)
   
    save_path=os.getcwd()
    new_path=os.path.join(MODULE_DIR, base_name)
 
    os.chdir(new_path)
    
    os.system(cmd)

    os.chdir(save_path)


def get_field(mesh, base_name='demag', Field='Demag'):
    new_path = os.path.join(MODULE_DIR, base_name)
    file_name = '%s-Oxs_%s-Field-00-0000001.ohf'%(base_name,Field)
    ovf_file = os.path.join(new_path, file_name)
    ovf = omf.OMF2(ovf_file)
    
    mx = np.zeros(mesh.nxyz)
    my = np.zeros(mesh.nxyz)
    mz = np.zeros(mesh.nxyz)
    
    for i in range(mesh.nx):
        for j in range(mesh.ny):
            for k in range(mesh.nz):
                id_n = mesh.index(i,j,k)
                mx[id_n] = ovf.get_mag(i,j,k,comp='x')
                my[id_n] = ovf.get_mag(i,j,k,comp='y')
                mz[id_n] = ovf.get_mag(i,j,k,comp='z')

    m = np.array([mx,my,mz])
    m.shape = (-1,)
    return m


def compute_demag_field(mesh, init_m0, Ms=8e5,base_name='demag'):
    gen_conf_demag(mesh, Ms=Ms, init_m0=init_m0, base_name=base_name)
    run_oommf(base_name)
    
    m = get_field(mesh,base_name=base_name)
    
    new_path=os.path.join(MODULE_DIR, base_name)
    
    command = ('rm',
           '-rf',
           new_path)
    cmd = ' '.join(command)
    os.system(cmd)
    
    return m
    

if __name__=="__main__":
    
    mesh=FDMesh(nx=5,ny=2,nz=1,dx=1.0,dy=1.0,dz=1.0)
    gen_conf_demag(mesh, 'return "1 0 0"')
    run_oommf()
    m = get_field(mesh)
    print m
    
    
