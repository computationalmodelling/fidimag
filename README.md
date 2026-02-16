# Fidimag


[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/computationalmodelling/fidimag/master)
[![Documentation Status](https://readthedocs.org/projects/fidimag/badge/?version=latest)](http://fidimag.readthedocs.org/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/computationalmodelling/fidimag/branch/master/graph/badge.svg)](https://codecov.io/gh/computationalmodelling/fidimag)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.167858.svg)](https://doi.org/10.5281/zenodo.167858)
[![Website](https://img.shields.io/website-up-down-green-red/http/shields.io.svg?label=Fidimag-Website)](http://computationalmodelling.github.io/fidimag/)

| Tests Type | Status |
|:-:|:-:|
| GHA Unit Tests | ![Build](https://github.com/computationalmodelling/fidimag/actions/workflows/build.yml/badge.svg?)

<img src="http://computationalmodelling.github.io/fidimag/figs/skyrmion.jpg" alt="Fidimag Image" width="400" align="right">

Fidimag solves finite-difference micromagnetic problems and supports atomistic simulations, using Python interface. The interface to both types of simulation is similar.

### Install

**Quick start** (using modern build system with uv):

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone and build
git clone https://github.com/computationalmodelling/fidimag.git
cd fidimag
uv sync

# Activate environment
source .venv/bin/activate
```

For detailed installation instructions and troubleshooting, see:
- **Build documentation**: [BUILD.md](BUILD.md)
- **Migration guide**: [MIGRATION.md](MIGRATION.md)
- **Online docs**: [https://fidimag.readthedocs.io/en/latest/install.html](https://fidimag.readthedocs.io/en/latest/install.html)

### Features
* Optimal LLG equation integration using modern [Sundial's v7](https://github.com/LLNL/sundials/) CVODE solver
* Offers LLG and LLG with spin torque terms (Zhang-Li and Sloncewski)
* Calculations using the Geodesic-Nudged-Elastic-Band and String methods to compute energy barriers.
* Exchange, Zeeman, Demagnetising, Uniaxial Anisotropy energy classes.
* Parallelised using OpenMP.
* Easily extensible to add new features.
* Cubic and Hexagonal Meshes in atomistic simulations.
* Open-source under the 2-clause BSD Licence.

### Example
Here we show how to relax a nanodisk system from an initial state. We have many more examples in the [documentation](http://fidimag.readthedocs.io/en/latest/?badge=latest)!

```python
import fidimag
from fidimag.common import CuboidMesh
from fidimag.micro import Sim, UniformExchange, Demag, DMI, UniaxialAnisotropy
mesh = CuboidMesh(nx=60, ny=60, nz=1, dx=2.0, dy=2.0, dz=2.0, unit_length=1e-9)

def Ms_init(position):
    """
    Set where the system has magnetic material
    Form a nanodisk shape
    """
    Ms = 8.6e5
    x, y, z = position
    if (x - 60)**2 + (y - 60)**2 < 60**2:
        return Ms
    else:
        return 0

def m_init(position):
    """
    Approximate skyrmion profile
    """
    x, y, z = position
    if (x - 60)**2 + (y - 60)**2 < 40**2:
        return (0, 0, 1)
    else:
        return (0, 0, -1)

sim = Sim(mesh, name='target_skyrmion')
sim.set_Ms(Ms_init)
sim.set_m(m_init)
sim.add(Demag())
sim.add(UniformExchange(A=1e-11))
sim.add(DMI(D=3e-3))
sim.add(UniaxialAnisotropy(Ku=4e5, axis=(0, 0, 1)))
sim.relax()
sim.save_vtk()
```
The results can be straightforwardly visualised from the outputted VTK files using programs such as Paraview:
<p align="center">
<img src="http://computationalmodelling.github.io/fidimag/figs/target.png" alt="Target Skyrmion State" width="250">
</p>

### Publications

A list tracking publications that cite FIDIMAG can be found in the [Publications](PUBLICATIONS.md) document.


### Attributions
The code is currently developed by Weiwei Wang, David Cortes, Ryan Pepper and Hans Fangohr at Anhui University, University of Southampton and European XFEL. Previous developers who have worked on Fidimag are Marc-Antonio Bisotti, Thomas Kluyver, Mark Vousden, Oliver Laslett and Rebecca Carey.

Contributions and pull requests to both the code and documentation are welcome.
If you use Fidimag, please cite as:

Bisotti, M.-A., Cortés-Ortuño, D., Pepper, R., Wang, W., Beg, M., Kluyver, T. and Fangohr, H., 2018. Fidimag – A Finite Difference Atomistic and Micromagnetic Simulation Package. Journal of Open Research Software, 6(1), p.22. DOI: http://doi.org/10.5334/jors.223


### Acknowledgements

We acknowledge financial support from EPSRC’s Centre for Doctoral Training in Next Generation Computational Modelling (EP/L015382/1),  EPSRC’s Doctoral Training Centre in Complex System Simulation (EP/G03690X/1), EPSRC Programme grant on Skyrmionics (EP/N032128/1) and OpenDreamKitHorizon 2020 European Research Infrastructure project (676541).
