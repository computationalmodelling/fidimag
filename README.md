# Fidimag


[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/computationalmodelling/fidimag/master)
[![Documentation Status](https://readthedocs.org/projects/fidimag/badge/?version=latest)](http://fidimag.readthedocs.org/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/computationalmodelling/fidimag/branch/master/graph/badge.svg)](https://codecov.io/gh/computationalmodelling/fidimag)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.167858.svg)](https://doi.org/10.5281/zenodo.167858)
[![Website](https://img.shields.io/website-up-down-green-red/http/shields.io.svg?label=Fidimag-Website)](http://computationalmodelling.github.io/fidimag/)

| Tests Type | Status |
|:-:|:-:|
| Tests | ![Tests](https://travis-matrix-badges.herokuapp.com/repos/computationalmodelling/fidimag/branches/master/1)
| Notebooks | ![Notebooks](https://travis-matrix-badges.herokuapp.com/repos/computationalmodelling/fidimag/branches/master/2)

<img src="http://computationalmodelling.github.io/fidimag/figs/skyrmion.jpg" alt="Fidimag Image" width="400" align="right">

Fidimag solves finite-difference micromagnetic problems and supports atomistic simulations, using Python interface. The interface to both types of simulation is similar.

### Features
* Offers LLG and LLG with spin torque terms (Zhang-Li and Sloncewski)
* Calculations using the Nudged-Elastic-Band method to compute energy barriers.
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




### Attributions
The code is currently developed by Weiwei Wang, David Cortes, Ryan Pepper and Hans Fangohr at Anhui University, University of Southampton and European XFEL. Previous developers who have worked on Fidimag are Marc-Antonio Bisotti, Thomas Kluyver, Mark Vousden, Oliver Laslett and Rebecca Carey.

Contributions and pull requests to both the code and documentation are welcome.
If you use Fidimag, please cite as:

Bisotti, M.-A., Cortés-Ortuño, D., Pepper, R., Wang, W., Beg, M., Kluyver, T. and Fangohr, H., 2018. Fidimag – A Finite Difference Atomistic and Micromagnetic Simulation Package. Journal of Open Research Software, 6(1), p.22. DOI: http://doi.org/10.5334/jors.223

### Publications

The following publications, in reverse chronological order, have used or cited Fidimag:

[14] [Learning Magnetization Dynamics](https://arxiv.org/abs/1903.09499), A. Kovacs, J. Fischbacher, H. Oezelt, M. Gusenbauer, L. Exl, F. Bruckner, D.Suess, T.Schrefl, arXiV (2019)

[13] [Computational micromagnetics with Commics](https://arxiv.org/abs/1812.05931), C.-M. Pfeiler, M. Ruggeri, B. Stiftner, L. Exl, M. Hochsteger, G. Hrkac, J. Schöberl, N. J. Mauser, D. Praetorius, arXiV (2018)

[12] [Micromagnetics and spintronics: Models and numerical methods](https://arxiv.org/abs/1810.12365), C. Abert, arXiV (2018)

[11] [Binding a hopfion in a chiral magnet nanodisk](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.98.174437), Y. Liu, R. K. Lake, and J. Zang, Physical Review B 98, 174437 (2018)

[10] [Proposal for a micromagnetic standard problem for materials with Dzyaloshinskii–Moriya interaction](http://iopscience.iop.org/article/10.1088/1367-2630/aaea1c), D. Cortés-Ortuño, M. Beg, V. Nehruji, L. Breth, R. Pepper, T. Kluyver, G. Downing, T. Hesjedal, P. Hatton, T. Lancaster, R. Hertel, O. Hovorka and H. Fangohr, New Journal of Physics, Volume 20 (2018)

[9] [Driving chiral domain walls in antiferromagnets using rotating magnetic fields](https://link.aps.org/doi/10.1103/PhysRevB.97.184418) K.Pan, L.Xing, H.Y.Yuan, and W.Wang, Physical Review B 97, 184418 (2018)

[8] [Fidimag - A Finite Difference Atomistic and Micromagnetic Simulation Package](http://doi.org/10.5334/jors.223), Bisotti, M.-A., Cortés-Ortuño, D., Pepper, R., Wang, W., Beg, M., Kluyver, T. and Fangohr, H., Journal of Open Research Software, 6(1), p.22. (2018)

[7] [Topological Spintronics in Confined Geometry](https://escholarship.org/uc/item/8wx626mw), Y. Liu, PhD Thesis, University of California Riverside (2017)

[6] [Thermal stability and topological protection of skyrmions in nanotracks](https://www.nature.com/articles/s41598-017-03391-8), D. Cortés-Ortuño, W. Wang, M. Beg, R.A. Pepper, M-A. Bisotti, R. Carey, M. Vousden, T. Kluyver, O. Hovorka & H. Fangohr, Scientific Reports 7, 4060 (2017)

[5] [Current-induced instability of domain walls in cylindrical nanowires](http://iopscience.iop.org/article/10.1088/1361-648X/aa9698/meta), W. Wang, Z. Zhang, R.A. Pepper, C. Mu, Y. Zhou & Hans Fangohr, Journal of Physics: Condensed Matter, 30, 1 (2017)

[4] [Magnonic analog of relativistic Zitterbewegung in an antiferromagnetic spin chain](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024430), W. Wang, C. Gu, Y. Zhou & H. Fangohr, Phys. Rev. B 96 024430 (2017)

[3] [Driving magnetic skyrmions with microwave fields](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.020403) W. Wang, M. Beg, B. Zhang, W. Kuch, and H. Fangohr, Phys. Rev. B 92, 020403 (2015).

[2] [Microwave-induced dynamic switching of magnetic skyrmion cores in nanodots](https://aip.scitation.org/doi/10.1063/1.4914496) B. Zhang, W. Wang, M. Beg, H. Fangohr, and W. Kuch, Applied Physics Letters 106, 102401 (2015).

[1] [Magnon-Driven Domain-Wall Motion with the Dzyaloshinskii-Moriya Interaction](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.087203) W. Wang, M. Albert, M. Beg, M-A. Bisotti, D. Chernyshenko, D. Cortés-Ortuño, I. Hawke & H. Fangohr, Phys. Rev. Lett. 114, 087203 (2015)

### Acknowledgements

We acknowledge financial support from EPSRC’s Centre for Doctoral Training in Next Generation Computational Modelling (EP/L015382/1),  EPSRC’s Doctoral Training Centre in Complex System Simulation (EP/G03690X/1), EPSRC Programme grant on Skyrmionics (EP/N032128/1) and OpenDreamKitHorizon 2020 European Research Infrastructure project (676541).
