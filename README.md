# Fidimag


[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/computationalmodelling/fidimag/master)
[![Documentation Status](https://readthedocs.org/projects/fidimag/badge/?version=latest)](http://fidimag.readthedocs.org/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/computationalmodelling/fidimag/branch/master/graph/badge.svg)](https://codecov.io/gh/computationalmodelling/fidimag)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.167858.svg)](https://doi.org/10.5281/zenodo.167858)
[![Website](https://img.shields.io/website-up-down-green-red/http/shields.io.svg?label=Fidimag-Website)](http://computationalmodelling.github.io/fidimag/)

| Tests Type | Status |
|:-:|:-:|
| GHA Unit Tests | ![Build](https://github.com/computationalmodelling/fidimag/actions/workflows/build.yml/badge.svg?)
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

[32] [Thermal Evolution of Skyrmion Formation Mechanism in Chiral Multilayer Films](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.17.044039) Phys. Rev. Applied 17, 044039 (2022)

[31] [Mutual conversion between a magnetic Néel hopfion and a Néel toron](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.174407) Phys. Rev. B 105, 174407 (2022)

[30] [Unveiling the emergent traits of chiral spin textures in magnetic multilayers](https://onlinelibrary.wiley.com/doi/abs/10.1002/advs.202103978) Advanced Science Vol 9, Iss 6 (2022)

[29] [The magnetic genome of two-dimensional van der waals materials](https://pubs.acs.org/doi/full/10.1021/acsnano.1c09150) ACS Nano 16, 5, 6960–7079 (2022)

[28] [L-shaped electrode design for high-density spin–orbit torque magnetic random access memory with perpendicular shape anisotropy](https://iopscience.iop.org/article/10.1088/1361-6463/abf61d/) J. Phys. D: Appl. Phys. 54 285002 (2021)

[27] [Periodically modulated skyrmion strings in Cu2OSeO3](https://doi.org/10.1038/s41535-021-00373-y), npj Quantum Materials volume 6, Article number: 73 (2021)

[26] [Speeding up explicit numerical evaluation methods for micromagnetic simulations using demagnetizing field polynomial extrapolation](https://ieeexplore.ieee.org/document/9737008), IEEE Transactions on Magnetics, Vol 58 Issue 5 (2022)

[25] [A Review of Modelling in Ferrimagnetic Spintronics](https://journals.jps.jp/doi/full/10.7566/JPSJ.90.081001) J. Phys. Soc. Jpn. 90, 081001 (2021)

[24] [Topological defect-mediated skyrmion annihilation in three dimensions](https://www.nature.com/articles/s42005-021-00675-4) Communications Physics 4, 175 (2021)

[23] [Stray Field Calculation for Micromagnetic Simulations Using True Periodic Boundary Conditions](https://doi.org/10.1038/s41598-021-88541-9) Scientific Reports  11, 9202 (2021)

[22] [Field-free spin–orbit torque perpendicular magnetization switching in ultrathin nanostructures](https://doi.org/10.1038/s41524-020-0347-0), npj Computational Materials volume 6, 78 (2020)

[21] [Hybrid FFT algorithm for fast demagnetization field calculations on non-equidistant magnetic layers](https://doi.org/10.1016/j.jmmm.2020.166592)
Journal of Magnetism and Magnetic Materials, Volume 503, 166592 (2020)

[20] [Review – Micromagnetic Simulation Using OOMMF and Experimental Investigations on Nano Composite Magnets](https://doi.org/10.1088/1742-6596/1172/1/012070)
Review – Micromagnetic Simulation Using OOMMF and Experimental Investigations on Nano Composite Magnets, J. Phys.: Conf. Ser. 1172 012070

[19] [Spin waves in thin films and magnonic crystals with Dzyaloshinskii-Moriya interactions](https://arxiv.org/abs/1903.04288), arxiv:1903.04288 (2019)

[18] [Tomorrow's Micromagnetic Simulations](https://doi.org/10.1063/1.5093730), Journal of Applied Physics 125, 180901 (2019)

[17] [Diameter-independent skyrmion Hall angle in the plastic flow regime observed in chiral magnetic
multilayers](https://arxiv.org/pdf/1908.04239.pdf](https://www.nature.com/articles/s41467-019-14232-9), Nature Communications volume 11, Article number: 428 (2020)

[16] [Efficient computation of demagnetising fields for magnetic multilayers using multilayered convolution](https://aip.scitation.org/doi/10.1063/1.5116754), Journal of Applied Physics 126, 103903 (2019)

[15] [Micromagnetics and spintronics: models and numerical methods](https://link.springer.com/article/10.1140%2Fepjb%2Fe2019-90599-6), Eur. Phys. J. B (2019) 92: 120

[14] [Nanoscale magnetic skyrmions and target states in confined geometries](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.214408), Physical Review B 99, 214408 (2019) 

[13] [Learning Magnetization Dynamics](https://www.sciencedirect.com/science/article/abs/pii/S0304885319307978?via%3Dihub), Journal of Magnetism and Magnetic Materials
Volume 491, (2019)

[12] [Computational micromagnetics with Commics](https://doi.org/10.1016/j.cpc.2019.106965), Computer Physics Communications, 248 (2020)

[11] [Binding a hopfion in a chiral magnet nanodisk](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.98.174437), Physical Review B 98, 174437 (2018)

[10] [Proposal for a micromagnetic standard problem for materials with Dzyaloshinskii–Moriya interaction](http://iopscience.iop.org/article/10.1088/1367-2630/aaea1c), New Journal of Physics, Volume 20 (2018)

[9] [Driving chiral domain walls in antiferromagnets using rotating magnetic fields](https://link.aps.org/doi/10.1103/PhysRevB.97.184418) Physical Review B 97, 184418 (2018)

[8] [Fidimag - A Finite Difference Atomistic and Micromagnetic Simulation Package](http://doi.org/10.5334/jors.223), Journal of Open Research Software, 6(1), p.22. (2018)

[7] [Topological Spintronics in Confined Geometry](https://escholarship.org/uc/item/8wx626mw), Y. Liu, PhD Thesis, University of California Riverside (2017)

[6] [Thermal stability and topological protection of skyrmions in nanotracks](https://www.nature.com/articles/s41598-017-03391-8), Scientific Reports 7, 4060 (2017)

[5] [Current-induced instability of domain walls in cylindrical nanowires](http://iopscience.iop.org/article/10.1088/1361-648X/aa9698/meta), Journal of Physics: Condensed Matter, 30, 1 (2017)

[4] [Magnonic analog of relativistic Zitterbewegung in an antiferromagnetic spin chain](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.96.024430), Phys. Rev. B 96 024430 (2017)

[3] [Driving magnetic skyrmions with microwave fields](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.020403) Phys. Rev. B 92, 020403 (2015).

[2] [Microwave-induced dynamic switching of magnetic skyrmion cores in nanodots](https://aip.scitation.org/doi/10.1063/1.4914496) Applied Physics Letters 106, 102401 (2015).

[1] [Magnon-Driven Domain-Wall Motion with the Dzyaloshinskii-Moriya Interaction](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.087203) Phys. Rev. Lett. 114, 087203 (2015)

### Acknowledgements

We acknowledge financial support from EPSRC’s Centre for Doctoral Training in Next Generation Computational Modelling (EP/L015382/1),  EPSRC’s Doctoral Training Centre in Complex System Simulation (EP/G03690X/1), EPSRC Programme grant on Skyrmionics (EP/N032128/1) and OpenDreamKitHorizon 2020 European Research Infrastructure project (676541).
