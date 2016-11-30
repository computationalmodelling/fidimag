Nudged Elastic Band Method (NEBM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NEBM is an algorithm to find minimum energy transitions between equilibrium
states. This method was developed by Henkelman and Jónsson [1] to solve
molecular problems in Chemistry. Later, Dittrich et al. [2] applied the NEBM to
magnetic systems based on the micromagnetic theory. Recently, Bessarab et al.
[3] have proposed an optimised version of the NEBM with improved control over
the behaviour of the algorithm. Their review of the method allows to apply the
technique to systems described by either a micromagnetic or atomistic model.

The algorithm is based on, firstly, generating :math:`N` copies of a magnetic system.
Every copy of the system, which is described by :math:`P` spins or magnetic
moments (arranged in a lattice or mesh), is called an *image*
:math:`\mathbf{Y}_{i}`. This sequence of images defines a *band*, where every
image is in a (ideally) different magnetic configuration. The energy of the
system depends on the magnetic configuration (e.g. a skyrmion, vortex, uniform
state,etc.), thus the energy is parametrised by the number of degrees of
freedom :math:`P`. Accordingly, every image will have a specific energy, i.e.
:math:`E=E(\mathbf{Y})`,  which determines the position of an image in an
energy landscape. The first and last images in a band are chosen as equilibrium
states of the magnetic system and they are kept fixed during the NEBM
evolution.

Secondly, to initiate the NEBM evolution, it is necessary to specify an initial
guess for the images, which means generating different magnetic configurations
between the extrema of the band. It is customary to use an interpolation of the
spins directions to achieve this.

And thirdly, the band is relaxed to find a path in energy space that costs less
energy. This path is characterised by passing through a saddle point, which is
a maximum in energy along certain directions in phase space, and this point
determines an energy barrier between the two equilibrium states. The energy
barrier is the energy necessary to drive one equilibrium state towards the
other. A first order saddle point is the one that is a maximum along a single
direction in phase space and will usually in the energy path that costs less
energy. It is possible that there are more than one energy paths.

To distinguish different images in the energy landscape it is necessary to
define a *distance*, and to keep the images equally spaced to avoid clustering
around saddle points or equilibrium states it is defined a spring force between
images.

For a thorough explanation of the method see references [1,4].
