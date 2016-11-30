Nudged Elastic Band Method (NEBM)
=================================

The NEBM is an algorithm to find minimum energy transitions between equilibrium
states. This method was developed by Henkelman and Jónsson [1] to solve
molecular problems in Chemistry. Later, Dittrich et al. [2] applied the NEBM to
magnetic systems based on the micromagnetic theory. Recently, Bessarab et al.
[3] have proposed an optimised version of the NEBM with improved control over
the behaviour of the algorithm. Their review of the method allows to apply the
technique to systems described by either a micromagnetic or atomistic model.

The algorithm is based on, firstly, generating :math:`N` copies of a magnetic
system.  Every copy of the system, which is described by :math:`P` spins or
magnetic moments (arranged in a lattice or mesh), is called an *image*
:math:`\mathbf{Y}_{i}`. This sequence of images defines a *band*, where every
image is in a (ideally) different magnetic configuration. The energy of the
system depends on the magnetic configuration (e.g. a skyrmion, vortex, uniform
state,etc.), thus the energy is parametrised by the number of degrees of
freedom, which is :math:`n\times P`, where :math:`n` is the number of
coordinates to describe a spin (e.g. :math:`n=3` for Cartesian coordinates and
:math:`n=2` in Spherical coordinates). Accordingly, every image will have a
specific energy, i.e.  :math:`E=E(\mathbf{Y})`,  which determines the position
of an image in an energy landscape. The first and last images in a band are
chosen as equilibrium states of the magnetic system and they are kept fixed
during the NEBM evolution.

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

It is also necessary to consider that:

* To distinguish different images in the energy landscape we need to define a
  *distance*

* To keep the images equally spaced to avoid clustering around saddle points or
  equilibrium states, we use a spring force between images.

* Spins or magnetic moments can be described in any coordinate system. The most
  commonly used are spherical and Cartesian coordinates. For the later we have
  to specify the constraint to fix the spin length.


For a thorough explanation of the method see references [3,4].


NEBM relaxation
---------------

Cartesian
^^^^^^^^^

When using Cartesian coordinates, every image of the band
:math:`\mathbf{Y}_{i}` is iterated with the following dynamical equation with a
fictitious time :math:`tau`

.. math::
    \frac{\partial \mathbf{Y}_i}{\partial \tau} = -\gamma \mathbf{Y}_{i} \times
    \mathbf{Y}_{i} \times \mathbf{G}_{i} + c \sqrt{ \left( \frac{\partial \mathbf{Y}_{i}}{\partial \tau} \right)^{2} }
			\left( 1 - \mathbf{Y}_{i}^{2} \right) \mathbf{Y}_{i}

In this equation, :math:`\gamma` is in units of :math:`\text{Hz T}^{-1}` and it
determines the time scale, which is irrelevant here so we set :math:`\gamma=1`.
The second term to the right is necessary to keep the spins/magnetisation length
equal to one, using an appropriate factor :math:`c` that we set to 6. The :math:`\mathbf{G}`
is the total force on an image, which in the atomistic theory is defined as

.. math::
    \mathbf{G}_{i} =  - \boldsymbol{\nabla}_{\boldsymbol{\mu}} E(\mathbf{Y}_{i})|_{\perp} +
                 \mathbf{F}(\mathbf{Y}_{i})|_{\parallel}

and in the micromagnetic theory as

.. math::
    \mathbf{G}_{i} =  - \boldsymbol{\nabla}_{\mathbf{M}} E(\mathbf{Y}_{i})|_{\perp} +
                 \mathbf{F}(\mathbf{Y}_{i})|_{\parallel}

We use the magnetic effective field definition to evaluate the gradients, i.e. 
:math:`\boldsymbol{\nabla}_{\boldsymbol{\mu}}E=\partial E / (\mu_{s}\partial\mathbf{s})=-\mathbf{H}_{\text{eff}}`
or :math:`\boldsymbol{\nabla}_{\mathbf{M}}E=\partial E / (M_{s}\partial\mathbf{m})=-\mathbf{H}_{\text{eff}}`.
The perpendicular component is with respect to the tangents :math:`\mathbf{t}` to the energy band, thus
for a vector :math:`\mathbf{A}`

.. math::
    \mathbf{A}|_{\perp} = \mathbf{A} - (\mathbf{A}\cdot\mathbf{t})\mathbf{t}

The tangents are defined according to the energies of the neighbouring images [3]. The second term
to the right hand side of the equation for :math:`\mathbf{G}` is the spring force that
tries to keep images at equal distance and is defined using the distance between neighbouring
images

.. math::
    \mathbf{F}(\mathbf{Y}_{i})|_{\parallel}=k\left(|\mathbf{Y}_{i+1}-\mathbf{Y}_{i}|-
        |\mathbf{Y}_{i}-\mathbf{Y}_{i-1}|\right)\mathbf{t}_{i}

which is parallel to the band, i.e. in the direction of the tangent.

Following the definitions from Bessarab et al. [3], the tangents and the total
force :math:`\mathbf{G}` must be *projected* into the spin/magnetisation
tangent space.

Vectors
-------

We usually refer to the images, forces, etc. as *vectors*. In Fidimag, for a
magnetic system of :math:`P` spins/magnetic moments, we define some vectors as
:math:`n\times P` arrays, where :math:`n` is the number of coordinates to
describe a spin (e.g. :math:`n=3` for Cartesian coordinates and :math:`n=2` in
Spherical coordinates). For example, the :math:`i` image of a band, in
Cartesian coordinates, is defined as

.. math::
    \mathbf{Y}_{i} = \left( s_{x,0}^{(i)}, s_{y,0}^{(i)}, s_{z,0}^{(i)}, s_{x,1}^{(i)}, 
                     s_{y,1}^{(i)},\ldots s_{y,P-1}^{(i)}, s_{z,P-1}^{(i)}   
                     \right)

Within the micromagnetic model we use :math:`\mathbf{m}` rather than :math:`\mathbf{s}`.
Notice that we have the three components of every spin in the system.

For the total force (or tangents, spring forces, etc.) is similar, every spin
will have associated a total force:

.. math::
    \mathbf{G}_{i} = \left( G_{x,0}^{(i)}, G_{y,0}^{(i)}, G_{z,0}^{(i)}, G_{x,1}^{(i)}, 
                     G_{y,1}^{(i)},\ldots G_{y,P-1}^{(i)}, G_{z,P-1}^{(i)}   
                     \right)

Projections
-----------

The projection of a vector into the spin/magnetisation tangent space simply
means projecting its components with the corresponding spin/magnetisation field
components. For example, for a vector :math:`\mathbf{A}`

.. math::
    \mathbf{A} = \left( \mathbf{A}_{0}, \ldots \mathbf{A}_{P-1}\right) = 
                 \left( A_{x,0}, A_{y,0}, A_{z,0}, A_{x,1}, A_{y,1},\ldots A_{y,P-1}, A_{z,P-1} \right)

the projection :math:`\mathcal{P}` is defined as

.. math::
    \mathcal{P}\mathbf{A} = \left( \mathcal{P}_{\mathbf{s}_{0}}\mathbf{A}_{0}, 
                                   \mathcal{P}_{\mathbf{s}_{1}}\mathbf{A}_{1},
                                   \ldots
                                   \mathcal{P}_{\mathbf{s}_{1}}\mathbf{A}_{P-1},
                            \right)

where

.. math::
   \mathcal{P}_{\mathbf{s}_{i}}\mathbf{A}_{i} =  \mathbf{A}_{i} - 
                        \left( \mathbf{A}_{i} \cdot \mathbf{s}_{i} \right) \mathbf{s}_{i}


Algorithm
---------

The algorithm can be summarised as:

1. Define a magnetic system and find two equilibrium states for which we want
   to find a minimum energy transition.

2. Set up a band of images and an initial sequence between the extrema. We can
   use linear interpolations on the spherical angles that define the spin
   directions [4] or Vicenty's formulae [3].

3. Evolve the system using the NEBM dynamical equation, which depends on the
   chosen coordinate system. This equation involves:
   
   I. Compute the effective field for every image (they are in different magnetic
      configurations) and the total energy of every image

   II. Compute the tangents according to the energies of the images and project them
       into the spin/magnetisation tangent space

   III. Compute the total force for every image in the band using the tangents
        and distances between neighbouring images. This allows to calculate the
        gradient (which uses the effective field) and the spring forces on the
        images

   IV. Project the total force into the spin/magnetisation tangent space

   V. Use the dynamical equation according to the coordinate system
