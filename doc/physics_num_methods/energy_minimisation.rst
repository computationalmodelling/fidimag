Energy minimisation
===================

To take a magnetic system towards an equilibrium state we do not necessarily
have to integrate the LLG equation in time. If we are only interested in the
final configuration, and not in the trajectory followed to reach it, we can
descend the energy directly, which is usually substantially cheaper. The
classes described here do this, and they are chosen through the ``driver``
argument of the simulation classes rather than through the integrators.

The problem is a constrained one. The energy :math:`E=E(\mathbf{m})` is a
function of the spin directions, and these have a fixed length, so what we
minimise is

.. math::
    E(\mathbf{m}) - \lambda \left( \mathbf{m}^{2} - 1 \right)

with :math:`\lambda` a Lagrange multiplier. Setting the derivative with respect
to :math:`\mathbf{m}` to zero leads to :math:`\mathbf{m}\times\mathbf{H}_{\text{eff}}=0`,
i.e. what vanishes at a minimum is the torque and not the energy gradient
itself. Notice that in spherical coordinates the constraint is implicit and the
conditions are simply :math:`\partial E/\partial\theta=\partial E/\partial\phi=0`,
but the classes in Fidimag work in Cartesian coordinates.


Gradient descent on the sphere
------------------------------

The gradient of the energy is obtained from the effective field, ignoring
scaling parameters such as :math:`\mu_0`, :math:`M_s` or :math:`\mu_s`,

.. math::
    \frac{\delta E}{\delta \mathbf{m}} = - \mathbf{H}_{\text{eff}}

This gradient has a component along :math:`\mathbf{m}`, which is the Lagrange
multiplier of the length constraint. That component is removed when we
re-normalise the spins after every update, so it tells us nothing about the
direction in which we should descend, and it is convenient to project it out
from the start,

.. math::
    \mathbf{g} = \frac{\delta E}{\delta \mathbf{m}}
                 - \left( \mathbf{m}\cdot\frac{\delta E}{\delta \mathbf{m}} \right) \mathbf{m}
               = - \mathbf{m}\times\left(\mathbf{m}\times\mathbf{H}_{\text{eff}}\right)

We will call :math:`\mathbf{g}` the tangential gradient, since it lies on the
tangent plane of the sphere at :math:`\mathbf{m}`. Its length at every site is
:math:`||\mathbf{m}\times\delta E/\delta\mathbf{m}||`, which is the residual
that has to tend to zero at a minimum, so the same quantity serves both as the
descent direction and as the stopping criterion.

Both minimisers described below step along :math:`-\mathbf{g}`. What
distinguishes them is how far they step, and this is the difficult part: the
gradient is computed from the effective field, whose magnitude has nothing to
do with the length of a spin, so the step length carries units and cannot be
guessed once and for all.


The configuration space is a manifold
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The projection above is worth stating in its proper setting, because it is
what tells us which vectors may legitimately appear inside an inner product,
and the step length rules of the next sections are built entirely out of inner
products.

A state of the system is :math:`N` unit vectors, each of which lives on a
sphere, so the set of all states is

.. math::
    \mathcal{M} = S^{2}\times S^{2}\times\cdots\times S^{2}
    \qquad (N \text{ copies})

a smooth manifold of dimension :math:`2N` embedded in
:math:`\mathbb{R}^{3N}`. It becomes a *Riemannian* manifold once each tangent
space is given an inner product, and the one we use is simply the restriction
of the Euclidean dot product of the ambient space, which is what every
``np.dot`` in the minimisers computes. The energy is a function on
:math:`\mathcal{M}`, and the minimiser has to walk on :math:`\mathcal{M}`,
not in :math:`\mathbb{R}^{3N}`.

The tangent space at :math:`\mathbf{m}` is the set of displacements
:math:`\delta\mathbf{m}` with
:math:`\mathbf{m}\cdot\delta\mathbf{m}=0` at every site. The Riemannian
gradient is defined as the unique tangent vector :math:`\mathbf{g}` with
:math:`\langle\mathbf{g},\mathbf{v}\rangle = \mathrm{d}E(\mathbf{v})` for
every tangent :math:`\mathbf{v}`, and because the metric is the restricted
Euclidean one, that vector is precisely the tangential part of
:math:`\delta E/\delta\mathbf{m}` written above. In other words, projecting
the field onto the tangent plane and taking the gradient of the energy on the
manifold are the same operation, and
:math:`-\mathbf{m}\times(\mathbf{m}\times\mathbf{H}_{\text{eff}})` is
just the projector :math:`\mathbb{1}-\mathbf{m}\mathbf{m}^{\mathsf{T}}`
written with cross products. This is what ``_project_gradient`` computes in
the Hubert class, and what ``mxmxH`` holds in the steepest descent one.

The part that is thrown away is the Lagrange multiplier of the first section:
a constraint force, not a descent direction.

The distinction is optional for the *direction* and mandatory for the *step
length*. For the direction it makes no difference, since the radial component
is annihilated by the re-normalisation and the iteration ends in the same
place. For the step length it is essential, because the Barzilai-Borwein
quotients are built from a difference of gradients taken at two different
points,

.. math::
    \mathbf{y} = \mathbf{g}_{k}-\mathbf{g}_{k-1}

and the radial components at those two points lie along two different
:math:`\mathbf{m}` vectors, so they do not cancel in the difference. They
leak into :math:`\mathbf{s}\cdot\mathbf{y}` and
:math:`\mathbf{y}\cdot\mathbf{y}` and corrupt the curvature estimate. The
secant condition the quotients solve is a statement about the Hessian *on the
manifold*, so both members of the pair have to be quantities on the manifold.

The same reasoning explains why :math:`\mathbf{s}` is the difference of the
spins after re-normalisation rather than the displacement
:math:`-\eta\mathbf{g}` that was attempted. Moving on a manifold means, in
principle, following a geodesic, and comparing tangent vectors at two
different points requires parallel transport. Fidimag does neither: it uses
the retraction "step in the ambient space, then project back", and compares
:math:`\mathbf{g}_{k}` with :math:`\mathbf{g}_{k-1}` as though they lived in
the same space. This is the standard approximation of Riemannian
optimisation [7]_, exact to first order in the step length, which is all the
secant model claims in any case. The update of the steepest descent class is
the more careful version of the same idea, being an exact rotation that stays
on the sphere by construction rather than leaving it and being pulled back.

The word carries the same meaning in the :doc:`nebm`, where the distance
between images is measured along the manifold rather than through the ambient
space, and where the two answers start to differ once the images are far
apart.


The Hubert minimiser
--------------------

This is the ``hubert_minimiser`` driver, based on the works of Berkov [1]_,
[2]_, and implemented in MERRILL [3]_. The spins are updated as

.. math::
    \mathbf{m}_{\text{new}} = \mathbf{m}_{\text{old}} - \eta \eta_{s} \mathbf{S}

where :math:`\mathbf{S}` is a Polak-Ribière conjugate direction built from the
gradient, :math:`\eta` is a scaling factor that the algorithm updates as it
goes, and :math:`\eta_{s}` is a fixed factor given by the user. The energy of
the last :math:`t` steps is kept in a trailing array, and the algorithm creeps
with a fixed :math:`\eta` for a number of steps: if the energy decreases,
:math:`\eta` is increased to accelerate the descent, and if it increases,
:math:`\eta` is decreased, with a minimum value below which the minimisation is
restarted.

The parameter that has to be tuned here is :math:`\eta_{s}`, i.e. the
``eta_scale`` argument, since it is what converts the units of the field into a
spin displacement. This is more delicate than it looks. For a one dimensional
domain wall in a micromagnetic sample, ``eta_scale = 1e-6`` and ``1e-7`` relax
the system in around 600 evaluations of the effective field, while ``1e-5``
does not converge at all, and for an atomistic system, where the field is in
different units, the useful range is somewhere else entirely.


Barzilai-Borwein step lengths
-----------------------------

We can avoid choosing the step length by estimating it from the curvature that
the system has already shown us. If we call

.. math::
    \mathbf{s} = \mathbf{m}_{k} - \mathbf{m}_{k-1} \qquad
    \mathbf{y} = \mathbf{g}_{k} - \mathbf{g}_{k-1}

the two Barzilai-Borwein quotients [4]_ are

.. math::
    \eta_{\text{BB1}} = \frac{\mathbf{s}\cdot\mathbf{s}}{\mathbf{s}\cdot\mathbf{y}}
    \qquad
    \eta_{\text{BB2}} = \frac{\mathbf{s}\cdot\mathbf{y}}{\mathbf{y}\cdot\mathbf{y}}

Both are inverse Rayleigh quotients of the Hessian along the step just taken,
i.e. they are the step length that a secant approximation of the curvature
would suggest. The two are not equivalent, BB1 tends to overshoot and BB2 to
undershoot, and alternating them is more robust than using either alone, which
is what we do by default.

The relevant observation for us is one of units. The quotients have the units
of :math:`[\mathbf{m}]/[\delta E/\delta\mathbf{m}]`, which is precisely the
product :math:`\eta\eta_{s}` of the previous section, so the conversion that
``eta_scale`` was providing by hand is now estimated at every step. Equivalently,
if the gradient is rescaled by a constant the quotients rescale by its inverse
and the update is unchanged, i.e. the method is invariant under the freedom that
``eta_scale`` was parametrising. This is why no equivalent parameter appears in
the BB path, and it is also why the same code works for the micromagnetic and
the atomistic classes without changes.

Note that :math:`\mathbf{s}` is the difference of the spins *after* they have
been re-normalised, and not the step :math:`-\eta\mathbf{g}` that we attempted.
This is deliberate: it is the step that was actually taken, and it is what the
theory of the projected version of the method requires.


Non-monotone acceptance
^^^^^^^^^^^^^^^^^^^^^^^

Barzilai-Borwein steps do not decrease the energy at every iteration. This is
not a defect but the reason they are fast, and it means that the creep
strategy of rejecting any step that increases the energy cannot be used here,
since it would discard exactly the behaviour we are after. What we ask instead
is that the energy decreases with respect to the largest of the last
:math:`t` energies, which are already being stored in the trailing array,

.. math::
    E(\mathbf{m}_{\text{new}}) \leq \max_{0\leq j < t} E_{j}
                                    - \gamma \lambda ||\mathbf{g}||^{2}

with :math:`\gamma` a small constant. If a trial step fails this test we
backtrack, i.e. we shorten it and try again from the same configuration. This
is the non-monotone line search of Grippo, Lampariello and Lucidi, and together
with the BB step lengths and the re-normalisation of the spins, which plays the
role of the projection onto the constraint set, it is the spectral projected
gradient method of Birgin, Martínez and Raydan [5]_. Setting ``nTrail = 1``
recovers a monotone line search.

There is a subtlety in the sufficient decrease term above, which compares a
gradient norm against an energy. In Fidimag the energy of the minimiser is
divided by the ``energyScale`` attribute while the gradient is the raw effective
field, so the two do not share units, and the missing constant is not the same
for the micromagnetic and the atomistic classes. Rather than hard-coding it, we
calibrate it from the decrease that the accepted steps actually produce, which
keeps the criterion meaningful without having to know the units.

A trust region completes the algorithm: no spin is allowed to move further than
``maxDeltaM`` in a single step, in units where the spin length is one. This
keeps the re-normalisation an accurate projection, and it also provides the
first step length, before any secant pair is available.

The BB path is selected with the ``stepControl`` argument, and takes no
``eta_scale``::

    sim = fidimag.micro.Sim(mesh, driver='hubert_minimiser')
    ...
    sim.driver.energyScale = Kd
    sim.driver.minimise(stepControl='BB', stopping_dE=1e-14, mXgradE_tol=1e-3)

The default is ``stepControl='hubert'``, which is the creep algorithm
described above. As a reference, the number of evaluations of the effective
field needed to reach the same minimum, with ``eta_scale`` tuned for each
system in the creep case:

+---------------------------------+------------------+---------------+
| System                          | creep            | BB            |
+=================================+==================+===============+
| 1D domain wall                  | 649              | 139           |
+---------------------------------+------------------+---------------+
| Skyrmion with demagnetising     | not converged    | 950           |
| field                           | in 6000          |               |
+---------------------------------+------------------+---------------+
| Atomistic skyrmion              | 1089             | 189           |
+---------------------------------+------------------+---------------+

The cost per evaluation is essentially the same for both, so these numbers
translate directly into computing time.


The steepest descent minimiser
------------------------------

The ``steepest_descent`` driver implements the algorithm of Exl et al. [6]_, in
which the spins are updated with

.. math::
    \mathbf{m}_{i+1} = \frac{\left(4 - \tau^{2} \mathbf{A}^{2}\right)\mathbf{m}_{i}
                             - 4\tau\, \mathbf{m}_{i}\times\mathbf{m}_{i}\times\mathbf{H}}
                            {4 + \tau^{2}\mathbf{A}^{2}}

with :math:`\mathbf{A}=\mathbf{m}_{i}\times\mathbf{H}` and :math:`\tau` a
fictitious time step. This update preserves the length of the spins by
construction, so it is a more natural way of moving on the sphere than stepping
along the tangent plane and re-normalising, and to first order in :math:`\tau`
the two agree.

It is worth noticing that :math:`\mathbf{m}\times\mathbf{m}\times\mathbf{H}` is,
up to a sign, the tangential gradient of the first section, and that the
:math:`\tau` of this method is chosen with the same Barzilai-Borwein quotients
described above, alternating between them. In other words the two minimisers
differ less than it appears: they use the same descent direction and the same
step length rule, and the differences are the update above, and the fact that
the steepest descent class does not evaluate the energy at all, stopping
instead when no spin moves further than ``stopping_dm`` in a step.

Not evaluating the energy makes each step cheaper, but it means that nothing
prevents a step from overshooting: the quotients are only an estimate, and when
the curvature they are extrapolating from is not representative the iteration
can end up far from where it should be. This is controlled by ``tmax``, which
bounds :math:`\tau` from above and which is set conservatively by default.


When the quotient carries no information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both quotients divide by a measure of the curvature along the step just taken,
:math:`\mathbf{s}\cdot\mathbf{y}` for BB1 and
:math:`\mathbf{y}\cdot\mathbf{y}` for BB2, and neither is guaranteed to be
usable. The code therefore reads

.. code-block:: c

    if (den == 0.0) {
        res = tmax;
    } else {
        res = num / den;
    }

A vanishing denominator says that the secant model sees a *flat* direction:
the gradient did not change over the step, so the quadratic it is fitting has
its minimum infinitely far away and there is no step to compute. In practice
this branch is the first iteration, where ``spin_last`` is still ``spin``, so
both :math:`\mathbf{s}` and :math:`\mathbf{y}` vanish identically and no
secant pair exists yet: Barzilai-Borwein needs two points and only one has
been visited. Falling back to the largest permitted step is a heuristic rather
than a result, taken from the MicroMagnum implementation, and it rests on the
observation that the next iteration will have a real secant pair to correct
it. The Hubert minimiser reaches the same point differently, using the
trust-region step :math:`\Delta m_{\text{max}}/\max||\mathbf{g}||` of its
``maxDeltaM`` argument, which is a length rather than a ceiling.

The exact comparison against zero only catches that degenerate case. A merely
*small* denominator is far more common and gives a large but finite
:math:`\tau`, which is handled not here but by the clamp on the following
line, ``tau = fmax(fmin(res, tmax), tmin)``. The two are the same safeguard
written twice.

A *negative* quotient is a third case, and a different one. It means negative
curvature along the last step, where the secant model is not uninformative but
wrong, and using it would reverse the sign of the update, driving the
iteration against :math:`-\mathbf{m}\times\mathbf{m}\times\mathbf{H}` and
up the energy. MicroMagnum keeps the sign; we fall back to ``tmax`` as above.
This is why the minimiser used to stall in configurations that were not
minima. The BB path of the Hubert class makes the same test as
:math:`\mathbf{s}\cdot\mathbf{y} > 0`, which is the standard condition for
accepting a secant pair.


The energy guard
^^^^^^^^^^^^^^^^

The same acceptance test used in the BB path can be enabled here with the
``energy_guard`` argument::

    sim.driver.minimise(stopping_dm=1e-9, energy_guard=True)

The energy is obtained from the same pass over the interactions as the
effective field, so the guard does not add any field computation, only the
occasional rejected step.

Whether this is worth its cost depends on ``tmax``. At the default value of
0.1 the steps are short enough that the test essentially never fails, the
guarded and unguarded runs give the same answer, and the guard is only
overhead, which is why it is off by default. It becomes useful when we raise
``tmax`` to take longer steps, which is what one would want to do to gain
speed. For a one dimensional domain wall with ``tmax = 10``, over eight
repetitions, the unguarded iteration took between 210 and 3000 steps and ended
in configurations whose scaled energies ranged from the correct 0.196 to 77,
while the guarded one took 130 steps with 12 rejected trial steps in every
single repetition. For comparison, the same wall at the default ceiling takes
around 305 steps.

The spread of the unguarded runs is not physics. The Barzilai-Borwein quotients
are accumulated with an OpenMP reduction, whose summation order is not fixed
from one run to the next, and an overshooting step amplifies those last-bit
differences into completely different trajectories. Rejecting those steps is
what makes the calculation repeatable.

.. [1] Berkov, D. *Numerical calculation of the energy barrier distribution in
   disordered many-particle systems: The path integral method*. J. Magn. Magn.
   Mater. 186, 199–213 (1998)

.. [2] Berkov, D. V., Ramstöck, K. & Hubert, A. *Solving Micromagnetic
   Problems. Towards an Optimal Numerical Method*. Phys. Stat. Sol. (a) 137,
   207–225 (1993)

.. [3] Ó Conbhuí, P. et al. *MERRILL: Micromagnetic earth related robust
   interpreted language laboratory*. Geochem. Geophys. Geosyst. 19, 1080–1106
   (2018)

.. [4] Barzilai, J. & Borwein, J. M. *Two-point step size gradient methods*.
   IMA J. Numer. Anal. 8, 141–148 (1988)

.. [5] Birgin, E. G., Martínez, J. M. & Raydan, M. *Nonmonotone spectral
   projected gradient methods on convex sets*. SIAM J. Optim. 10, 1196–1211
   (2000)

.. [6] Exl, L. et al. *LaBonte's method revisited: An effective steepest
   descent method for micromagnetic energy minimization*. J. Appl. Phys. 115,
   17D118 (2014)

.. [7] Absil, P.-A., Mahony, R. & Sepulchre, R. *Optimization Algorithms on
   Matrix Manifolds*. Princeton University Press (2008)
