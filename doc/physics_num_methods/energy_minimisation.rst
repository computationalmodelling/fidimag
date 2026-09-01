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

The same acceptance test used in the BB path can be turned on here with the
``energy_guard`` argument::

    sim.driver.minimise(stopping_dm=1e-9, energy_guard=True)

The energy is obtained from the same pass over the interactions as the
effective field, so the guard costs no field computation, only the occasional
rejected step.

What it buys is a higher ``tmax``. The ceiling decides how fast the method
can be, and it cannot be raised far without the guard. On the one dimensional
domain wall of ``tests/test_steepest_descent.py``, counting every evaluation
of the effective field and the mean error of the relaxed wall profile:

=============  ===================  =================
``tmax``       guard off            guard on
=============  ===================  =================
0.1            346, err 0.0018      350, err 0.0018
1              142, err 0.0018      173, err 0.0018
3              157, err 0.0018      177, err 0.0018
10             223, **err 0.90**    157, err 0.0018
=============  ===================  =================

Unguarded, the wall collapses somewhere between 3 and 10. The default ceiling
is 1, an order of magnitude below that, which on the standard problem 4
s-state takes 2498 evaluations against the 12318 of the 0.1 it used to be. A
ceiling of 3 takes 1478 and of 10, guarded, would take fewer still.

The guard is nonetheless off by default, because of how it interacts with the
stopping criterion. ``stopping_dm`` bounds how far a spin moves in a step, and
a rejected trial step is a short one, so backtracking can satisfy the
criterion while the residual is still large. On the s-state, the guarded
iteration stopped at :math:`4\times10^{-6}` A/m with ``tmax = 3`` and at
:math:`5\times10^{-6}` with ``tmax = 10``, never reaching the
:math:`10^{-6}` that the unguarded runs reach; the configuration is right, its
energy agreeing to ten digits, but the convergence is looser than asked for.
A criterion on the torque rather than on the displacement would fix this, and
is what the guard needs before it can be turned on by default.

The spread of the unguarded runs is not physics. The effective field is not
bit-reproducible from one run to the next, for the reason given under
:ref:`reproducibility` below, and an overshooting step amplifies those last-bit
differences into completely different trajectories. Rejecting those steps is
what makes the calculation repeatable.

How this compares with OOMMF
----------------------------

OOMMF minimises with ``Oxs_CGEvolve``, a nonlinear conjugate gradient with a
proper line search: it brackets a minimum along the search direction and
refines it by interpolation, using the directional derivative as well as the
energy. That is a more careful step rule than anything here, so the comparison
worth making is not per iteration but per *effective field evaluation*, which
is what dominates the cost and which counts a line search honestly.

The problem below is the standard problem 4 film, 500 x 125 x 3 nm of
permalloy with :math:`M_{s} = 8\times10^{5}` A/m and
:math:`A = 1.3\times10^{-11}` J/m, exchange and demagnetising field, relaxed
from a uniform :math:`\mathbf{m}=(1,1,1)` to the s-state. Both codes report the
same convergence measure, the largest
:math:`||\mathbf{m}\times(\mathbf{m}\times\mathbf{H})||` over the mesh, in
A/m; OOMMF calls it ``Max mxHxm`` and publishes the evaluation count as
``Energy calc count``. All three minimisers reach the same s-state, at the same
energy to ten digits.

Evaluations needed to first reach a given torque, on 2500 cells of 5 nm:

===================  =====  =====  =====  =====  =====  =====
torque (A/m)          1e-1   1e-2   1e-3   1e-4   1e-5   1e-6
===================  =====  =====  =====  =====  =====  =====
OOMMF CG               343    409    457    507    551    614
Fidimag BB             226    259    305    321    377    413
Fidimag SD            1281   1522   1767   2003   2248   2498
===================  =====  =====  =====  =====  =====  =====

and on 10000 cells of 2.5 nm:

===================  =====  =====  =====  =====  =====  =====
torque (A/m)          1e-1   1e-2   1e-3   1e-4   1e-5   1e-6
===================  =====  =====  =====  =====  =====  =====
OOMMF CG               728    826    936   1050   1156   1278
Fidimag BB             731    792    894    929    989   1038
Fidimag SD            2528   3013   3507   3991   4490   5000
===================  =====  =====  =====  =====  =====  =====

The Barzilai-Borwein path of the Hubert class needs fewer evaluations than the
conjugate gradient at every tolerance, by about a third on the coarse mesh and
a fifth on the fine one, and it keeps going well past the last row, to
:math:`2\times10^{-10}` A/m on the coarse mesh. It carries a conjugate
direction, a trust region and a non-monotone line search that the steepest
descent does not.

The steepest descent is the slowest of the three here, by a factor of four
against the conjugate gradient at its default step ceiling, though far from
the factor of twenty that the old ceiling of ``tmax = 0.1`` cost: that took
12318 evaluations on the coarse mesh where the current default of 1 takes
2498, and a ceiling of 3 takes 1478 and of 10 takes 886. It implements the
method of Exl et al. [6]_, which is also what MuMax3 minimises with, so what
is being compared is a step length rule and the ceiling put on it rather than
two codes, and the ceiling here is deliberately conservative. The next section
is about choosing it.

Wall clock does not change the picture: measured per evaluation, 0.23 ms
against 0.80 ms for OOMMF on the coarse mesh and 0.99 ms against 1.07 ms on
the fine one, both codes being bound by the same transform.

The counts move by about ten per cent from one run to the next, for the reason
under :ref:`reproducibility`.


The energy precision floor
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The numbers above are what the Barzilai-Borwein path does now. Until recently
it stopped at :math:`1.5\times10^{-6}` A/m on the coarse mesh and
:math:`5\times10^{-4}` A/m on the fine one, with
``Could not decrease the energy along the gradient``, and the reason is worth
recording because it was not the step rule.

The acceptance test compares the energy of a trial step against the trailing
window. Near a minimum the decrease it has to detect becomes smaller than the
gap between neighbouring doubles at the size of the total energy:

.. math::
    E \approx 6.31\times10^{-19}\,\text{J} \qquad
    \varepsilon\,|E| \approx 1.40\times10^{-34}\,\text{J} \qquad
    \Delta E \approx 1.93\times10^{-35}\,\text{J}

The decrease is a seventh of the smallest representable difference. Every
trial step then looks like a failure, the restarts are used up and the
iteration gives up, at a torque that depends on the last bits of the field and
so is not even repeatable.

The cure is to sum the change rather than subtract two totals. Every
interaction reports the energy of each cell in ``energy``, in joules, with
``total_energy`` their sum, so the change between two configurations can be
taken site by site,

.. math::
    \Delta E = \sum_{i} \left( E_{i}^{\text{new}} - E_{i}^{\text{ref}} \right)

Each term is a difference of two numbers of the size of a single cell energy,
not of the whole sample, so nothing large is subtracted and the result is
accurate relative to :math:`\Delta E` rather than to :math:`E`. This is what
``Oxs_CGEvolve`` does, and the reason the interaction classes were made to
agree on what ``energy`` means.

Scaling the energy does not help, and it is worth saying why, since
``energyScale`` looks like exactly the knob for the job. What decides whether
the difference survives is :math:`\Delta E/E` against :math:`\varepsilon`, and
dividing both by a constant leaves that ratio alone: here
:math:`\Delta E/E = 3.1\times10^{-17}` is below
:math:`\varepsilon = 2.2\times10^{-16}`, at any scale. The subtraction
``(E - dE) - E`` returns exactly zero for :math:`E` of
:math:`6\times10^{-19}`, of one, and of :math:`10^{12}` alike. Nor is the
algorithm supposed to care: the gradient is never divided by ``energyScale``,
and ``gradScale`` is calibrated from the decrease, so the scale cancels out of
the acceptance test in exact arithmetic.

What a different scale does change is the rounding, and so the sequence of
accepted steps and the point at which the iteration gives up. On the problem
above, five values of ``energyScale`` that all put the total within an order
of magnitude of one stopped at torques of
:math:`4.5\times10^{-8}`, :math:`2.9\times10^{-6}`,
:math:`9.5\times10^{-6}`, :math:`4.5\times10^{-6}` and
:math:`3.1\times10^{-7}` A/m. That is a spread of two hundred with no trend:
a good value is luck rather than a cure, and there is no way to recognise one
in advance.

``_minimise_BB`` therefore carries :math:`E-E_{0}`, accumulated from the
accepted steps, instead of reading totals.

With that in place the same problem reaches :math:`2.4\times10^{-10}` A/m
instead of :math:`1.5\times10^{-6}`, in fewer evaluations than OOMMF needs to
reach :math:`10^{-6}`, and it returns the same answer to the last bit for
every ``energyScale`` from one to :math:`10^{-31}`, the scale having dropped
out of a subtraction that is no longer performed. The steepest descent never
had the problem, since it does not evaluate the energy at all.


A second route, not taken
""""""""""""""""""""""""""

There is another way to the same number, worth recording because it needs no
per cell energy at all, only the effective field, which the minimiser has
already computed. Writing the energy of the interactions that are quadratic in
the magnetisation as
:math:`E=-\tfrac{1}{2}\sum_{i}w_{i}\,\vec{m}_{i}\cdot\vec{H}_{i}`, with
:math:`w_{i}=\mu_{0}VM_{s,i}` in the micromagnetic case and
:math:`\mu_{s,i}` in the atomistic one, and labelling the reference and the
trial configurations :math:`r` and :math:`n`, the field is linear in
:math:`\vec{m}` with a symmetric operator, so
:math:`\vec{m}_{n}\cdot\vec{H}_{r}=\vec{m}_{r}\cdot\vec{H}_{n}` and the two
mixed terms of
:math:`(\vec{m}_{n}-\vec{m}_{r})\cdot(\vec{H}_{n}+\vec{H}_{r})` cancel,
leaving

.. math::
    \Delta E = -\frac{1}{2}\sum_{i} w_{i}\,
               (\vec{m}_{n}-\vec{m}_{r})_{i}\cdot
               (\vec{H}_{n}+\vec{H}_{r})_{i}

a trapezoid that the reciprocity makes exact rather than first order, so it
holds at any step length. It carries the small displacement in every term, so
it does not cancel either, and a constant Zeeman field is covered by the same
expression because its two fields are equal.

Measured on the problem above it reached the same convergence as the sum over
sites, to within the run to run spread. It was not kept because it is worse on
three counts: it is exact only while the energy is quadratic, a cubic
anisotropy being the case to watch, where it degrades to the usual first order
estimate; it needs the field at both ends of the step, :math:`3n` stored
against :math:`n`; and it gives the change only up to a constant, which
differs between the micromagnetic and atomistic classes and has to be
calibrated. It is recorded here because none of those objections is fatal, and
it is the route to take if per cell energies are ever unavailable.


What OOMMF does about it
""""""""""""""""""""""""

None of this is new, and ``Oxs_CGEvolve`` is worth reading on the point. It
solved the same problem long ago, and in the same way, summing the *per cell*
energy differences against the best state so far,

.. code-block:: c

    work_etemp.Accum((stenergy[j] - sbenergy[j]));

The comment nearby attributes the analysis to notes from 2002. Fidimag could
not do this until its interaction classes were made to agree on what
``energy`` means: it was an energy density in most of the micromagnetic ones,
already weighted in ``micro/demag.py``, and never written at all by
``micro/zeeman.py``.

OOMMF is also more careful than Fidimag is. It
accumulates in extended precision rather than in doubles, it carries an error
bar with every bracket point,
``E_error_estimate = fabs(relenergy)*OC_REAL8m_EPSILON*8`` plus a density
term, and ``EstimateEnergySlack`` turns those into a threshold below which,
in its own words, two energies "should be considered equal". When a bracket
falls inside that slack ``BadPrecisionTest`` stops trusting the energies and
brackets on the slope instead, which comes from the field and does not
cancel.

The sum used here has a floor of its own, since each site difference carries
an error of order :math:`\varepsilon` times its own site energy, so the total
is good to about :math:`\varepsilon E/\sqrt{N}` rather than
:math:`\varepsilon E`. That was not reached on any problem tried here, but a
compensated sum with a fallback on the slope, as OOMMF does, is the more
robust arrangement if it ever matters.

Note that ``stopping_dE`` is an absolute energy, so its default of
:math:`10^{-6}` is meaningless against a micromagnetic energy of
:math:`10^{-17}` J and will stop the iteration immediately. Set
``energyScale`` to bring the energy near unity, or pass a sensible
``stopping_dE``, and lean on ``mXgradE_tol`` instead. This is long standing
and unrelated to the above.

.. _reproducibility:

Reproducibility
^^^^^^^^^^^^^^^

Repeating one of these minimisations does not always give the same number of
steps or the same final torque. The minimisers are deterministic, and contain
no random numbers; what varies is the demagnetising field, in its last bits.
Over six identical runs on the mesh above the largest difference between any
two evaluations of the same field was :math:`2.9\times10^{-10}` A/m against a
field of :math:`4.5\times10^{5}` A/m, a relative :math:`6.4\times10^{-16}`, or
three units in the last place.

The cause is the FFTW planner. The demagnetising transforms are created with
``FFTW_MEASURE``, which benchmarks several algorithms at plan time and keeps
whichever ran fastest, so timing noise on a loaded machine selects a different
algorithm and therefore a different summation order. Each is equally correct.
The plan is cached for the lifetime of the process, which is why repeated runs
inside one process agree exactly while the first run of a fresh process may
not.

Given a fixed plan the whole minimisation is reproducible to the last digit,
so what one actually observes is a handful of discrete trajectories rather
than a continuous spread. On the coarse mesh above, three plans gave 848, 1906
and 844 evaluations and final torques of
:math:`1.5`, :math:`4.2` and :math:`2.8 \times 10^{-6}` A/m. Note what does
*not* vary: every run finds the same s-state at the same energy, and the work
needed to reach any useful tolerance agrees to within about 20%. What varies is
where the iteration gives up, which is the energy precision floor of the
previous section being decided by those last bits.

For a calculation that must be repeatable to the bit, build against a planner
that does not measure, or arrange for fixed FFTW wisdom to be loaded. For a
scientific result the distinction does not usually matter, since the minimum
is well defined and the difference is far below the accuracy at which a
micromagnetic answer is compared with an analytical prediction.


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
