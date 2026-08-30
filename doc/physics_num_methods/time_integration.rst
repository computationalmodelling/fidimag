Time integration
================

The LLG equation, as written in the previous section, conserves the length of
the magnetisation exactly: the right hand side is a cross product with
:math:`\vec{m}`, so it is perpendicular to :math:`\vec{m}` and cannot change
:math:`|\vec{m}|`. A numerical integrator does not inherit that property. Over
a step of finite size the length drifts, and the drift accumulates.

Fidimag deals with this in two ways, which are independent of each other and
can be used together. The first is built into the equation and is described
below; the second is to project the spins back onto the unit sphere after
every step, which only the explicit integrators can do.


The stabilisation term
----------------------

The right hand side that Fidimag integrates is not the LLG equation alone, but

.. math::
    \frac{\partial \vec{m}}{\partial t} \longrightarrow
    \frac{\partial \vec{m}}{\partial t} + c \left( 1 - \vec{m}^{2} \right)
    \vec{m}

There is no physics in the added term. It is a numerical stabilisation of the
constraint :math:`|\vec{m}| = 1`, of the kind introduced by Baumgarte for
constrained mechanics, and its effect is easiest to see in the length itself.
Writing :math:`\rho = |\vec{m}|`, and using that the LLG part contributes
nothing to :math:`\partial \rho / \partial t`,

.. math::
    \frac{\partial \rho}{\partial t} = c \, \rho \left( 1 - \rho^{2} \right)

which is the normal form of the supercritical pitchfork bifurcation:
:math:`\rho = 0` is an unstable fixed point and :math:`\rho = 1` a stable one,
so a length error is pulled back towards one rather than left to grow.
Linearising about :math:`\rho = 1`, with :math:`\rho = 1 + \epsilon`,

.. math::
    \frac{\partial \epsilon}{\partial t} = -2 c \, \epsilon

so :math:`c` is nothing more than the rate at which length errors decay, with
a time constant of :math:`1 / (2c)`. That also fixes the two limits between
which it has to sit: too small and the error is not removed within a step, too
large and the equation is stiffened by an eigenvalue :math:`-2c`, which costs
step size.

The coefficient is the ``default_c`` attribute of the driver:

``default_c < 0``
    Use :math:`c = 6 \, |\partial \vec{m} / \partial t|`, so that length errors
    relax twelve times faster than the magnetisation rotates. Writing it
    against :math:`|\partial \vec{m} / \partial t|` rather than as an absolute
    rate is what makes it scale free. This is the default of the atomistic
    driver.

``default_c > 0``
    Use the value as it stands, as a rate in inverse seconds. The
    micromagnetic driver defaults to :math:`10^{11}`, a time constant of 5 ps,
    against a precession period of about 35 ps at :math:`\gamma M_s`.

``default_c = 0``
    Turn the term off.


How much does the value matter?
-------------------------------

Less than its arbitrary appearance suggests. Measured on standard problem 4
against the OOMMF reference in ``examples/micromagnetic/std4``, at
``rtol = atol = 1e-10``.

Throughout this section, *deviation* means the largest difference, taken over
all sample times and over the three components, between the spatially averaged
magnetisation :math:`\langle \vec{m} \rangle` and the reference value at the
same time. It is a worst case over the trajectory, not an average along it, and
it is a difference of the averaged magnetisation, which is what the standard
problem 4 reference data reports; a per cell error would be larger. The drift
column is the largest :math:`||\vec{m}_i| - 1|` over every cell and every
sample time.

+-------------------+---------------------+----------------------+
| ``default_c``     | max :math:`||m|-1|` | deviation from OOMMF |
+===================+=====================+======================+
| 0 (off)           | 5.80e-09            | 3.20e-05             |
+-------------------+---------------------+----------------------+
| 1e9               | 5.74e-09            | 3.20e-05             |
+-------------------+---------------------+----------------------+
| 1e11 (default)    | 2.28e-09            | 3.20e-05             |
+-------------------+---------------------+----------------------+
| 1e13              | 5.21e-11            | 3.20e-05             |
+-------------------+---------------------+----------------------+
| :math:`6|dm/dt|`  | 5.09e-10            | 3.20e-05             |
+-------------------+---------------------+----------------------+

The length drift falls as :math:`1/c`, exactly as the time constant above
predicts, and at :math:`10^{13}` the step count rises by 15%, which is the
stiffening. The accuracy of the solution, however, does not move at all: the
deviation is the same to three figures for every value, including with the
term switched off, because at this tolerance the integrator's own error
control is already holding the drift to :math:`6 \times 10^{-9}` unaided.

The term earns its keep at loose tolerances. At ``rtol = atol = 1e-6`` with
``cvode_bdf`` it nearly halves the error, 6.36e-05 with the term off against
3.58e-05 at :math:`10^{11}`. For the explicit ARKODE methods it does nothing
useful at any tolerance, and at 1e-8 it makes the drift worse, 3.0e-10 with
the term off against 1.8e-08 with :math:`6|dm/dt|`, since the term is itself a
perturbation of the solution. Those integrators can project instead.

So the value is not a knob to tune for accuracy. Anything from about
:math:`10^{11}` to :math:`10^{12}`, or the :math:`6|dm/dt|` rule, behaves the
same, and it only has to be large enough to remove drift within a step and
small enough not to stiffen the equation. The same term, and the same
pitchfork analysis, are in [1]_, which argues that the coefficient should be
unity once the equation is written in a suitable dimensionless form, and which
attributes the :math:`|\partial \vec{m}/\partial t|` rule used here to Sec. 2.4
of [2]_.


Choosing an integrator
----------------------

The integrator is chosen with the ``integrator`` argument of the simulation
classes. The names say which solver is used and which method it runs:

``cvode_bdf``
    CVODE with backward differentiation formulae and a Newton iteration,
    solved with restarted GMRES. The default, and the one to keep on a fine
    mesh: the stiffness of the LLG equation comes from the exchange
    interaction, whose fastest timescale grows as :math:`1/\Delta x^{2}`, so
    halving the cell size quadruples it.

``cvode_bdf_diag``
    The same, with CVODE's diagonal approximate Jacobian in place of GMRES.

``cvode_adams``
    CVODE with Adams-Moulton and a fixed point iteration: the non-stiff arm of
    the same solver, with no linear solve at all.

``arkode_dopri5``, ``arkode_rkf45``
    Explicit Runge-Kutta 5(4) pairs from ARKODE. ``arkode_dopri5`` uses the
    Dormand and Prince tableau, which is the method OOMMF calls ``rkf54m``, so
    the two codes can be compared directly. Note that OOMMF's own default,
    ``rkf54``, is RK5(4)7FC, another member of the same family rather than the
    Fehlberg tableau its name suggests; the genuine Fehlberg one is
    ``arkode_rkf45``.

``arkode_dopri5_normalised``, ``arkode_rkf45_normalised``
    The same, rescaling every spin to unit length after each accepted step,
    through ARKODE's post-step hook. This is the projection mentioned above,
    and it is independent of the stabilisation term: either, both or neither
    can be used. CVODE has no equivalent hook.

``euler``, ``rk4``
    Fixed step, implemented in Fidimag rather than taken from a library, and
    mostly useful for debugging.

Any of the SUNDIALS names also takes an ``_openmp`` suffix, which uses the
OpenMP N_Vector. That threads the integrator's own vector arithmetic only; the
effective field is already parallel whichever vector is used.

The names used before Fidimag 4.0, ``sundials``, ``sundials_diag``,
``sundials_openmp`` and ``sundials_diag_openmp``, still work and raise a
deprecation warning naming their replacement.

On standard problem 4, at ``rtol = atol = 1e-10`` and with four threads, all of
these agree with each other to within :math:`2\times10^{-8}`, in the same sense
as above, and put the crossing of :math:`\langle m_x \rangle = 0` in the same
place. What differs is the cost:

+---------------------------------+-----------+---------+
| integrator                      | wall time | steps   |
+=================================+===========+=========+
| ``cvode_bdf``                   | 8.63 s    | 6196    |
+---------------------------------+-----------+---------+
| ``cvode_bdf_openmp``            | 6.67 s    | 6196    |
+---------------------------------+-----------+---------+
| ``arkode_dopri5``               | 4.05 s    | 2417    |
+---------------------------------+-----------+---------+
| ``arkode_dopri5_openmp``        | 3.74 s    | 2417    |
+---------------------------------+-----------+---------+

The explicit methods evaluate the effective field more often, but a step costs
only the stages of the tableau, with no Newton iteration and no Krylov solve.
That ordering is particular to a problem of this size and should not be taken
as general: on a finer mesh the exchange stiffness grows and the implicit
methods recover their advantage. To find an equilibrium state rather than to
follow the dynamics, use the minimisers described in
:doc:`energy_minimisation` instead of integrating at all.

The analytical Jacobian
-----------------------

The implicit methods solve a linear system at every Newton iteration, and
GMRES needs only the product of the Jacobian with a vector, never the matrix.
By default CVODE forms that product by a difference quotient of the right hand
side. Passing ``use_jac=True`` to the simulation supplies it analytically
instead, through ``sundials_jtimes``, which differentiates

.. math::
    \frac{\partial \vec{m}}{\partial t} = \frac{-\gamma}{1 + \alpha^{2}}
    \left( \vec{m}\times\vec{H}_{\perp}
            - \alpha \vec{H}_{\perp} \right)
    + c \left( 1 - \vec{m}^{2} \right) \vec{m}

in a direction :math:`\vec{m}'`, with
:math:`\vec{H}_{\perp} = (\vec{m}\cdot\vec{m})\vec{h} -
(\vec{m}\cdot\vec{h})\vec{m}`. Because the effective field is linear in
:math:`\vec{m}` for exchange and for the demagnetising field, its derivative
in that direction is just the field evaluated at :math:`\vec{m}'`, which is
what ``compute_effective_field_jac`` computes.

This did not work for a decade, and the failure was not obvious. It looked
like a hang: the product being wrong, GMRES never converged and spent its
whole iteration budget on every solve, half a million evaluations for a single
picosecond of simulated time. There were three faults, all now fixed. The
routine was given ``dm/dt`` where it wanted the effective field. It indexed
the damping and the pinned sites per component rather than per site, reading
past the end of both arrays. And the precession term omitted the factors of
:math:`\vec{m}\cdot\vec{m}` and a whole term: the derivative of
:math:`\vec{m}\times\vec{H}_{\perp}` is

.. math::
    (\vec{m}\cdot\vec{m})(\vec{m}'\times\vec{h})
    + (\vec{m}\cdot\vec{m})(\vec{m}\times\vec{h}')
    + 2(\vec{m}\cdot\vec{m}')(\vec{m}\times\vec{h})

the two :math:`(\vec{m}\cdot\vec{h})(\vec{m}\times\vec{m}')` terms having
cancelled. Only the first two were present, and without their prefactors,
which is correct only where :math:`|\vec{m}|` is exactly one.

The product also needs the effective field at :math:`\vec{m}`, and
:math:`\vec{m}` does not change through the GMRES iterations of one Newton
solve, so that field is computed once for the whole batch and reused, by
``effective_field_at``. The cache is keyed on the time and on the state
itself, and is dropped at the start of every ``run_until`` in case the
interactions have been changed in between. It is what makes the analytical
product worth using: on a 30x30x10 mesh of 1 nm cells with exchange and the
demagnetising field, 20 ps of dynamics take 2.21 s with the difference
quotient and 2.42 s with the analytical product recomputing the field every
call, but 1.70 s once it is reused, for the same 361 steps. The field is
evaluated 1798, 1814 and 913 times respectively.

``tests/test_jacobian.py`` compares the product against a finite difference
of the right hand side, which is the one reference that cannot be wrong in
the same way, and it agrees to about :math:`10^{-8}`.

There is no preconditioner
^^^^^^^^^^^^^^^^^^^^^^^^^^

GMRES runs unpreconditioned. There used to be a preconditioner, and it was
presumably meant to speed up the convergence that the analytical Jacobian was
supposed to improve, but it never did anything: its solve was
:math:`\vec{z} = \vec{r}`, the identity, so it only added a copy per
iteration, and its setup function was attached with the wrong signature, so
CVODE called it with every argument shifted by one. Both are removed. A real
preconditioner is still worth having, and would be the next thing to try for
the stiff case.

Two known gaps
^^^^^^^^^^^^^^

Both are recorded in the tests rather than hidden:

- Interactions are included in the Jacobian only if they set ``jac = True``.
  Nothing has to be derived for the linear ones: exchange, the DMI and the
  demagnetising field are all linear in :math:`\vec{m}`, so the derivative in
  a direction is simply the field evaluated at that direction, which is what
  ``compute_effective_field_jac`` computes. The demagnetising field is linear
  to machine precision, :math:`3\times10^{-16}` on a check of
  :math:`\vec{H}(a\vec{m}_1 + b\vec{m}_2)`. A constant Zeeman field is
  correctly excluded, its derivative being zero.

  What the flag really trades is cost against exactness, and CVODE tolerates
  an approximate Jacobian. Exchange is local and cheap, and it is the stiff
  term, so leaving it out is expensive: on a 30x30x4 atomistic lattice,
  including it took the integration from 150 steps to 18. The demagnetising
  field is the opposite case, costing a transform per evaluation: on a
  30x30x10 mesh of 1 nm cells, the run above takes 1.70 s with it left out
  against 2.84 s with it in, for the same 361 steps. So exchange is in for
  both models and the demagnetising field is out, which is a deliberate trade
  rather than an oversight.
- The stabilisation term is differentiated only when ``default_c`` is
  positive. A negative value selects :math:`c = 6\,|dm/dt|`, which the right
  hand side applies and the Jacobian then ignores. That is the default of the
  atomistic driver.

.. [1] Botha, A. E. *Stabilisation of the Landau-Lifshitz-Gilbert equation for
   numerical solution via standard methods*. Sci. Rep. 15, 15775 (2025)

.. [2] Chernyshenko, D. *Computational Methods in Micromagnetics*. PhD thesis,
   University of Southampton (2016)
