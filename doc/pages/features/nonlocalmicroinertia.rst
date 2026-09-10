Non-local micro-inertia (hyperbolic gradient damage)
====================================================

Theory
------

The implicit gradient-enhanced formulation regularises softening by solving a Helmholtz
equation for the non-local variable :math:`\bar\varepsilon` alongside the momentum
balance,

.. math::

   \bar\varepsilon - c\, \nabla^2 \bar\varepsilon = \tilde\varepsilon ,
   \qquad c = l^2 ,

with :math:`\tilde\varepsilon` the local driving variable and :math:`l` the internal
length. An explicit solver cannot solve that equation -- it has no time derivative to
march -- so it is made parabolic by adding a viscosity :math:`\eta`,

.. math::

   \eta\, \dot{\bar\varepsilon} + \bar\varepsilon - c\, \nabla^2 \bar\varepsilon
     = \tilde\varepsilon ,

and integrated with forward Euler. That works, and it costs dearly on a fine mesh. The
largest eigenvalue of the discrete operator is

.. math::

   \lambda_\mathrm{max} \approx \frac{1}{\eta}\left(1 + C\,\frac{l^2}{h^2}\right),

so the stable increment is

.. math::

   \Delta t \le \frac{2\eta}{1 + C\,l^2/h^2}
             \;\xrightarrow{\;h \ll l\;}\;
             \frac{2\eta\,h^2}{C\,l^2} \;\propto\; h^2 ,

where :math:`C` is a discretisation constant of order 4 in one dimension and roughly 12
to 24 for a hexahedron. **The limit falls off with the square of the element size.** Under
adaptive refinement that is the binding constraint long before the mechanical one is, and
the only knob it offers is :math:`\eta` itself -- which is the artificial delay between
the damage front and the strain concentration driving it. Buying stability with
:math:`\eta` buys it by making the regularisation lag.

Adding a **micro-inertia** :math:`m_k` breaks that coupling:

.. math::

   m_k\, \ddot{\bar\varepsilon} + \eta\, \dot{\bar\varepsilon} + \bar\varepsilon
     - c\, \nabla^2 \bar\varepsilon = \tilde\varepsilon .

The equation is now a damped wave equation -- a telegraph, or damped Klein--Gordon,
equation -- for a damage front travelling at

.. math::

   c_k = \sqrt{\frac{c}{m_k}} = \frac{l}{\sqrt{m_k}} ,

and :math:`m_k` has units of seconds squared. Its viscosity has not changed meaning; it
has changed *role*. What was the coefficient of the highest time derivative is now the
damping.

Stability
^^^^^^^^^

With a lumped micro-inertia and central differences the limit is

.. math::

   \omega_\mathrm{max} = \sqrt{\frac{1 + C\,l^2/h^2}{m_k}} ,
   \qquad
   \zeta = \frac{\eta/m_k}{2\,\omega_\mathrm{max}} ,
   \qquad
   \Delta t \le \frac{2}{\omega_\mathrm{max}}
                 \left( \sqrt{1+\zeta^2} - \zeta \right) .

For :math:`h \ll l` this is :math:`\Delta t \approx 2\sqrt{m_k}\,h / (\sqrt{C}\, l)`:
**linear in the element size**, the same scaling elastodynamics has. On a mesh coarser
than the internal length the reaction term takes over instead and the limit saturates at
:math:`2\sqrt{m_k}`, so the field can never demand a smaller increment than that
regardless of how coarse the mesh is.

Two things are easy to get wrong here and both cost a factor of two or more:

- The continuum estimate :math:`\Delta t \le h/c_k` **overestimates** the limit. The
  discrete Laplacian's eigenvalue constant contributes :math:`2/\sqrt{C} \approx 0.4` to
  :math:`0.6`, not 1.
- Damping **lowers** the central-difference limit; it never raises it.

Both are accounted for in
:cpp:func:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeCriticalTimeStepForExplicitDynamics`,
which bounds :math:`\lambda_\mathrm{max}` with a Gershgorin estimate over the assembled
operator rather than with a wave speed, and takes the minimum against the mechanical
limit.

Choosing the two parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^

They are not independent. The zeroth-order reaction term is an oscillator of frequency
:math:`\omega_0 = 1/\sqrt{m_k}`, damped by :math:`\eta`; requiring that it not ring gives

.. math::

   m_k \le \frac{\eta^2}{4} ,

and the largest admissible micro-inertia is the best one, because :math:`\Delta t` grows
with :math:`\sqrt{m_k}`. So **there is one free parameter, not two**: pick the viscosity
:math:`\eta` for the lag you are willing to accept -- exactly as before -- and take
:math:`m_k = \eta^2/4`. The resulting damage wave speed :math:`c_k = 2l/\eta` should then
be checked against the loading rate, not chosen from it.

At that choice the gain over the parabolic scheme *at the same artificial lag* is

.. math::

   \frac{\Delta t_\mathrm{hyperbolic}}{\Delta t_\mathrm{parabolic}}
     = \frac{\sqrt{1 + C\,l^2/h^2}}{2}
     \;\approx\; \frac{\sqrt{C}\,l}{2h} ,

a factor that **grows with every refinement level** -- which is the point: the finer the
mesh, the more the second-order form is worth.

Keeping the mechanical problem in charge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is a second criterion on :math:`\eta`, and it is the one that decides whether any of
the usual quasi-static-explicit tricks work at all. The two limits the element takes the
minimum of are

.. math::

   \Delta t_\mathrm{mech} \approx k_c\,\frac{h}{c_d} ,
   \qquad
   \Delta t_\mathrm{nl} \approx \frac{2 k_c}{\sqrt{C}}\,\frac{\sqrt{m_k}\,h}{l}
   \qquad (h \ll l)

and **both are linear in** :math:`h`. Their ratio therefore contains no mesh size at all,

.. math::

   \frac{\Delta t_\mathrm{nl}}{\Delta t_\mathrm{mech}}
     = \frac{\eta\,c_d}{l\,\sqrt{C}} ,

the courant number cancelling because it is applied to both. Which of the two is in charge
is thus a property of the *parameters*, not of the discretisation: fix it once and it holds
at every refinement level. Requiring the mechanical limit to be the binding one gives

.. math::

   \eta \;\ge\; \frac{\sqrt{C}\,l}{c_d}
   \qquad\Longleftrightarrow\qquad
   m_k \;\ge\; \frac{C}{4}\left(\frac{l}{c_d}\right)^{\!2} ,

i.e. **the micro-inertia must exceed the square of the time a mechanical wave needs to
cross one non-local length** -- the only natural time scale the coupled problem offers.
Together with :math:`m_k \le \eta^2/4` this is a window rather than a conflict: the
non-ringing bound caps :math:`m_k` from above, this one floors it from below.

Why it matters for mass scaling
"

Mass scaling -- multiplying the density by :math:`f` to buy a larger increment, the
standard device for a quasi-static explicit run -- lowers :math:`c_d` as
:math:`1/\sqrt{f}`, so the floor above **rises as** :math:`\sqrt{f}`. The two knobs have
to move together. Moving either one alone does nothing, which is easy to verify and
surprising the first time:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - change
     - what it raises
     - what still binds
     - measured :math:`\Delta t`
   * - :math:`f\!:\,1 \to 10^4` at :math:`\eta = 10^{-5}`
     - mechanical, 100-fold
     - non-local
     - 8.852981e-08, *unchanged*
   * - :math:`\eta\!:\,10^{-5} \to 10^{-4}` at :math:`f = 1`
     - non-local, tenfold
     - mechanical
     - 6.063391e-08, *unchanged*
   * - both, :math:`f = 70` and :math:`\eta = 10^{-4}`
     - both
     - mechanical
     - 5.072997e-07, the full :math:`\sqrt{70}`

(:math:`h = 2.5` mm, :math:`l = 5` mm, `GC3D20R`, for which the ratio above gives
:math:`C \approx 40`.)

Two caveats. For :math:`h \sim l` the asymptotic forms do not hold -- the non-local limit
saturates at the reaction term, as noted above -- so on a mesh that coarse the ratio has to
be read off rather than predicted; it comes out *below* its asymptote. And a **parabolic**
non-local field cannot be brought into this regime at all: its limit
:math:`\Delta t \le 2\eta h^2/(C l^2)` is quadratic in :math:`h` and contains no density,
so no mass scaling reaches it and refinement always wins in the end. Mass scaling a
parabolic gradient-damage model therefore lifts the increment straight through a bound that
nothing checks, and the run goes to NaN while still reporting that it finished. Making the
non-local field second order in time is what makes "keep the mechanical problem in charge"
an achievable state.

When it is not achievable, the element says so:
:cpp:func:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeCriticalTimeStepForExplicitDynamics`
warns once per element type when the increment it returns came from the non-local limit
rather than the mechanical one, since the caller receives a single number and cannot
otherwise tell which knob would help it.

What this does not do
^^^^^^^^^^^^^^^^^^^^^

- **It does not remove the artificial lag.** Below :math:`\omega_0` the damped hyperbolic
  equation *is* the parabolic one, with :math:`\eta` in the same place and the same
  meaning. What it removes is the need to raise :math:`\eta` in order to afford the time
  step.
- **It does not damp short-wavelength ringing.** The damping is mass-proportional, so the
  damping ratio :math:`\zeta = \eta/(2\sqrt{m_k(1+C l^2/h^2)})` falls as the frequency
  rises: the reaction mode is critically damped at :math:`m_k = \eta^2/4` and the
  shortest-wavelength modes are barely touched. Where damage is driven by an accumulating
  internal variable, an overshoot of :math:`\bar\varepsilon` above
  :math:`\tilde\varepsilon` is written in irreversibly, so this is worth measuring rather
  than assuming: track :math:`\max(\bar\varepsilon - \tilde\varepsilon)` and the internal
  variable it drives. Reaching the short modes needs a stiffness-proportional term, which
  costs an extra evaluation of the non-local operator and is not implemented.
- **It is a numerical device, not a model.** The physical formulation is the
  :math:`m_k = 0` one. Two meshes of the same material may want different values.

References
^^^^^^^^^^

- Askes, H. & Sluys, L. J. (2002). *Explicit and implicit gradient series in damage
  mechanics*. European Journal of Mechanics A/Solids 21(3), 379--390.
- Askes, H., Bennett, T. & Aifantis, E. C. (2007). *A new formulation and C0
  implementation of dynamically consistent gradient elasticity*. International Journal
  for Numerical Methods in Engineering 72(1), 111--126.
- Peerlings, R. H. J., de Borst, R., Brekelmans, W. A. M. & de Vree, J. H. P. (1996).
  *Gradient enhanced damage for quasi-brittle materials*. International Journal for
  Numerical Methods in Engineering 39(19), 3391--3403.

Usage
-----

The micro-inertia is an **element** property for the same reason the artificial bulk
viscosity is one: it is a numerical device, the physical model is the one without it, and
the value that pays off depends on the mesh. It is assigned through the named-property
interface under the name ``nonlocal micro inertia``, which takes one value per non-local
variable:

.. code-block:: cpp

   const double microInertia = 2.5e-11; // seconds squared, = eta^2 / 4 for eta = 1e-5 s
   element->assignProperty( "nonlocal micro inertia", &microInertia, 1 );

From EdelweissFE the property is reached through the ``*elementproperty`` keyword, and the
field has to be moved from the first-order to the second-order integration scheme and
declared as carrying a micro-inertia rather than a mass:

.. code-block:: none

   *elementproperty, elSet=concrete, propertyName=nonlocal micro inertia
   2.5e-11

   *solver, solver=NEDParallel, name=theSolver
   second-order-fields="displacement, nonlocal damage"
   micro-inertia-fields="nonlocal damage"

Both directions of that declaration are checked: a micro-inertia the elements carry but
the solver was not told about, and one the solver expects but the elements do not
assemble, are both refused with a message rather than integrated.

**Unset, the feature is completely inert.** An element that was never given the property
reports a zero micro-inertia, its non-local field stays first order in time, and the
stable increment is the mechanical one exactly as before.

Notes and limitations
---------------------

- **The lumping matches the mass.** The non-local block is weighted exactly as
  ``computeLumpedInertia`` weights it, because the stable increment is read off that
  distribution and the lightest node sets the highest frequency.
- **The non-local interaction parameter is a material response.** Bounding the eigenvalue
  needs :math:`c`, which the material interface only exposes through a stress evaluation.
  It is queried once per element per step on a scratch copy of the state variables and
  with a zero strain increment, so it leaves no trace; the stable increment is computed
  when the solver builds its system, not per increment.
- **Explicit only.** An implicit solve has no stability limit to relieve and needs the
  elliptic form.

Implementation
--------------

.. doxygenfunction:: Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeLumpedNonlocalMicroInertia

.. doxygenfunction:: Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeCriticalTimeStepForExplicitDynamics
