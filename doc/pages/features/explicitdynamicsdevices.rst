Explicit dynamics devices
=========================

Two numerical devices for explicit dynamics, developed and calibrated together: artificial bulk
viscosity, which damps volumetric ringing in the mechanical field, and non-local micro-inertia,
which makes a gradient-enhanced field second order in time so its stable increment scales linearly
rather than quadratically with the element size. Neither is a material model; both are element
properties, inert unless assigned.

Artificial bulk viscosity
-------------------------

Theory
^^^^^^

Explicit dynamics carries no numerical dissipation of its own. The central-difference integrator
is non-dissipative by construction, so any energy released inside the mesh -- most sharply when an
element loses its stiffness on cracking -- is converted to kinetic energy and rings at the highest
frequency the mesh can carry until something removes it. Artificial bulk viscosity is the standard
device for removing it.

A viscous stress is added to the volumetric part of the stress at every quadrature point:

.. math::

   \sig_\mathrm{bv} = \left[ b_1\, \rho\, c_d\, L_e\, \dot{\eps}_\mathrm{vol}
                    - \rho\, (b_2 L_e)^2\, \dot{\eps}_\mathrm{vol}^2\,
                      H(-\dot{\eps}_\mathrm{vol}) \right] \mathbf{I}

where :math:`\dot{\eps}_\mathrm{vol} = \mathrm{tr}(\dot{\boldsymbol{\eps}})` is the volumetric
strain rate (negative in compression), :math:`\rho` the mass density, :math:`c_d` the dilatational
wave speed, :math:`L_e` the element's smallest physical extent, :math:`H` the Heaviside function,
and :math:`b_1`, :math:`b_2` two dimensionless coefficients.

The **linear** (Landshoff) term damps element-level ringing and acts on both signs of the
volumetric rate. The **quadratic** (Von Neumann--Richtmyer) term spreads a shock over a few
elements and is restricted to compression, so that it cannot resist the opening of a crack or a
void. Both terms are dissipative by construction -- the power density
:math:`\sig_\mathrm{bv} : \dot{\boldsymbol{\eps}} \ge 0` for either sign of the rate -- so the term
can only remove energy, and negative coefficients are rejected at assignment for the same reason.

This is a **numerical device, not a material model**. It adds no state, the constitutive law never
sees it, and it does not enter the stored or reported stress -- it is added only to the stress that
is integrated into the internal force vector.

References
""""""""""

- VonNeumann, J. & Richtmyer, R. D. (1950). *A method for the numerical calculation of hydrodynamic
  shocks*. Journal of Applied Physics 21(3), 232--237. https://doi.org/10.1063/1.1699639
- Landshoff, R. (1955). *A numerical method for treating fluid flow in the presence of shocks*. Los
  Alamos Scientific Laboratory report LA-1930.

Usage
^^^^^

Bulk viscosity is an **element** property rather than a material one: the same material integrated
implicitly needs none of it, and two meshes of the same material may want different amounts. It is
assigned through the named-property interface under the name ``bulk viscosity``, which takes
exactly two values, :math:`b_1` and :math:`b_2`, and defaults to zero -- inactive -- unless
assigned:

.. code-block:: cpp

   const std::vector< double > coefficients = { 0.06, 1.2 };
   element->assignProperty( "bulk viscosity", coefficients.data(), coefficients.size() );

From EdelweissFE the same property is reached through the ``*elementproperty`` keyword:

.. code-block:: none

   *elementproperty, elSet=concrete, propertyName=bulk viscosity
   0.06, 1.2

``0.06`` and ``1.2`` are the values Abaqus/Explicit applies, a common choice rather than the
default: :math:`b_1` is sized to damp the highest resolvable frequency of the mesh, not to
represent any physical dissipation, and a problem dominated by a single sharp energy release may
need considerably more.

Degradation with damage
^^^^^^^^^^^^^^^^^^^^^^^^

The linear term acts on both signs of the volumetric rate, so in an element whose material has
already failed it is a viscous resistance to the crack **opening** -- one that never relaxes, and
whose work is charged to the dissipated energy, inflating any fracture energy measured from the
load-displacement response.

Measured on a calibrated gradient-enhanced damage-plasticity bar in tension, at the coefficients
above, the dissipated work per unit fracture area rose by 27 to 55 percent depending on the mesh.
Applying bulk viscosity everywhere **except** the elements that damaged changed the same quantity
by :math:`-1.6` percent: essentially all of the error is generated inside the damaged elements, and
slowing the loading does not remove it (a ten-fold slower ramp only brought 38 percent down to 10
percent, because the strain rate inside a localising band is set by the softening, not by the
imposed rate).

The optional degradation addresses this by scaling the viscous stress with the material's
remaining stiffness,

.. math::

   \sigma_\mathrm{bv} \leftarrow \left( \frac{c}{c_0} \right)^{n} \sigma_\mathrm{bv},

with :math:`c` the wave speed at the current state, :math:`c_0` the cached undamaged reference,
and the ratio clamped to :math:`[0,1]`. This is the same argument that restricts the quadratic term
to compression, applied to the **state** rather than to the sign of the rate. Because a wave speed
is the square root of a stiffness, the exponent selects what the stress follows for a model whose
tangent degrades as :math:`(1-\omega)\,\mathbb{C}_0`:

.. list-table::
   :header-rows: 1
   :widths: 10 40

   * - :math:`n`
     - Effect
   * - ``0``
     - Off, and the default. The current wave speed is never evaluated.
   * - ``1``
     - Scales with the wave speed, i.e. with :math:`\sqrt{1-\omega}`.
   * - ``2``
     - Scales with the tangent stiffness, i.e. with :math:`1-\omega`.

It is requested through a second named property, taking exactly one value:

.. code-block:: none

   *elementproperty, elSet=concrete, propertyName=bulk viscosity damage degradation
   2.0

.. warning::

   The stress is degraded with the current **tangent**, not with a damage variable: no material
   interface here reports damage, and adding one would change the vtable of every element. For a
   quasi-brittle material in tension the two coincide, because the softening *is* the damage. For a
   model that merely yields, the algorithmic tangent also drops and the viscous stress is degraded
   by plastic flow rather than by cracking, which is not what the device is for.

.. note::

   This is why the degradation is opt-in rather than the default: it needs the material's current
   tangent, which costs a full constitutive evaluation per quadrature point per increment -- exactly
   the cost the cached reference wave speed exists to avoid.

Notes and limitations
^^^^^^^^^^^^^^^^^^^^^^

- **Volumetric only.** Ringing that is predominantly deviatoric or flexural is not reached, which
  is a property of the device and not of this implementation.
- **The wave speed is cached.** Querying a material for its current wave speed costs a full
  constitutive evaluation, unaffordable per quadrature point per explicit increment, so the
  **undamaged** wave speed is cached on first use. That leaves the damping slightly stronger than a
  current-stiffness value would as the material softens -- the safe direction for a device whose
  purpose is to remove energy -- and it is the reference the optional degradation above measures
  against.
- **The characteristic length is the element's smallest physical extent**, twice the smallest
  singular value of the Jacobian -- the same length the stable time increment is computed from,
  shared through one function so the two cannot drift apart.
- **Plane stress is approximate.** The out-of-plane strain is not carried by the element's
  kinematics there, so the trace is taken over the in-plane components only.
- **Explicit only.** Evaluated in ``computeKernelsExplicit``, not part of the implicit kernels,
  where it has no purpose and would need a consistent tangent contribution.

Supported elements
^^^^^^^^^^^^^^^^^^

- :cpp:class:`Marmot::Elements::DisplacementFiniteElement`
- :cpp:class:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement`

Not yet supported by ``DisplacementFiniteStrainULElement``: a finite-strain formulation needs the
volumetric rate taken as :math:`\mathrm{tr}(\mathbf{D}) = \dot{J}/J` and the characteristic length
measured in the current configuration, a separate derivation rather than a reuse of the small-strain
one.

Non-local micro-inertia (hyperbolic gradient damage)
-----------------------------------------------------

Theory
^^^^^^

The implicit gradient-enhanced formulation regularises softening by solving a Helmholtz equation
for the non-local variable :math:`\bar\varepsilon` alongside the momentum balance,

.. math::

   \bar\varepsilon - c\, \nabla^2 \bar\varepsilon = \tilde\varepsilon ,
   \qquad c = l^2 ,

with :math:`\tilde\varepsilon` the local driving variable and :math:`l` the internal length. An
explicit solver cannot march an equation with no time derivative, so it is made parabolic by adding
a viscosity :math:`\eta`,

.. math::

   \eta\, \dot{\bar\varepsilon} + \bar\varepsilon - c\, \nabla^2 \bar\varepsilon
     = \tilde\varepsilon ,

and integrated with forward Euler. That works, but the largest eigenvalue of the discrete operator
is :math:`\lambda_\mathrm{max} \approx \tfrac{1}{\eta}(1 + C l^2/h^2)`, so the stable increment

.. math::

   \Delta t \le \frac{2\eta}{1 + C\,l^2/h^2}
             \;\xrightarrow{\;h \ll l\;}\;
             \frac{2\eta\,h^2}{C\,l^2} \;\propto\; h^2

falls off with the **square** of the element size (:math:`C` a discretisation constant, order 4 in
1D and roughly 12 to 24 for a hexahedron). Under refinement this is the binding constraint long
before the mechanical one is, and the only knob it offers is :math:`\eta` itself -- the artificial
delay between the damage front and the strain concentration driving it. Buying stability with
:math:`\eta` buys it by making the regularisation lag more.

Adding a **micro-inertia** :math:`m_k` breaks that coupling:

.. math::

   m_k\, \ddot{\bar\varepsilon} + \eta\, \dot{\bar\varepsilon} + \bar\varepsilon
     - c\, \nabla^2 \bar\varepsilon = \tilde\varepsilon .

This is now a damped wave (telegraph, or damped Klein--Gordon) equation for a damage front
travelling at :math:`c_k = \sqrt{c/m_k} = l/\sqrt{m_k}`, with :math:`m_k` in seconds squared. The
viscosity has not changed meaning; it has changed *role*. What was the coefficient of the highest
time derivative is now the damping.

Stability
"""""""""

With a lumped micro-inertia and central differences the limit is

.. math::

   \omega_\mathrm{max} = \sqrt{\frac{1 + C\,l^2/h^2}{m_k}} ,
   \qquad
   \zeta = \frac{\eta/m_k}{2\,\omega_\mathrm{max}} ,
   \qquad
   \Delta t \le \frac{2}{\omega_\mathrm{max}}
                 \left( \sqrt{1+\zeta^2} - \zeta \right) .

For :math:`h \ll l` this is :math:`\Delta t \approx 2\sqrt{m_k}\,h / (\sqrt{C}\, l)`: **linear in
the element size**, the same scaling elastodynamics has. On a mesh coarser than the internal length
the reaction term takes over and the limit saturates at :math:`2\sqrt{m_k}`, so the field can never
demand a smaller increment than that regardless of how coarse the mesh is.

Two things are easy to get wrong here, each costing a factor of two or more: the continuum estimate
:math:`\Delta t \le h/c_k` **overestimates** the limit, because the discrete Laplacian's eigenvalue
constant contributes :math:`2/\sqrt{C} \approx 0.4` to :math:`0.6`, not 1; and damping **lowers**
the central-difference limit, never raises it. Both are accounted for in
:cpp:func:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeCriticalTimeStepForExplicitDynamics`,
which bounds :math:`\lambda_\mathrm{max}` with a Gershgorin estimate over the assembled operator
rather than a wave speed, and takes the minimum against the mechanical limit.

Choosing the two parameters
"""""""""""""""""""""""""""

They are not independent. The zeroth-order reaction term is an oscillator of frequency
:math:`\omega_0 = 1/\sqrt{m_k}`, damped by :math:`\eta`; requiring that it not ring gives
:math:`m_k \le \eta^2/4`, and the largest admissible micro-inertia is the best one, because
:math:`\Delta t` grows with :math:`\sqrt{m_k}`. So **there is one free parameter, not two**: pick
:math:`\eta` for the lag you accept -- exactly as before -- and take :math:`m_k = \eta^2/4`. The
resulting damage wave speed :math:`c_k = 2l/\eta` should then be checked against the loading rate,
not chosen from it.

At that choice the gain over the parabolic scheme *at the same artificial lag* is

.. math::

   \frac{\Delta t_\mathrm{hyperbolic}}{\Delta t_\mathrm{parabolic}}
     = \frac{\sqrt{1 + C\,l^2/h^2}}{2}
     \;\approx\; \frac{\sqrt{C}\,l}{2h} ,

a factor that **grows with every refinement level** -- the finer the mesh, the more the
second-order form is worth.

Keeping the mechanical problem in charge
""""""""""""""""""""""""""""""""""""""""

There is a second criterion on :math:`\eta`, deciding whether the usual quasi-static-explicit
tricks work at all. The two limits the element takes the minimum of,

.. math::

   \Delta t_\mathrm{mech} \approx k_c\,\frac{h}{c_d} ,
   \qquad
   \Delta t_\mathrm{nl} \approx \frac{2 k_c}{\sqrt{C}}\,\frac{\sqrt{m_k}\,h}{l}
   \qquad (h \ll l),

are **both linear in** :math:`h`, so their ratio contains no mesh size at all:

.. math::

   \frac{\Delta t_\mathrm{nl}}{\Delta t_\mathrm{mech}} = \frac{\eta\,c_d}{l\,\sqrt{C}} ,

the courant number cancelling because it applies to both. Which of the two is in charge is thus a
property of the *parameters*, not of the discretisation: fix it once and it holds at every
refinement level. Requiring the mechanical limit to bind gives

.. math::

   \eta \;\ge\; \frac{\sqrt{C}\,l}{c_d}
   \qquad\Longleftrightarrow\qquad
   m_k \;\ge\; \frac{C}{4}\left(\frac{l}{c_d}\right)^{\!2} ,

i.e. **the micro-inertia must exceed the square of the time a mechanical wave needs to cross one
non-local length**. Together with :math:`m_k \le \eta^2/4` this is a window, not a conflict: the
non-ringing bound caps :math:`m_k` from above, this one floors it from below.

Why it matters for mass scaling
"""""""""""""""""""""""""""""""

Mass scaling -- multiplying the density by :math:`f` to buy a larger increment -- lowers
:math:`c_d` as :math:`1/\sqrt{f}`, so the floor above **rises as** :math:`\sqrt{f}`. The two knobs
have to move together; moving either alone does nothing:

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

(:math:`h = 2.5` mm, :math:`l = 5` mm, ``GC3D20R``, for which the ratio above gives
:math:`C \approx 40`.)

Two caveats. For :math:`h \sim l` the asymptotic forms do not hold -- the non-local limit saturates
at the reaction term -- so the ratio has to be read off rather than predicted there. And a
**parabolic** non-local field can never be brought into this regime: its limit is quadratic in
:math:`h` and density-free, so mass scaling lifts the increment straight through a bound nothing
checks, and the run goes to NaN while still reporting that it finished. Making the field second
order in time is what makes "keep the mechanical problem in charge" achievable.

When it is not achievable, the element says so:
:cpp:func:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeCriticalTimeStepForExplicitDynamics`
warns once per element type when the increment it returns came from the non-local limit rather
than the mechanical one, since the caller receives a single number and cannot otherwise tell which
knob would help.

What this does not do
^^^^^^^^^^^^^^^^^^^^^^

- **It does not remove the artificial lag.** Below :math:`\omega_0` the damped hyperbolic equation
  *is* the parabolic one, with :math:`\eta` in the same place and the same meaning. What it removes
  is the need to raise :math:`\eta` to afford the time step.
- **It does not damp short-wavelength ringing.** The damping is mass-proportional, so
  :math:`\zeta = \eta/(2\sqrt{m_k(1+C l^2/h^2)})` falls as the frequency rises: the reaction mode is
  critically damped at :math:`m_k = \eta^2/4` and the shortest-wavelength modes are barely touched.
  Where damage rides an accumulating internal variable, an overshoot of :math:`\bar\varepsilon`
  above :math:`\tilde\varepsilon` is written in irreversibly -- worth measuring
  (:math:`\max(\bar\varepsilon - \tilde\varepsilon)` and the internal variable it drives) rather
  than assuming. Reaching the short modes needs a stiffness-proportional term, not implemented.
- **It is a numerical device, not a model.** The physical formulation is the :math:`m_k = 0` one.

References
""""""""""

- Askes, H. & Sluys, L. J. (2002). *Explicit and implicit gradient series in damage mechanics*.
  European Journal of Mechanics A/Solids 21(3), 379--390.
- Askes, H., Bennett, T. & Aifantis, E. C. (2007). *A new formulation and C0 implementation of
  dynamically consistent gradient elasticity*. International Journal for Numerical Methods in
  Engineering 72(1), 111--126.
- Peerlings, R. H. J., de Borst, R., Brekelmans, W. A. M. & de Vree, J. H. P. (1996). *Gradient
  enhanced damage for quasi-brittle materials*. International Journal for Numerical Methods in
  Engineering 39(19), 3391--3403.

Usage
^^^^^

The micro-inertia is an **element** property for the same reason artificial bulk viscosity is one:
a numerical device, inert unless assigned, whose useful value depends on the mesh. It is assigned
under the name ``nonlocal micro inertia``, one value per non-local variable:

.. code-block:: cpp

   const double microInertia = 2.5e-11; // seconds squared, = eta^2 / 4 for eta = 1e-5 s
   element->assignProperty( "nonlocal micro inertia", &microInertia, 1 );

From EdelweissFE the field also has to be moved from the first-order to the second-order
integration scheme and declared as carrying a non-mechanical inertia rather than a mass:

.. code-block:: none

   *elementproperty, elSet=concrete, propertyName=nonlocal micro inertia
   2.5e-11

   *solver, solver=NEDParallel, name=theSolver
   second-order-fields="displacement, nonlocal damage"
   non-mechanical-inertia-fields="nonlocal damage"

Both directions of that declaration are checked: a micro-inertia the elements carry but the solver
was not told about, and one the solver expects but the elements do not assemble, are both refused
with a message rather than integrated. **Unset, the feature is completely inert**: an element that
was never given the property reports zero, its non-local field stays first order in time, and the
stable increment is the mechanical one exactly as before.

Notes and limitations
^^^^^^^^^^^^^^^^^^^^^^

- **The lumping matches the mass.** Both are assembled by
  :cpp:func:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeLumpedInertia`,
  using the same weights on the non-local block as on the displacement one, because the stable
  increment is read off that distribution and the lightest node sets the highest frequency.
- **The non-local interaction parameter is a material response.** Bounding the eigenvalue needs
  :math:`c`, exposed only through a stress evaluation -- queried once per element per step on a
  scratch copy of the state variables and a zero strain increment, so it leaves no trace.
- **Explicit only.** An implicit solve has no stability limit to relieve and needs the elliptic
  form.

Implementation
--------------

.. doxygennamespace:: Marmot::FiniteElement::BulkViscosity
   :members:

.. doxygenfunction:: Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeLumpedInertia

.. doxygenfunction:: Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeLumpedDamping

.. doxygenfunction:: Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement::computeCriticalTimeStepForExplicitDynamics
