Artificial bulk viscosity
=========================

Theory
------

Explicit dynamics carries no numerical dissipation of its own. The central-difference
integrator is non-dissipative by construction, so any energy released inside the mesh --
most sharply when an element loses its stiffness on cracking -- is converted to kinetic
energy and rings at the highest frequency the mesh can carry until something removes it.
Artificial bulk viscosity is the standard device for removing it.

A viscous stress is added to the volumetric part of the stress at every quadrature point:

.. math::

   \sig_\mathrm{bv} = \left[ b_1\, \rho\, c_d\, L_e\, \dot{\eps}_\mathrm{vol}
                    - \rho\, (b_2 L_e)^2\, \dot{\eps}_\mathrm{vol}^2\,
                      H(-\dot{\eps}_\mathrm{vol}) \right] \mathbf{I}

where :math:`\dot{\eps}_\mathrm{vol} = \mathrm{tr}(\dot{\boldsymbol{\eps}})` is the
volumetric strain rate (negative in compression), :math:`\rho` the mass density,
:math:`c_d` the dilatational wave speed, :math:`L_e` the element's smallest physical
extent, :math:`H` the Heaviside function, and :math:`b_1`, :math:`b_2` two dimensionless
coefficients.

The **linear** (Landshoff) term damps element-level ringing and acts on both signs of the
volumetric rate. The **quadratic** (Von Neumann--Richtmyer) term spreads a shock over a
few elements and is restricted to compression, so that it cannot resist the opening of a
crack or a void.

Both terms are dissipative by construction. The power density is

.. math::

   \sig_\mathrm{bv} : \dot{\boldsymbol{\eps}}
     = b_1 \rho c_d L_e \dot{\eps}_\mathrm{vol}^2
     - \rho (b_2 L_e)^2 \dot{\eps}_\mathrm{vol}^3 H(-\dot{\eps}_\mathrm{vol})
     \;\ge\; 0,

with both contributions non-negative, so the term can only remove energy from the system
and never add it. Negative coefficients are rejected at assignment for the same reason.

This is a **numerical device, not a material model**. It adds no state, the constitutive
law never sees it, and it does not enter the stored or reported stress -- it is added only
to the stress that is integrated into the internal force vector.

References
^^^^^^^^^^

- VonNeumann, J. & Richtmyer, R. D. (1950). *A method for the numerical calculation of
  hydrodynamic shocks*. Journal of Applied Physics 21(3), 232--237.
  https://doi.org/10.1063/1.1699639
- Landshoff, R. (1955). *A numerical method for treating fluid flow in the presence of
  shocks*. Los Alamos Scientific Laboratory report LA-1930.

Usage
-----

Bulk viscosity is an **element** property rather than a material one: the same material
integrated implicitly needs none of it, and two meshes of the same material may want
different amounts. It is assigned through the named-property interface under the name
``bulk viscosity``, which takes exactly two values, :math:`b_1` and :math:`b_2`:

.. code-block:: cpp

   const std::vector< double > coefficients = { 0.06, 1.2 };
   element->assignProperty( "bulk viscosity", coefficients.data(), coefficients.size() );

From EdelweissFE the same property is reached through the ``*elementproperty`` keyword:

.. code-block:: none

   *elementproperty, elSet=concrete, propertyName=bulk viscosity
   0.06, 1.2

The values ``0.06`` and ``1.2`` are the defaults Abaqus/Explicit applies. They are
defaults rather than recommendations: :math:`b_1` is sized to damp the highest resolvable
frequency of the mesh, not to represent any physical dissipation, and a problem dominated
by a single sharp energy release may need considerably more.

**Unset, the feature is completely inert.** Nothing is evaluated when both coefficients
are zero, so a model that does not ask for bulk viscosity produces bit-identical results
to one run before the feature existed.

Degradation with damage
-----------------------

The linear term acts on both signs of the volumetric rate, which means that in an element
whose material has already failed it is a viscous resistance to the crack **opening** --
and one that never relaxes. Its work is charged to the dissipated energy, so it inflates
any fracture energy measured from the load-displacement response.

That effect is not small and it is not spread out. Measured on a calibrated
gradient-enhanced damage-plasticity bar in tension, at the coefficients above, the
dissipated work per unit fracture area rose by 27 to 55 percent depending on the mesh, and
falling only as about :math:`h^{0.4}` it was still a quarter of the total at the finest
mesh. Applying bulk viscosity everywhere **except** the elements that damaged changed the
same quantity by :math:`-1.6` percent: essentially the whole of the error is generated
inside the damaged elements, and slowing the loading does not remove it (a ten-fold slower
ramp only brought 38 percent down to 10 percent, because the strain rate inside a
localising band is set by the softening and not by the imposed rate).

The optional degradation addresses this by scaling the viscous stress with the material's
remaining stiffness:

.. math::

   \sigma_\mathrm{bv} \leftarrow \left( \frac{c}{c_0} \right)^{n} \sigma_\mathrm{bv}

with :math:`c` the wave speed at the current state, :math:`c_0` the cached undamaged
reference, and the ratio clamped to :math:`[0,1]` so that a stiffer-than-virgin tangent
cannot amplify the damping. This is the same argument that restricts the quadratic term to
compression, applied to the **state** rather than to the sign of the rate.

Because a wave speed is the square root of a stiffness, the exponent selects what the
stress follows for a model whose tangent degrades as :math:`(1-\omega)\,\mathbb{C}_0`:

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

It is requested through a second named property, which takes exactly one value:

.. code-block:: none

   *elementproperty, elSet=concrete, propertyName=bulk viscosity
   0.06, 1.2

   *elementproperty, elSet=concrete, propertyName=bulk viscosity damage degradation
   2.0

.. warning::

   The stress is degraded with the current **tangent**, not with a damage variable: no
   material interface here reports damage, and adding one would change the vtable of every
   element. For a quasi-brittle material in tension the two coincide, because the softening
   *is* the damage. For a model that merely yields, the algorithmic tangent also drops and
   the viscous stress is then degraded by plastic flow rather than by cracking, which is
   not what the device is for.

.. note::

   This is why the degradation is opt-in rather than the default: it needs the material's
   current tangent, which costs a full constitutive evaluation per quadrature point per
   increment -- exactly the cost the cached reference wave speed exists to avoid.

Notes and limitations
---------------------

- **Volumetric only.** The term damps the volumetric mode. Ringing that is predominantly
  deviatoric or flexural is not reached by it, which is a property of the device and not
  of this implementation.
- **The wave speed is cached.** Querying a material for its current wave speed costs a
  full constitutive evaluation -- for a gradient-enhanced damage-plasticity model, a
  complete return mapping -- which is not affordable per quadrature point per explicit
  increment. The wave speed of the **undamaged** material is therefore cached on first
  use. The term exists to damp the highest frequency the *mesh* can carry, which the
  undamaged material sets, and holding it fixed as the material softens leaves the damping
  slightly stronger than a current-stiffness value would: the safe direction for a device
  whose purpose is to remove energy. Where that is the wrong trade -- above all when the
  quantity of interest is a fracture energy -- the cached value becomes the reference of
  the opt-in degradation described above, which does pay for the current wave speed.
- **The characteristic length is the element's smallest physical extent**, twice the
  smallest singular value of the Jacobian -- the same length the stable time increment is
  computed from, shared through one function so the two cannot drift apart.
- **Plane stress is approximate.** The out-of-plane strain is not carried by the element's
  kinematics there, so the trace is taken over the in-plane components only. In 3D and in
  plane strain the term uses the exact volumetric strain increment.
- **Explicit only.** The viscous stress is evaluated in ``computeKernelsExplicit`` and is
  not part of the implicit kernels, where it has no purpose and would need a consistent
  tangent contribution.

Supported elements
------------------

- :cpp:class:`Marmot::Elements::DisplacementFiniteElement`
- :cpp:class:`Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement`

Not yet supported by ``DisplacementFiniteStrainULElement``: a finite-strain formulation
needs the volumetric rate taken as :math:`\mathrm{tr}(\mathbf{D}) = \dot{J}/J` and the
characteristic length measured in the current configuration, which is a separate
derivation rather than a reuse of the small-strain one.

Implementation
--------------

.. doxygennamespace:: Marmot::FiniteElement::BulkViscosity
   :members:
