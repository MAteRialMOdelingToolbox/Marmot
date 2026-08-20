.. _displacementpressurefinitestrainelement:

Displacement-pressure mixed finite-strain element
==================================================

Theory
------

Incompressible hyperelastic materials (e.g. :ref:`ogden`) carry no volumetric energy term, so a
plain displacement element would leave the volumetric response entirely unconstrained. Following
the classical mixed formulation of Simo, Taylor and Pister, incompressibility is instead enforced
weakly by introducing an independent pressure field :math:`p`, conjugate to the constraint
:math:`J-1=0` (:math:`J=\det\boldsymbol F`), via the two-field functional

.. math::

   \Pi(\boldsymbol u, p)
   =
   \int_{\Omega_0} \Psi_{\rm iso}(\boldsymbol F) \, \mathrm dV_0
   +
   \int_{\Omega_0} p\,(J-1) \, \mathrm dV_0 .

Taking variations with respect to :math:`\boldsymbol u` and :math:`p` gives the total Kirchhoff
stress

.. math::

   \boldsymbol\tau = \boldsymbol\tau_{\rm iso}(\boldsymbol F) + p\,J\,\boldsymbol I,

and the two coupled weak-form equations

.. math::

   \boldsymbol R_u = \int_{\Omega_0} (\nabla_x \boldsymbol N_u)^{\mathsf T}\,\boldsymbol\tau
   \, \mathrm dV_0 = \boldsymbol 0,
   \qquad
   R_p = \int_{\Omega_0} \boldsymbol N_p^{\mathsf T}\,(J-1)\, \mathrm dV_0 = 0,

where :math:`\boldsymbol\tau_{\rm iso}` and its algorithmic tangent are supplied by any attached
purely isochoric ``MarmotMaterialFiniteStrain`` instance.

Interpolation and inf-sup stability
------------------------------------

Because the energy is linear in :math:`p`, the discrete stiffness matrix has a vanishing pressure
block :math:`K_{pp}=0`: the system is a genuine saddle-point (mixed) problem, and the
displacement/pressure interpolation pair must satisfy the discrete inf-sup (Ladyzhenskaya-Babuska-
Brezzi) condition, or the discrete pressure field develops spurious checkerboard oscillations.

This element uses a **quadratic displacement field on all nodes** combined with a **continuous
linear pressure field on the corner-node subset** - a P2/P1-type pairing analogous to the classical
Taylor-Hood element, and inf-sup stable. Concretely:

.. list-table::
   :header-rows: 1
   :align: left

   * - Element
     - Displacement nodes
     - Pressure nodes
   * - Tetra10 (``C3D10MP``)
     - 10 (quadratic)
     - 4 (corner, linear)
   * - Hexa20 (``C3D20MP``)
     - 20 (quadratic serendipity)
     - 8 (corner, linear)

Both fields share the same full Gauss quadrature rule - no selective/reduced integration is needed,
since inf-sup stability follows from the interpolation pair itself rather than from
under-integrating the pressure.

.. note::

   Restricted to 3D solid elements. A 2D plane-strain variant would additionally require the
   finite-strain plane-strain 3D-embedding machinery used by
   :ref:`displacementfinitestrainelement` and is not implemented here.

Implementation
--------------

Per quadrature point, the deformation gradient :math:`\boldsymbol F` is computed from the total
displacement dofs, the attached material supplies :math:`\boldsymbol\tau_{\rm iso}` and
:math:`\partial\boldsymbol\tau_{\rm iso}/\partial\boldsymbol F`, and the pressure :math:`p` is
interpolated from the pressure dofs. The tangent blocks follow the standard nonlinear (updated
Lagrangian style) split into a material part and a geometric part for :math:`K_{uu}`, using the
*total* stress and tangent :math:`\boldsymbol\tau=\boldsymbol\tau_{\rm iso}+pJ\boldsymbol I` and
:math:`\partial\boldsymbol\tau/\partial\boldsymbol F = \partial\boldsymbol\tau_{\rm
iso}/\partial\boldsymbol F + p\,\partial J/\partial\boldsymbol F \otimes \boldsymbol I` in exactly
the same assembly used for a single-field finite-strain element (see
:ref:`displacementfinitestrainelement`) - no second derivative of :math:`J` is required. The
coupling blocks follow directly as

.. math::

   K_{up} = K_{pu}^{\mathsf T} = \int_{\Omega_0} J\,(\nabla_x \boldsymbol N_u) \otimes
   \boldsymbol N_p \, \mathrm dV_0 .

.. doxygenclass:: Marmot::Elements::DisplacementPressureFiniteStrainElement
   :allow-dot-graphs:
