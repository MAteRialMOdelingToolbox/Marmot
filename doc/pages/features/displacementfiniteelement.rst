Displacement Finite Element
===========================

Theory
------

This element implements a small-strain displacement-based Finite Element formulation.
The incremental strain is computed from the displacement increment via the kinematic
matrix B, and internal forces and consistent tangents are evaluated by Gauss quadrature.

.. math::

   \Delta \eps = \mathbf{B}\, \Delta \mathbf{u}, \qquad
   \mathbf{K}_e = \sum_{qp} \mathbf{B}^\mathsf{T}\, \mathbf{C}\, \mathbf{B}\, J_0 w, \qquad
   \mathbf{P}_e = \sum_{qp} \mathbf{B}^\mathsf{T}\, \sig\, J_0 w.

where :math:`\Delta \eps` is the linearized strain tensor in Voigt notation,
:math:`\mathbf{B}` is the strain–displacement matrix, :math:`\Delta \mathbf{u}` the nodal displacement increment;
:math:`\mathbf{K}_e` the element tangent stiffness, :math:`\mathbf{C}` the consistent material tangent;
:math:`J_0 = \det \mathbf{J}` the Jacobian determinant and :math:`w` the quadrature weight;
:math:`\mathbf{P}_e` the element internal force vector, :math:`\sig` the Cauchy stress tensor in Voigt notation,
and :math:`\sum_{qp}` denotes summation over quadrature points.

The mass matrix and body-force vector are obtained from the expressions:

.. math::

   \mathbf{M}_e = \sum_{qp} \rho\, \mathbf{N}^\mathsf{T} \mathbf{N}\, J_0 w, \qquad
   \mathbf{P}_e^{(b)} = \sum_{qp} \mathbf{N}^\mathsf{T} \mathbf{f}\, J_0 w.

where :math:`\mathbf{M}_e` is the consistent mass matrix, :math:`\rho` the mass density,
:math:`\mathbf{N}` the shape-function matrix of the displacement field, and
:math:`\mathbf{f}` the body-force vector per unit volume.

The consistent mass is integrated with the **full** Gauss rule of the element's shape (3x3x3 for a
20-node hexahedron, 3x3 for an 8-node quadrilateral), also for a reduced-integration element: its
own rule cannot produce a mass matrix of full rank -- 8 points give a 20-node hexahedron a mass of
rank 24 out of 60 -- and a singular mass leaves the modes of its null space without inertia and
makes an equilibrium solve for the initial acceleration, :math:`\mathbf{M}\,\mathbf{a}_0 = \mathbf{R}`,
ill-posed. The density at a point of that rule is taken from the nearest quadrature point of the
element, which is exact for a density that is constant over the element. The same applies to the
finite-strain (UL) and the gradient-enhanced displacement elements, for the latter to the
micro-inertia of the non-local block as well.

Surface tractions and pressures on a boundary face :math:`\Gamma_e` are integrated as

.. math::

   \mathbf{P}_e^{(t)} = \int_{\Gamma_e} \mathbf{N}^\mathsf{T} \, \mathbf{t} \, \mathrm{d}\Gamma, \qquad
   \mathbf{P}_e^{(p)} = - \int_{\Gamma_e} p \, \mathbf{N}^\mathsf{T} \, \mathbf{n} \, \mathrm{d}\Gamma.

where :math:`\Gamma_e` is the element boundary, :math:`\mathbf{n}` its outward unit normal,
:math:`\mathbf{t}` the surface traction vector, :math:`p` the pressure magnitude, and
:math:`\mathbf{N}` the shape-function matrix. The minus sign for pressure follows
an outward-normal convention.

For 2D and 1D problems, the integration measure :math:`J_0 w` is scaled by thickness
:math:`t` or cross-sectional area :math:`A`, respectively: :math:`J_0 w \leftarrow J_0 w\, t`
(2D), :math:`J_0 w \leftarrow J_0 w\, A` (1D).
Here, :math:`t` denotes thickness and :math:`A` denotes cross-sectional area.

Constitutive updates are carried out per quadrature point using a Marmot material model:

- 3D solid

- Plane stress

- Plane strain

- 1D uniaxial stress

Each quadrature point stores stress, strain and a material state vector; the element
accumulates the history incrementally and may suggest a reduced time step if required
by the material.


Implementation
--------------

.. doxygenclass:: Marmot::Elements::DisplacementFiniteElement
   :allow-dot-graphs:
