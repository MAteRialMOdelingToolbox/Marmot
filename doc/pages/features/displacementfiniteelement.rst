Displacement Finite Element
===========================

Theory
------

This element employs a small-strain, displacement-based formulation.  The incremental
strain is obtained from the displacement increment using the kinematic
(strain–displacement) matrix :math:`B`.  Internal forces and consistent tangents are
evaluated by Gauss quadrature.

.. math::

   \Delta \boldsymbol{\varepsilon} = \mathbf{B}\, \Delta \mathbf{u}, \qquad
   \mathbf{K}_e = \sum_{qp} \mathbf{B}^\mathsf{T}\, \mathbf{C}\, \mathbf{B}\, J_0 w, \qquad
   \mathbf{P}_e = \sum_{qp} \mathbf{B}^\mathsf{T}\, \boldsymbol{\sigma}\, J_0 w.

The mass matrix and body-force vector are obtained from the expressions:

.. math::

   \mathbf{M}_e = \sum_{qp} \rho\, \mathbf{N}^\mathsf{T} \mathbf{N}\, J_0 w, \qquad
   \mathbf{P}_e^{(b)} = \sum_{qp} \mathbf{N}^\mathsf{T} \mathbf{f}\, J_0 w.

Surface tractions and pressures on a boundary face :math:`\Gamma_e` are computed as

.. math::

   \mathbf{P}_e^{(t)} = \int_{\Gamma_e} \mathbf{N}^\mathsf{T} \, \mathbf{t} \, \mathrm{d}\Gamma, \qquad
   \mathbf{P}_e^{(p)} = - \int_{\Gamma_e} p \, \mathbf{N}^\mathsf{T} \, \mathbf{n} \, \mathrm{d}\Gamma.

For 2D and 1D problems, the integration measure :math:`J_0 w` is scaled by thickness
:math:`t` or cross-sectional area :math:`A`, respectively: :math:`J_0 w \leftarrow J_0 w\, t`
(2D), :math:`J_0 w \leftarrow J_0 w\, A` (1D).

Constitutive updates are carried out per quadrature point using a Marmot material model.

- 3D solid: full 3D update returning :math:`\boldsymbol{\sigma}` and :math:`\mathbf{C}`.
  
- Plane stress: in-plane update consistent with :math:`\sigma_{zz}=0`.
  
- Plane strain: 3D update with reduction to a plane-strain tangent.
  
- 1D uniaxial stress: scalar update consistent with axial behavior.

Each quadrature point stores stress, strain and a material state vector; the element
accumulates the history incrementally and may suggest a reduced time step if required
by the material.

Implementation
--------------

.. doxygenclass:: Marmot::Elements::DisplacementFiniteElement
   :allow-dot-graphs:
