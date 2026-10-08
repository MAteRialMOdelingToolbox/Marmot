Gradient-Enhanced Finite-Strain Cell
====================================

Theory
------
The gradient-enhanced finite-strain cell is the MPM background cell of the implicit-gradient finite-strain formulation.
It carries the displacement field :math:`\boldsymbol u` (``nDim`` dofs per node) and a scalar nonlocal field
:math:`\bar N` (1 dof per node) and assembles the contributions of the
:doc:`gradientenhancedfinitestrainmaterialpoint` instances it currently hosts, which in turn drive a material of the
interface ``MarmotMaterialGradientEnhancedFiniteStrain`` (e.g. :doc:`gradientenhancedcompressibleneohookedamage`,
:doc:`gradientenhancedfinitestraindruckerprager`). It is the MPM counterpart of
:doc:`gradientenhancedfinitestraindisplacementelement`; the cell geometries and the MPM workflow are described in
:doc:`meshfreecells`.

Kinematics
^^^^^^^^^^
As usual in MPM, the grid dofs are the *increments* :math:`\Delta\boldsymbol q` of the current step, while the material
points carry the accumulated state. When the material points are assigned, the cell evaluates its shape functions
:math:`N_A` and their gradients :math:`\nabla_Y N_A` at the position of each material point in the intermediate
reference configuration :math:`\boldsymbol Y` (the last accepted configuration). The increment is interpolated to each
point (``interpolateFieldsToMaterialPoints``),

.. math::

   \Delta\boldsymbol u_p = N_A\,\Delta\boldsymbol q^U_A, \qquad
   \frac{\partial\Delta u_i}{\partial Y_j} = \Delta q^U_{Ai}\,\frac{\partial N_A}{\partial Y_j}, \qquad
   \Delta\bar N_p = N_A\,\Delta q^N_A,

and the spatial and the undeformed gradients follow with the increment :math:`\Delta\boldsymbol F` and the accepted
deformation gradient :math:`\boldsymbol F_n` of the material point,

.. math::

   \frac{\partial N_A}{\partial x_i} = \frac{\partial N_A}{\partial Y_j}\,\Delta F^{-1}_{ji}, \qquad
   \frac{\partial N_A}{\partial X_i} = \frac{\partial N_A}{\partial Y_j}\,F_{n,ji} .

Weak forms
^^^^^^^^^^
Summing over the material points :math:`p` with their undeformed volumes :math:`V^0_p`, the momentum balance is
assembled with the Kirchhoff stress,

.. math::

   r_{U,Aj} = \sum_p \frac{\partial N_A}{\partial x_i}\,\tau_{ij}\,V^0_p ,

and the nonlocal balance :math:`\bar N - c\,\nabla_X^2\bar N = L` in the undeformed configuration, in increment form,

.. math::

   r_{N,A} = \sum_p \Bigl( N_A\,\Delta\bar N_p + c\,\nabla_X N_A\cdot\nabla_X\Delta\bar N_p - N_A\,\Delta L_p \Bigr)
   V^0_p ,

with :math:`\nabla_X\Delta\bar N_p = \nabla_X N_B\,\Delta q^N_B`, the change :math:`\Delta L_p` of the local driving force
reported by the material point and :math:`c = R^2` from the material. The undeformed gradient :math:`\nabla_X` does not
depend on the current increment. The residual is the internal force vector; body forces and pressures are added with
a negative sign.

Tangent
^^^^^^^
With the tangents of the material point with respect to :math:`\Delta\boldsymbol F` and :math:`\bar N`, the blocks of
the consistent tangent are

.. math::

   K^{UU}_{jAkB} &= \sum_p \Bigl( \frac{\partial N_A}{\partial x_i}\,\frac{\partial\tau_{ij}}{\partial\Delta F_{kL}}\,
   \frac{\partial N_B}{\partial Y_L} - \frac{\partial N_A}{\partial x_k}\,\tau_{ij}\,
   \frac{\partial N_B}{\partial x_i} \Bigr) V^0_p ,\\
   K^{UN}_{jAB} &= \sum_p \frac{\partial N_A}{\partial x_i}\,\frac{\partial\tau_{ij}}{\partial\bar N}\,N_B\,V^0_p ,\\
   K^{NU}_{AkB} &= -\sum_p N_A\,\frac{\partial L}{\partial\Delta F_{kL}}\,\frac{\partial N_B}{\partial Y_L}\,V^0_p ,\\
   K^{NN}_{AB} &= \sum_p \Bigl( N_A N_B\,\bigl(1 - \tfrac{\partial L}{\partial\bar N}\bigr)
   + c\,\nabla_X N_A\cdot\nabla_X N_B \Bigr) V^0_p ,

where the second term of :math:`K^{UU}` is the geometric stiffness from the dependence of :math:`\nabla_x N_A` on
:math:`\Delta\boldsymbol F`. The interaction :math:`c` is treated as constant.

Loads and inertia
^^^^^^^^^^^^^^^^^

- ``BODYFORCE``: :math:`f_{Ai} \mathrel{-}= \sum_p N_A\,b_i\,V^0_p` with the body force :math:`\boldsymbol b` per unit
  undeformed volume, no load stiffness.
- ``PRESSURE``: a follower load acting on a single material point (selected by its label). The undeformed load vector
  :math:`\boldsymbol f_0 = p\,\boldsymbol N\,dA_0` is pushed forward by Nanson's formula with the total deformation
  gradient :math:`\boldsymbol F = \Delta\boldsymbol F\,\boldsymbol F_n`,
  :math:`\boldsymbol f = J\,\boldsymbol F^{-\mathsf T}\boldsymbol f_0`, and applied as
  :math:`f_{Ai} \mathrel{-}= N_A f_i`, with the consistent load stiffness.
- Inertia: the consistent mass :math:`M_{AiBi} = \sum_p N_A N_B\,\rho_0 V^0_p` of the displacement field and its row-sum
  lumped version. The nonlocal field carries no inertia.

The dofs are ordered field by field inside the cell (all displacements, then all nonlocal dofs);
``getDofIndicesPermutationPattern`` maps them to the node-by-node layout of the host.

Registered types
----------------

The cells are registered in the ``MarmotCellFactory`` (Lagrangian, ``createCell``) and as B-spline cells
(``createBSplineCell``). The nodal fields are ``displacement`` and ``nonlocal damage``.

.. list-table::
   :header-rows: 1

   * - Name
     - Geometry
     - Dimension
     - Nodes
   * - ``GradientEnhancedFiniteStrain/Quad4``
     - Lagrangian, bilinear quadrilateral
     - 2 (plane strain)
     - 4
   * - ``GradientEnhancedFiniteStrain/Hexa8``
     - Lagrangian, trilinear hexahedron
     - 3
     - 8
   * - ``GradientEnhancedFiniteStrain/BSpline/1``
     - B-spline, order 1
     - 2 (plane strain)
     - 4
   * - ``GradientEnhancedFiniteStrain/BSpline/2``
     - B-spline, order 2
     - 2 (plane strain)
     - 9
   * - ``GradientEnhancedFiniteStrain/BSpline/3``
     - B-spline, order 3
     - 2 (plane strain)
     - 16
   * - ``GradientEnhancedFiniteStrain/BSpline/3D/1``
     - B-spline, order 1
     - 3
     - 8
   * - ``GradientEnhancedFiniteStrain/BSpline/3D/2``
     - B-spline, order 2
     - 3
     - 27
   * - ``GradientEnhancedFiniteStrain/BSpline/3D/3``
     - B-spline, order 3
     - 3
     - 64

A cell accepts only material points of the matching dimension, ``GradientEnhancedFiniteStrain/PlaneStrain`` or
``GradientEnhancedFiniteStrain/3D`` (see :doc:`gradientenhancedfinitestrainmaterialpoint`); any other material point
is rejected by ``assignMaterialPoints``.

Implementation
--------------

.. doxygenclass:: Marmot::Cells::GradientEnhancedFiniteStrainCell
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Cells::LagrangianGradientEnhancedFiniteStrainCell
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Cells::BSplineGradientEnhancedFiniteStrainCell
   :allow-dot-graphs:
