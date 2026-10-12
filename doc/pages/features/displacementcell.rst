Displacement Cell
=================

Theory
------

The displacement cell is the MPM background cell of the finite-strain displacement formulation. It carries the
displacement field (``nDim`` dofs per node) and assembles the residual and the tangent of the
:doc:`displacementmaterialpoint` instances it currently hosts. Cell geometries (Lagrangian and B-spline), the cell
interface and the MPM workflow are described in :doc:`meshfreecells`; the RKPM counterpart of this formulation is the
:doc:`displacementparticle` (see :doc:`meshfreeparticles` and :doc:`meshfreeapproximation`).

As usual in MPM, the grid dofs :math:`\Delta q_{Bk}` are the **increments** of the current step, and the material points
carry the accumulated state. The cell nodes define the intermediate configuration :math:`\boldsymbol{Y}` (the
configuration at the beginning of the increment). When the material points are assigned, the cell evaluates at the
position :math:`\boldsymbol{Y}_p` of each material point the shape functions :math:`N_A` and their gradients
:math:`\partial N_A/\partial Y_J`, and caches them. Each material point then receives

.. math::

   \Delta\boldsymbol{u}_p = N_B\,\Delta\boldsymbol{q}_B, \qquad
   \Delta F_{iJ} = \delta_{iJ} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial Y_J}.

The momentum balance is assembled with the Kirchhoff stress :math:`\boldsymbol{\tau}` of the material points, the
spatial gradients :math:`\partial N_A/\partial x_i = \Delta F^{-1}_{Ji}\,\partial N_A/\partial Y_J` and the undeformed
volumes :math:`V_p^0` of the material points (the weak form in the undeformed configuration, evaluated at the material
points),

.. math::

   r_{Aj} = \sum_p \frac{\partial N_A}{\partial x_i}\,\tau_{ij}\,V_p^0 ,

with the consistent tangent, consisting of the material part and the geometric stiffness,

.. math::

   \frac{\partial r_{Aj}}{\partial \Delta q_{Bk}} = \sum_p \left(
     \frac{\partial N_A}{\partial x_i}\,\frac{\partial \tau_{ij}}{\partial \Delta F_{kL}}\,
     \frac{\partial N_B}{\partial Y_L}
     - \frac{\partial N_A}{\partial x_k}\,\tau_{ij}\,\frac{\partial N_B}{\partial x_i} \right) V_p^0 ,

where :math:`\partial\tau_{ij}/\partial\Delta F_{kL}` is provided by the material point. The cell does not contain an
inertia term in the residual; it provides the consistent mass matrix
:math:`M_{AiBi} = \sum_p N_A\,N_B\,\rho_0\,V_p^0` and the lumped mass (its row sums) to the host.

Loads
-----

.. list-table::
   :header-rows: 1

   * - Name
     - Type
     - Contribution
   * - ``BODYFORCE``
     - body load
     - :math:`r_{Aj} \mathrel{-}= \sum_p N_A\,f_j\,V_p^0`, body force :math:`\boldsymbol{f}` per undeformed volume
       (no tangent)
   * - ``PRESSURE``
     - distributed load at one material point
     - the given load vector :math:`\boldsymbol{f}_0` (e.g. :math:`p\,\boldsymbol{N}\,dA_0`) is transformed by
       Nanson's formula with the total deformation gradient,
       :math:`\boldsymbol{f} = J\,\boldsymbol{F}^{-\mathsf T}\boldsymbol{f}_0`, and assembled as
       :math:`r_{Aj} \mathrel{-}= N_A\,f_j` with its load stiffness

Registered cells
----------------

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension
     - Geometry
   * - ``Displacement/Quad4``
     - ``LagrangianDisplacementCell<2, 4>``
     - 2D, plane strain
     - bilinear quadrilateral
   * - ``Displacement/Hexa8``
     - ``LagrangianDisplacementCell<3, 8>``
     - 3D
     - trilinear hexahedron
   * - ``Displacement/BSpline/1``, ``/2``, ``/3``
     - ``BSplineDisplacementCell<2, 4|9|16, 1|2|3>``
     - 2D, plane strain
     - B-spline of order 1, 2, 3
   * - ``Displacement/BSpline/3D/1``, ``/2``, ``/3``
     - ``BSplineDisplacementCell<3, 8|27|64, 1|2|3>``
     - 3D
     - B-spline of order 1, 2, 3

The node field is ``displacement``; the dofs are ordered node by node. The cells have no properties; the material is
assigned to the material points.

Implementation
--------------

.. doxygenclass:: Marmot::Cells::DisplacementCell
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Cells::LagrangianDisplacementCell
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Cells::BSplineDisplacementCell
   :allow-dot-graphs:
