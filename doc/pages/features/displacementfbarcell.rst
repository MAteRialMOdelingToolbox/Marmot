Displacement F-bar Cell
=======================

Theory
------

The displacement F-bar cell is the :doc:`displacementcell` with an F-bar projection of the volumetric deformation,
against volumetric locking of nearly incompressible materials and against the spurious volumetric response of cells
that host few, unevenly distributed material points. It hosts the same :doc:`displacementmaterialpoint` instances, has
the same dofs, loads and mass, and is a separate module (``cells/DisplacementFBarCell``, which requires
``DisplacementCell``).

The volume change of the increment is evaluated at the geometric center of the hosted material points, the mean
:math:`\bar{\boldsymbol{\xi}} = \frac{1}{n_p}\sum_p \boldsymbol{\xi}_p` of their parametric coordinates,

.. math::

   \Delta\boldsymbol{F}_c = \boldsymbol{I} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial \boldsymbol{Y}}\Big|_{\bar{\boldsymbol{\xi}}},
   \qquad
   \Delta\bar{\boldsymbol{F}}_p = s_p\,\Delta\boldsymbol{F}_p ,
   \qquad
   s_p = \left(\frac{\det\Delta\boldsymbol{F}_c}{\det\Delta\boldsymbol{F}_p}\right)^{1/n_\mathrm{dim}} ,

and each material point receives :math:`\Delta\bar{\boldsymbol{F}}_p` instead of :math:`\Delta\boldsymbol{F}_p`, so that
its accumulated deformation gradient is :math:`\bar{\boldsymbol{F}} = \Delta\bar{\boldsymbol{F}}\,\bar{\boldsymbol{F}}_n`.
The projection is incremental, since the cell has no history of its own. In plane strain, the exponent is :math:`1/2`
and the out-of-plane stretch remains 1.

The residual is that of the :doc:`displacementcell` with :math:`\boldsymbol{\tau}(\bar{\boldsymbol{F}})` and the spatial
gradients of the unmodified kinematics, :math:`\partial N_A/\partial x_i = \Delta F^{-1}_{p,Ji}\,\partial N_A/\partial Y_J`.
The tangent is exact: in its material part,
:math:`\delta_{km}\,\partial N_B/\partial Y_L` is replaced by

.. math::

   \frac{\partial \Delta\bar F_{kL}}{\partial \Delta q_{Bm}}
     = s_p\,\delta_{km}\,\frac{\partial N_B}{\partial Y_L}
     + \frac{1}{n_\mathrm{dim}}\,\Delta\bar F_{kL}
       \left(\frac{\partial N^c_B}{\partial x_m} - \frac{\partial N_B}{\partial x_m}\right),
   \qquad
   \frac{\partial N^c_B}{\partial x_m} = \Delta F^{-1}_{c,Lm}\,\frac{\partial N_B}{\partial Y_L}\Big|_{\bar{\boldsymbol{\xi}}} ,

and the geometric part is unchanged. The ``PRESSURE`` load is transformed with :math:`\bar{\boldsymbol{F}}`, with the same
chain rule in its load stiffness.

With a single material point in a cell, the center is the point itself, :math:`s_p \equiv 1`, and the cell coincides
with the :doc:`displacementcell`; the F-bar projection cannot stabilize such a cell (see the aggregated MPM and the basis
extrapolation in :doc:`meshfreecells`). For a homogeneous increment the residual coincides with that of the
:doc:`displacementcell`.

Registered cells
----------------

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension
     - Geometry
   * - ``DisplacementFBar/Quad4``
     - ``LagrangianDisplacementFBarCell<2, 4>``
     - 2D, plane strain
     - bilinear quadrilateral
   * - ``DisplacementFBar/Hexa8``
     - ``LagrangianDisplacementFBarCell<3, 8>``
     - 3D
     - trilinear hexahedron
   * - ``DisplacementFBar/BSpline/1``, ``/2``, ``/3``
     - ``BSplineDisplacementFBarCell<2, 4|9|16, 1|2|3>``
     - 2D, plane strain
     - B-spline of order 1, 2, 3
   * - ``DisplacementFBar/BSpline/3D/1``, ``/2``, ``/3``
     - ``BSplineDisplacementFBarCell<3, 8|27|64, 1|2|3>``
     - 3D
     - B-spline of order 1, 2, 3

The cells are combined with the material points ``Displacement/PlaneStrain`` (2D) and ``Displacement/3D``.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Cells::DisplacementFBarCell
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Discretizations::Cells::LagrangianDisplacementFBarCell
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Discretizations::Cells::BSplineDisplacementFBarCell
   :allow-dot-graphs:
