MPM Cells and Material Points
=============================

Overview
--------

In the material point method (MPM), the continuum is represented by **material points**, which carry mass, position
and the complete history of the material, and are advected through a fixed background grid of **cells**, which carry
the nodal fields and the interpolation. In every increment, each cell integrates the residual and the tangent over the
material points it currently hosts, and the grid is reset afterwards. Marmot provides the interfaces of both building
blocks, the cell geometries, and the factories through which a host framework (e.g. EdelweissMeshfree) creates them by
name; the concrete formulations are described in :doc:`displacementcell`, :doc:`displacementmaterialpoint`,
:doc:`gradientenhancedfinitestraincell` and :doc:`gradientenhancedfinitestrainmaterialpoint`. For the overall
architecture and the naming convention see :doc:`meshfree`.

All classes of this page are part of the module ``core/MarmotMeshfreeCore``.

Cells
-----

A cell (``MarmotCell``) is a background element with :math:`n_\mathrm{nodes}` nodes and shape functions
:math:`N_A`. It owns no quadrature points: the host assigns to it the material points :math:`p` located inside it, and
the cell evaluates its integrals as sums over those material points. For the finite-strain cells of Marmot, the grid
is in the configuration :math:`\boldsymbol{Y}` at the beginning of the increment, the nodal unknowns are the increments
:math:`\Delta\boldsymbol{Q}_A` of the current step, and the momentum balance reads, e.g. for
:doc:`displacementcell`,

.. math::

   \boldsymbol{r}_{A} = \sum_{p} \boldsymbol{\tau}_p\,\nabla_x N_A(\boldsymbol{Y}_p)\; V^0_p ,
   \qquad \nabla_x N_A = \Delta\boldsymbol{F}^{-\mathsf T}_p\,\nabla_Y N_A ,

with the Kirchhoff stress :math:`\boldsymbol{\tau}_p`, the incremental deformation gradient
:math:`\Delta\boldsymbol{F}_p` and the reference volume :math:`V^0_p` of the material point.

**Grid-to-particle map.** ``interpolateFieldsToMaterialPoints`` evaluates the increments at the material points,

.. math::

   \Delta\boldsymbol{u}_p = \sum_A N_A(\boldsymbol{\xi}_p)\,\Delta\boldsymbol{Q}_A ,
   \qquad
   \frac{\partial\Delta\boldsymbol{u}}{\partial\boldsymbol{Y}}\Big|_p = \sum_A \Delta\boldsymbol{Q}_A \otimes \nabla_Y N_A(\boldsymbol{\xi}_p) ,

and passes them to the material point, which accumulates them into its kinematic increment. The parametric
coordinates :math:`\boldsymbol{\xi}_p` and the shape function values and gradients at each material point are computed
once in ``assignMaterialPoints``.

**Dof layout.** A cell has ``getNDofPerCell()`` dofs. ``getNodeFields()`` lists for each node the names of the fields
it carries (e.g. ``displacement``, ``nonlocal damage``). Internally, the cells of Marmot use a *blocked* layout: all dofs
of the first field (node by node, component fastest), then all dofs of the next field.
``getDofIndicesPermutationPattern()`` maps the node-major order of the host to this internal order.

**Loads.** Body loads (``BODYFORCE``) are integrated over the hosted material points; distributed loads
(``PRESSURE``) are carried by an individual material point, selected by its number. The supported load names are
reported by ``getSupportedBodyLoadTypes()`` and ``getSupportedDistributedLoadTypes()``. Lumped and consistent inertia
are available for dynamic analyses.

Call sequence
~~~~~~~~~~~~~

Per time increment, the host calls:

.. list-table::
   :header-rows: 1

   * - Phase
     - Material points
     - Cells
   * - once, at the start
     - ``assignMaterial``, ``assignStateVars``, ``initializeYourself``
     -
   * - connectivity update
     - ``prepareYourself``
     - ``isCoordinateInCell`` / ``getBoundingBox`` (search), ``assignMaterialPoints`` (active cells)
   * - every iteration
     - ``prepareYourself``, then (after the cells) ``computeYourself``
     - ``interpolateFieldsToMaterialPoints``, then (after the material points) ``computeMaterialPointKernels``,
       ``computeBodyLoad``, ``computeDistributedLoad``, inertia
   * - after convergence
     - ``acceptStateAndPosition``
     -

``prepareYourself`` resets the accumulated kinematic increment, so it must precede the interpolation in every
iteration. A cell and a material point exchange formulation-specific data (the kinematic increment, the stress and its
tangent) that is not part of the abstract interfaces; a cell therefore accepts only material points of its own
formulation and throws ``std::invalid_argument`` in ``assignMaterialPoints`` otherwise.

Cell geometries
---------------

The formulation cells are templates over a *geometry policy*, a class that satisfies the concept
``GeometryCellPolicy``: it provides the types ``XiSized``, ``NSized``, ``dNdXSized`` and the methods
``findReferenceCoordinate`` (inverse map), ``N``, ``dNdX`` and ``detJ``; the cells additionally use
``isCoordinateInCell``, ``getBoundingBox`` and ``getElementShape``. (The virtual interface ``MarmotCellGeometry``
describing the same contract is deprecated and unused.) Two policies are available.

**Lagrangian cells** (``MarmotLagrangianCellGeometry``): isoparametric Quad4 and Hexa8 cells with Abaqus node order,

.. math::

   \boldsymbol{X}(\boldsymbol{\xi}) = \sum_A N_A(\boldsymbol{\xi})\,\boldsymbol{X}_A , \qquad
   \boldsymbol{\xi}\in[-1,1]^{n_\mathrm{dim}}, \qquad
   \nabla_X N_A = \boldsymbol{J}^{-\mathsf T}\,\nabla_\xi N_A, \quad J_{ij} = \frac{\partial X_i}{\partial \xi_j} .

The point location test is the half-open bounding box test :math:`X_{\min,i} \le x_i < X_{\max,i}`, so that a point on
a shared face belongs to exactly one cell.

**B-spline cells** (``MarmotBSplineCellGeometry`` on top of ``MarmotBSplineGeometryElement``): a cell is one knot span
:math:`[u_p, u_{p+1}]` per direction of a tensor-product B-spline of degree :math:`p`. The
:math:`(p+1)^{n_\mathrm{dim}}` basis functions that are nonzero on the span follow from the :math:`2p+2` knots
:math:`u_0,\dots,u_{2p+1}` per direction by the Cox--de Boor recursion (see :doc:`meshfreeapproximation`),

.. math::

   N_a(\boldsymbol{\xi}) = B_{i,p}(\xi_1)\,B_{j,p}(\xi_2)\,B_{k,p}(\xi_3), \qquad a = i + j\,(p+1) + k\,(p+1)^2 .

The knot vectors are given in physical coordinates, and the parametric coordinates are the physical coordinates:
``findReferenceCoordinate`` is the identity and ``dNdX`` returns the parametric derivatives. This corresponds to a
background grid whose B-spline geometry map is the identity (e.g. a uniform, axis-aligned grid). The cell's bounding
box is the knot span. Compared to Lagrangian cells, B-spline cells of degree :math:`p \ge 2` have
:math:`C^{p-1}`-continuous shape functions across cell boundaries, which mitigates the cell-crossing error of MPM.

``MarmotLagrangeCell`` (aliases ``Marmot::Meshfree::Quad4`` and ``Hex8``, header ``MarmotMeshfreeQuadHexCell.h``) is a
self-contained Quad4/Hexa8 geometry with centroid, volume, second moments, face centers, boundary surface vectors and
uniform subdivision, evaluated with a :math:`2^{n_\mathrm{dim}}`-point Gauss rule. It is used for the integration
domains of particles (see :doc:`meshfreeparticles`), not by the MPM cells.

Cell elements
-------------

A cell element (``MarmotCellElement``) is a cell that defines its own material points, as a finite element defines its
quadrature points: the host queries ``getNMaterialPoints()``, ``getRequestedMaterialPointCoordinates()`` and
``getRequestedMaterialPointVolumes()``, creates one material point per request and assigns them back with
``assignMaterialPoints``. The quadrature rule is passed on creation; a cell element may carry state variables of its
own. No concrete cell element is currently contained in Marmot.

Material points
---------------

A material point (``MarmotMaterialPoint``) owns its constitutive model (``assignMaterial`` with a
``MarmotMaterialSection``) and works on an externally owned state variable array (``getNumberOfRequiredStateVars``,
``assignStateVars``), which holds the state of the formulation (e.g. displacement, velocity, acceleration, deformation
gradient and its increment, stress) followed by the state variables of the material; entries are accessible by name via ``getStateView``. Its
kinematic update is incremental: ``prepareYourself`` sets :math:`\Delta\boldsymbol{u} = \boldsymbol{0}` and
:math:`\Delta\boldsymbol{F} = \boldsymbol{I}`, the cells add their contributions, ``computeYourself`` evaluates the
material at :math:`\boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n` and provides the stress and the tangent
to the cells, and ``acceptStateAndPosition`` commits
:math:`\boldsymbol{u}_{n+1} = \boldsymbol{u}_n + \Delta\boldsymbol{u}`,
:math:`\boldsymbol{F}_{n+1} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n`.

The position reported by ``getCoordinatesAtCenter`` and ``getVertexCoordinates`` is the last accepted one,
:math:`\boldsymbol{X} + \boldsymbol{u}_n`; the cells locate the material points with it. For mass and inertia, a
material point reports its reference volume :math:`V^0` and density :math:`\rho_0`.

Factories and registered types
------------------------------

Cells, cell elements and material points register themselves at static initialization with a factory function in
their module's ``*Registration.cpp``, and the host creates them by name (case-insensitive, stored in upper case):

- ``MarmotLibrary::MarmotMaterialPointFactory::createMaterialPoint(name, number, vertexCoordinates, size, volume)``,
- ``MarmotLibrary::MarmotCellFactory::createCell(name, number, nodeCoordinates, size)`` for Lagrangian cells and
  ``createBSplineCell(name, number, controlPointCoordinates, size, knotVectors, sizeKnotVectors)`` for B-spline cells
  (two separate registries),
- ``MarmotLibrary::MarmotCellElementFactory::createCellElement(name, number, nodeCoordinates, size, quadratureRule,
  quadratureOrder)``.

An unknown name throws ``std::invalid_argument``; the caller owns the returned object. The cells map (do not copy) the
coordinate array passed on creation, which must therefore outlive the cell.

Registered cells (``<Formulation>`` is ``Displacement`` or ``GradientEnhancedFiniteStrain``):

.. list-table::
   :header-rows: 1

   * - Name
     - Factory
     - Dimension
     - Geometry
     - Nodes
   * - ``<Formulation>/Quad4``
     - ``createCell``
     - 2
     - Lagrangian Quad4
     - 4
   * - ``<Formulation>/Hexa8``
     - ``createCell``
     - 3
     - Lagrangian Hexa8
     - 8
   * - ``<Formulation>/BSpline/1``, ``/2``, ``/3``
     - ``createBSplineCell``
     - 2
     - B-spline, degree 1, 2, 3
     - 4, 9, 16
   * - ``<Formulation>/BSpline/3D/1``, ``/2``, ``/3``
     - ``createBSplineCell``
     - 3
     - B-spline, degree 1, 2, 3
     - 8, 27, 64

Registered material points:

.. list-table::
   :header-rows: 1

   * - Name
     - Dimension
   * - ``Displacement/PlaneStrain``
     - 2 (plane strain)
   * - ``Displacement/3D``
     - 3
   * - ``GradientEnhancedFiniteStrain/PlaneStrain``
     - 2 (plane strain)
   * - ``GradientEnhancedFiniteStrain/3D``
     - 3

A cell must be combined with a material point of the same formulation and dimension.

Limitations
-----------

- The inverse map of the Lagrangian cells (``MarmotLagrangianCellGeometry::findReferenceCoordinate``) uses the affine
  map of the cell's bounding box and has no Newton update yet; together with the bounding-box point location test, the
  Lagrangian cells are therefore only valid for **axis-aligned, box-shaped** cells. For a distorted cell the inverse
  map throws ``std::runtime_error``.
- The B-spline cells assume an identity geometry map (see above); the control point coordinates enter only ``detJ``.
- Cubic B-spline cells have no output shape name: ``getCellShape()`` returns an empty string for
  ``.../BSpline/3`` and ``.../BSpline/3D/3`` (``Quad16`` and ``Hexa64`` are missing in the shape table of
  ``MarmotBSplineGeometryElement``).
- ``MarmotBSplineGeometryElement`` is usable in 2D and 3D only: in 1D, ``dNdXi`` would not compile if instantiated,
  and the constructor throws for degrees above 1.

Implementation
--------------

.. doxygenclass:: MarmotCell
   :allow-dot-graphs:

.. doxygenclass:: MarmotCellElement
   :allow-dot-graphs:

.. doxygenclass:: MarmotMaterialPoint
   :allow-dot-graphs:

.. doxygenclass:: MarmotLagrangianCellGeometry
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Cells::MarmotBSplineCellGeometry
   :allow-dot-graphs:

.. doxygenclass:: Marmot::FiniteElement::MarmotBSplineGeometryElement
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::MarmotLagrangeCell
   :allow-dot-graphs:

.. doxygenclass:: MarmotCellGeometry
   :allow-dot-graphs:

.. doxygenclass:: MarmotLibrary::MarmotMaterialPointFactory
   :allow-dot-graphs:

.. doxygenclass:: MarmotLibrary::MarmotCellFactory
   :allow-dot-graphs:

.. doxygenclass:: MarmotLibrary::MarmotCellElementFactory
   :allow-dot-graphs:
