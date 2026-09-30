Meshfree methods
================

Marmot provides the constitutive and discretization kernels for two meshfree methods, which a host framework
(e.g. `EdelweissMeshfree <https://github.com/Edelweiss-Numerics/EdelweissMeshfree>`_) combines into a simulation:

- the **material point method (MPM)**: *material points* carry the material state and move through a background grid
  of *cells*; each cell assembles the residual and the tangent of the material points it currently hosts;
- the **reproducing kernel particle method (RKPM)**: *particles* carry the material state and a domain of integration;
  their shape functions follow from a *reproducing kernel approximation* built on *kernel functions* attached to the
  nodes.

Both methods reuse Marmot's material interfaces (e.g. ``MarmotMaterialFiniteStrain`` and
``MarmotMaterialGradientEnhancedFiniteStrain``), so that every registered material of such an interface can be used in
them.

Architecture
------------

The meshfree layer consists of a core module and of formulation modules in the dedicated module categories:

.. list-table::
   :header-rows: 1

   * - Module category
     - Content
     - Documentation
   * - ``core/MarmotMeshfreeCore``
     - interfaces (``MarmotCell``, ``MarmotMaterialPoint``, ``MarmotParticle``, ``MarmotCellElement``), cell geometries,
       the reproducing kernel approximation, kernel functions, generic particle implementations and the factories
     - :doc:`meshfreeapproximation`, :doc:`meshfreecells`, :doc:`meshfreeparticles`
   * - ``materialpoints/``
     - MPM material points of a formulation
     - :doc:`displacementmaterialpoint`, :doc:`gradientenhancedfinitestrainmaterialpoint`
   * - ``cells/``
     - MPM cells of a formulation, which assemble the contributions of their material points
     - :doc:`displacementcell`, :doc:`gradientenhancedfinitestraincell`
   * - ``particles/``
     - RKPM particles of a formulation, with their integration schemes
     - :doc:`displacementparticle`, :doc:`gradientenhancedfinitestrainparticle`

Two formulations are available, each for MPM (material point and cell) and RKPM (particle):

- **Displacement**: finite-strain displacement formulation, for materials of the ``MarmotMaterialFiniteStrain``
  interface;
- **GradientEnhancedFiniteStrain**: finite-strain displacement formulation coupled to a scalar nonlocal field
  (implicit-gradient regularization), for materials of the ``MarmotMaterialGradientEnhancedFiniteStrain`` interface,
  e.g. :doc:`gradientenhancedcompressibleneohookedamage` and :doc:`gradientenhancedfinitestraindruckerprager`. It is
  the meshfree counterpart of :doc:`gradientenhancedfinitestraindisplacementelement`.

Registration and naming
-----------------------

As materials and elements, all meshfree types register themselves by name with a factory of
``MarmotMeshfreeCore`` (``MarmotMaterialPointFactory``, ``MarmotCellFactory``, ``MarmotCellElementFactory``,
``MarmotParticleFactory``); the host creates them by that name. The names follow the pattern

.. list-table::
   :header-rows: 1

   * - Type
     - Pattern
     - Examples
   * - material point
     - ``<Formulation>/<Dimension>``
     - ``Displacement/PlaneStrain``, ``GradientEnhancedFiniteStrain/3D``
   * - cell (Lagrangian)
     - ``<Formulation>/<Shape>``
     - ``Displacement/Quad4``, ``GradientEnhancedFiniteStrain/Hexa8``
   * - cell (B-spline)
     - ``<Formulation>/BSpline[/3D]/<Order>``
     - ``Displacement/BSpline/2``, ``GradientEnhancedFiniteStrain/BSpline/3D/1``
   * - particle
     - ``<Formulation><IntegrationScheme>/<Dimension>/<Shape>``
     - ``GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad``, ``DisplacementSQCNI/3D/Hexa``,
       ``Displacement/PlaneStrain/Point``

The integration schemes of the particles (SQCNI, SNNI, NSNI, SDI) and the corrections VCI and CWF are defined in
:doc:`meshfreeparticles`; the registered names of each formulation are listed on its pages.

.. note::
   The displacement NSNI particles were formerly registered as ``Displacement/SQCNIxNSNI/...`` (with a slash after
   the formulation) and ``Displacement/R-SNNIxNSNI/...``; these names remain as aliases, see
   :doc:`displacementparticle`.

Usage with EdelweissMeshfree
----------------------------

EdelweissMeshfree creates the Marmot types by their registered names:

- **MPM**: the cells through the cell provider ``LagrangianMarmotCell`` (or ``BSplineMarmotCell``) with, e.g.,
  ``cellType="GradientEnhancedFiniteStrain/Quad4"``, the material points through the provider ``marmot`` with, e.g.,
  ``mpType="GradientEnhancedFiniteStrain/PlaneStrain"``, and the material by its name and properties;
- **RKPM**: the kernel functions through ``MarmotMeshfreeKernelFunctionWrapper`` (e.g. ``"BSplineBoxed"`` with a support
  radius; its ``continuityOrder`` argument selects the 2nd- or 3rd-order kernel, see :doc:`meshfreeapproximation`), the approximation through ``MarmotMeshfreeApproximationWrapper``, and the particles
  through ``MarmotParticleWrapper`` with the particle name, e.g.
  ``"GradientEnhancedFiniteStrainSQCNI/PlaneStrain/Quad"``.

The EdelweissMeshfree examples ``160_gradient_enhanced_finite_strain_mpm_test`` (MPM) and
``161_gradient_enhanced_finite_strain_rkpm_test`` (RKPM) are complete, verified input scripts for both methods.

.. toctree::
  :maxdepth: 1

  meshfreeapproximation
  meshfreecells
  meshfreeparticles
  displacementmaterialpoint
  displacementcell
  displacementparticle
  gradientenhancedfinitestrainmaterialpoint
  gradientenhancedfinitestraincell
  gradientenhancedfinitestrainparticle
