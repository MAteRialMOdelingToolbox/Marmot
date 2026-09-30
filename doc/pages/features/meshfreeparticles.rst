Meshfree Particles
==================

Overview
--------

In the reproducing kernel particle method (RKPM), a *particle* is the integration entity: it represents a piece of the
body, carries the material state, and evaluates the weak form on its piece. Its shape functions are those of a
reproducing kernel approximation built on the *kernel functions* of the nodes (see :doc:`meshfreeapproximation`).
Unlike a finite element, a particle has no fixed nodes: the host framework (e.g., EdelweissMeshfree) determines the
kernel functions whose support covers the particle and hands them over; their number, and hence the size of the
particle's residual vector and stiffness matrix, may change from increment to increment.

``MarmotMeshfreeCore`` provides the physics-independent part of the particles:

- :cpp:class:`Marmot::Meshfree::MarmotParticle`, the abstract interface called by the host framework;
- :cpp:class:`Marmot::Meshfree::ParticleDomain`, the geometry of a cell-shaped particle (quadrilateral or hexahedron)
  with its smoothing domain and its subdivision;
- :cpp:class:`Marmot::Meshfree::GenericParticle`, the base of the point particles, including the variationally
  consistent integration (VCI);
- :cpp:class:`Marmot::Meshfree::GenericSDIParticle`, the base of the particles with subdomain integration (SDI);
- :cpp:class:`MarmotLibrary::MarmotParticleFactory`, the registry by which particles are created by name.

The concrete particles of a formulation live in the ``particles/`` module category, see :doc:`displacementparticle`
and :doc:`gradientenhancedfinitestrainparticle`; the overall architecture is described in :doc:`meshfree`.

Throughout, :math:`\boldsymbol{X}` denotes the undeformed configuration and :math:`\boldsymbol{Y}` the *reference
intermediate configuration*, i.e., the configuration of the last accepted (converged) increment. All particles work
in this updated Lagrangian setting: shape function gradients are taken with respect to :math:`\boldsymbol{Y}`, and at
the end of every increment the particle is moved to the new configuration.

The particle interface
----------------------

The following methods of :cpp:class:`Marmot::Meshfree::MarmotParticle` are called by the host framework, in this
order:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Method
     - Purpose
   * - ``MarmotParticleFactory::createParticle``
     - create the particle by name from its vertex coordinates (the center for point particles), its volume (point
       particles; zero for cell-shaped particles, whose volume follows from the vertices), the material and the
       meshfree approximation
   * - ``getPropertyNames``, ``setProperties``, ``setProperty``
     - particle properties, e.g., ``VCI order`` and the Newmark parameters ``newmark-beta beta``,
       ``newmark-beta gamma``; ``setProperties`` expects the values in the order of ``getPropertyNames``
   * - ``getNumberOfRequiredStateVars``, ``assignStateVars``, ``initializeYourself``
     - the host owns the state variable block; the particle maps its states into it
   * - ``assignMeshfreeKernelFunctions``
     - assign the nodes (kernel functions) and evaluate the trial functions :math:`N_A`, their (possibly smoothed)
       gradients :math:`\partial N_A/\partial Y_i` and the test functions :math:`T_A` (initially :math:`T_A = N_A`)
   * - ``vci_...``
     - optional: local contributions to the VCI correction of the test function gradients (see below)
   * - ``computePhysicsKernels``
     - internal force vector and its derivative with respect to the dof increment :math:`\Delta\boldsymbol{q}` since
       the last accepted state
   * - ``computeBodyLoad``, ``computeDistributedLoad``
     - external loads and their derivatives; the load type is an integer from ``getSupportedBodyLoadTypes`` /
       ``getSupportedDistributedLoadTypes`` (maps from the upper-case name, e.g., ``PRESSURE``, ``CWFCORRECTION``)
   * - ``acceptStateAndPosition``
     - after convergence: accept the state and move the particle (position, volume, geometry, VCI basis) to the new
       reference intermediate configuration
   * - ``getStateView``
     - a view (pointer and size) on a named state for output
   * - geometry queries
     - ``getParticleShape`` (Ensight Gold name, e.g., ``point``, ``quad4``, ``hexa8``), ``getNumberOfVertices``,
       ``getVertexCoordinates``, ``getVisualizationVertexCoordinates``, ``getCenterCoordinates``,
       ``getFaceCoordinates``, ``getEvaluationCoordinates`` and ``getNumberOfEvaluationPoints`` (the points at which
       the shape functions are evaluated, e.g., for checking the coverage by kernel supports),
       ``getInterpolationVector``

For explicit dynamics, ``updatePhysicsExplicit``, ``computePhysicsKernelsExplicit``, ``computeBodyLoadExplicit`` and
``computeDistributedLoadExplicit`` are provided, together with ``computeLumpedInertia`` and ``computeLumpedMomentum``.
The explicit methods default to doing nothing and the lumped methods throw, unless a particle overrides them.

**Dof layout.** Since the number of nodes of a particle is not fixed, no field-blocked layout (no
``dofIndicesPermutationPattern``) is used: the residual and the stiffness matrix are stored node-wise, in the order in
which the kernel functions were assigned, e.g., ``[node_1_displacement, node_1_temperature, node_2_displacement, ...]``,
with ``getNBaseDof()`` dofs per node and the fields listed by ``getFields()``. The stiffness matrix is column-major.
Contributions are accumulated (``+=``) into the arrays passed by the host.

**Loads.** Body loads are currently not available for the particles, although ``BODYFORCE`` is listed as a
supported body load type: the point particles and the SQCNI/SNNI/NSNI particles derived from them add nothing, and the
subdomain-integration particles throw ``std::runtime_error``. Distributed loads (``PRESSURE`` and, where available,
``CWFCORRECTION``) act on a face of a cell-shaped particle;
point particles have no faces and add nothing.

**State views.** Particles with a single material point forward ``getStateView`` to it. Cell-shaped particles
additionally provide ``vertex displacements`` and ``smoothing vertex displacements`` (``nDim * nVertices`` values
each, the displacements of the vertices of the deformed geometry and of the smoothing domain with respect to the
undeformed geometry), and the subdomain-integration particles forward all other names to the material point of
subdomain ``qp``.

The particle domain
-------------------

A cell-shaped particle is described by a :cpp:class:`Marmot::Meshfree::ParticleDomain`: a 4-node quadrilateral or an
8-node hexahedron (a ``MarmotLagrangeCell``, see :doc:`meshfreecells`), with faces numbered from 1 in the Abaqus
convention. It is kept in three versions:

- the **undeformed geometry**, with vertices :math:`\boldsymbol{X}_v` and centroid :math:`\boldsymbol{X}_c`, which
  also gives the undeformed volume;
- the **deformed geometry**, the particle in the reference intermediate configuration, used for the particle
  position, the face centers and boundary surface vectors :math:`\boldsymbol{N}\,dA_Y` of distributed loads, and the
  second moments :math:`\int (\boldsymbol{Y}-\boldsymbol{Y}_c)\otimes(\boldsymbol{Y}-\boldsymbol{Y}_c)\,dV` about the
  centroid;
- the **smoothing domain** :math:`\Omega_s`, over whose boundary the smoothed shape function gradients are
  integrated.

Both deformed versions are rebuilt from the undeformed geometry by a homogeneous deformation about the centroid plus
the center displacement :math:`\boldsymbol{u}_c`,

.. math::

   \boldsymbol{x}_v = \boldsymbol{X}_c + \boldsymbol{F}\,(\boldsymbol{X}_v - \boldsymbol{X}_c) + \boldsymbol{u}_c ,

with the total deformation gradient :math:`\boldsymbol{F}` of the particle for the deformed geometry, and a tensor
:math:`\boldsymbol{F}_s` chosen by the ``SmoothingDomainUpdateType`` for the smoothing domain:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - ``SmoothingDomainUpdateType``
     - :math:`\boldsymbol{F}_s`
     - name tag
   * - ``DeformationGradient``
     - :math:`\boldsymbol{F}`
     - ``SQCNI``
   * - ``None``
     - :math:`\boldsymbol{I}` (translated only)
     - ``SNNI``
   * - ``RotationOnly``
     - :math:`\boldsymbol{R}` of :math:`\boldsymbol{F} = \boldsymbol{R}\boldsymbol{U}` (from an SVD)
     - ``SQCNI_R``, ``R-SNNI``
   * - ``RotationAndPrincipalStretch``
     - :math:`\boldsymbol{R}\,\mathrm{diag}(\boldsymbol{R}^T\boldsymbol{F})`, i.e., the rotation and the diagonal
       entries of :math:`\boldsymbol{U}` in the global basis
     - ``SQCNI_RU``, ``RS-SNNI``

``ParticleDomain::uniformSubdivided`` splits the undeformed domain once in each direction (4 quadrilaterals or 8
hexahedra); every subdomain is again a ``ParticleDomain`` with the same update type. This is the basis of the subdomain
integration.

Integration schemes
-------------------

The particle type names combine the formulation with the integration scheme (see :doc:`meshfree` for the naming
pattern). The abbreviations denote what the code does:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Abbreviation
     - Meaning in Marmot
   * - (none), ``Point``
     - **Direct nodal integration** of a point particle (e.g., :cpp:class:`Marmot::Meshfree::GenericParticle`): shape
       functions and their gradients are evaluated at the particle center, which is the only integration point, with
       the particle volume as weight.
   * - SQCNI
     - **Stabilized quasi-conforming nodal integration**: one integration point at the particle center, with the
       gradients smoothed over the smoothing domain by the divergence theorem (one point per face, the face center),

       .. math::

          \frac{\partial N_A}{\partial Y_i} \approx \frac{1}{|\Omega_s|} \sum_f N_A(\boldsymbol{Y}_f)\, n_{f,i}\,dA_f ,

       where the smoothing domain deforms with the deformation gradient of the particle (``DeformationGradient``),
       so that neighboring smoothing domains remain (approximately) conforming. The suffixes ``_R`` and ``_RU``
       select the ``RotationOnly`` and ``RotationAndPrincipalStretch`` updates instead.
   * - SNNI
     - **Stabilized non-conforming nodal integration**: the same smoothed gradient, but the smoothing domain keeps its
       undeformed shape and is only translated with the particle (``None``). The prefixes ``R-`` and ``RS-`` select
       the ``RotationOnly`` and ``RotationAndPrincipalStretch`` updates.
   * - NSNI
     - **Naturally stabilized nodal integration**: in addition to the smoothed first derivatives, smoothed second
       derivatives of the shape functions are computed (by the divergence theorem applied to the first derivatives on
       the smoothing domain boundary). They enter a stabilization term built from a Taylor expansion of the stress
       gradient about the particle center, weighted with the second moments of the deformed particle geometry
       (``ParticleDomain::getGeometrySecondMoments``). Combined with SQCNI or SNNI smoothing, e.g.,
       ``SQCNIxNSNI``; see :doc:`displacementparticle`.
   * - SDI
     - **Subdomain integration** (:cpp:class:`Marmot::Meshfree::GenericSDIParticle`): the particle domain is
       subdivided into :math:`2^{n_\text{dim}}` subdomains, each an integration point with its own material point
       and its own smoothed gradients; the weak form is the sum over the subdomains. Combined with SQCNI or SNNI
       smoothing, e.g., ``SQCNIxSDI``.
   * - VCI
     - **Variationally consistent integration** (Chen, Hillman and Rüter, 2013): a Petrov-Galerkin correction of
       the test function gradients such that the integration by parts holds exactly under the numerical integration
       for a monomial basis of the chosen order; enabled through the property ``VCI order`` and computed by the host
       framework with the particles' ``vci_...`` contributions (see below). It is not part of the type names.
   * - CWF
     - **Consistent weak form** correction: the distributed load type ``CWFCORRECTION`` of the SQCNI-type particles
       adds, on a particle face at an essential boundary, the boundary term of the weak form that is otherwise
       dropped (the RK test functions do not vanish on the boundary). It is evaluated with the current stress,
       :math:`\boldsymbol{f}^{ext}_A \mathrel{-}= T_A(\boldsymbol{Y}_f)\,
       \boldsymbol{\tau}\,\Delta\boldsymbol{F}^{-T}\boldsymbol{N}\,dA_Y / J_Y`, with the Kirchhoff stress
       :math:`\boldsymbol{\tau}`, the deformation gradient increment :math:`\Delta\boldsymbol{F}` and the Jacobian
       :math:`J_Y` of the reference intermediate configuration, including its linearization. See
       :doc:`displacementparticle` and :doc:`gradientenhancedfinitestrainparticle`.

The complete list of registered names, with their shapes and dimensions, is given on the pages of the formulations,
:doc:`displacementparticle` and :doc:`gradientenhancedfinitestrainparticle`. The names are case-insensitive. Note the
two spellings of the NSNI names of the displacement formulation (``Displacement/SQCNIxNSNI/PlaneStrain/Quad``, with a
slash after the formulation) and of the gradient-enhanced formulation
(``GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad``).

GenericParticle: point particles and VCI
----------------------------------------

:cpp:class:`Marmot::Meshfree::GenericParticle` implements the physics-independent part of a point particle: a single
vertex and evaluation point at the center :math:`\boldsymbol{Y}_c` (shape ``point``), the evaluation of
:math:`N_A` and :math:`\partial N_A/\partial Y_i` at the center in ``assignMeshfreeKernelFunctions`` (which also resets
the test functions to :math:`T_A = N_A`), and the property ``VCI order``. Derived classes (e.g., the displacement
particles) keep the center, the volume :math:`V_Y` in the reference intermediate configuration and the VCI basis up to
date in ``acceptStateAndPosition``; the SQCNI particles derived from it replace the gradient by the smoothed one.

**Variationally consistent integration.** Meshfree shape functions integrated by nodal or smoothed quadrature do not,
in general, satisfy the integration by parts that the Galerkin exactness of the method relies on. VCI corrects the test
function gradients,

.. math::

   \frac{\partial T_A}{\partial Y_i} = \frac{\partial N_A}{\partial Y_i} + \chi_A \sum_{C} \eta_{AiC}\, P_C ,

with the complete monomial basis :math:`\boldsymbol{P}(\boldsymbol{Y})` of degree :math:`\le k` (the ``VCI order``),
which has :math:`n_C = \binom{k + n_\text{dim}}{n_\text{dim}}` entries, and with :math:`\chi_A \in \{0,1\}` indicating
whether the particle center lies in the support of kernel function :math:`A`. The coefficients :math:`\eta_{AiC}` are
chosen such that, under the particle quadrature, for every node :math:`A`

.. math::

   \int_\Omega \frac{\partial T_A}{\partial Y_i} P_C\,dV + \int_\Omega T_A \frac{\partial P_C}{\partial Y_i}\,dV
   = \int_{\partial\Omega} T_A P_C\, n_i\,dA .

The correction needs global integrals over all particles and the boundary, so it is split between the particles and
the host framework:

- each particle accumulates its contributions, for a point particle with the weight :math:`V_Y` and the basis at the
  center: ``vci_compute_TestGradient_P_Integral`` adds :math:`\partial_i T_A\,P_C\,V_Y`,
  ``vci_compute_Test_PGradient_Integral`` adds :math:`T_A\,\partial_i P_C\,V_Y`, and ``vci_compute_MMatrix`` adds the
  moment matrix :math:`M_{ACD} = \chi_A\,P_C\,P_D\,V_Y`; the particles at the boundary add the boundary integral
  (``vci_compute_Test_P_BoundaryIntegral``, which depends on the kinematics of the physics and is therefore
  implemented by the derived classes; the generic versions throw);
- the host assembles these per node :math:`A` and solves

  .. math::

     \sum_D M_{ACD}\,\eta_{AiD} = R_{AiC}, \qquad
     R_{AiC} = \int_{\partial\Omega} T_A P_C n_i\,dA - \int_\Omega \partial_i T_A\,P_C\,dV
     - \int_\Omega T_A\,\partial_i P_C\,dV ;

- each particle receives its :math:`\eta_{AiC}` in ``vci_assignTestFunctionCorrectionTerms``.

All VCI arrays are row-major, of shape :math:`n_\text{nodes} \times n_\text{dim} \times n_C` (index :math:`AiC`) or
:math:`n_\text{nodes} \times n_C \times n_C` (index :math:`ACD`). ``setVCIOrder`` computes :math:`n_C` for the given
dimension and evaluates the basis at the current center right away, since the VCI may run before the first accepted
increment updates it. In ``GenericParticle`` the correction is *added* to the current test function gradients, so it
is meant to be applied once after ``assignMeshfreeKernelFunctions``.

GenericSDIParticle: subdomain integration
-----------------------------------------

:cpp:class:`Marmot::Meshfree::GenericSDIParticle` is the physics-independent base of the subdomain-integration
particles. It owns the main :cpp:class:`Marmot::Meshfree::ParticleDomain` and its :math:`2^{n_\text{dim}}`
subdomains. For every subdomain :math:`s`, ``assignMeshfreeKernelFunctions`` evaluates :math:`N_A` at the subdomain
center, the smoothed gradient over the subdomain's smoothing domain (formula above), the test functions and the VCI
basis at the subdomain center. The derived class attaches a material point to every subdomain and integrates the weak
form as a sum over the subdomains (``computePhysicsKernelsOnSubdomains``).

The particle as a whole moves with its center. ``computePhysicsKernels`` first evaluates the shape functions of the
main domain and updates the center kinematics,

.. math::

   \Delta\boldsymbol{u}_c = \sum_B N_B\,\Delta\boldsymbol{q}_B, \qquad
   \Delta\boldsymbol{F}_c = \boldsymbol{I} + \sum_B \Delta\boldsymbol{q}_B \otimes \frac{\partial N_B}{\partial
   \boldsymbol{Y}} ,

(:math:`\boldsymbol{u}_c \mathrel{+}= \Delta\boldsymbol{u}_c`, which relies on the host restoring the state block of
the last accepted state before every iteration, as EdelweissMeshfree does), and then calls the physics on the
subdomains. ``acceptStateAndPosition`` sets
:math:`\boldsymbol{F}_c \leftarrow \Delta\boldsymbol{F}_c\,\boldsymbol{F}_c`, resets the increment to identity and moves
the main domain and all subdomains with :math:`\boldsymbol{F}_c` and :math:`\boldsymbol{u}_c`.

The VCI contributions are sums over the subdomains, with the subdomain volume :math:`V_s` (``getSubdomainVolume``) as
weight and the basis at the subdomain centers; :math:`\chi_A` is evaluated at the center of the main domain. Unlike in
``GenericParticle``, ``vci_assignTestFunctionCorrectionTerms`` sets
:math:`\partial_i T_A|_s = \partial_i N_A|_s + \chi_A\sum_C\eta_{AiC}P_C|_s`, so repeated calls do not accumulate, and
``setVCIOrder`` only sets :math:`n_C`; the basis is evaluated in ``assignMeshfreeKernelFunctions``, which must
therefore follow the setting of the order.

**State variable layout.** The block passed to ``assignStateVars`` is

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Offset
     - Size
     - Content
   * - 0
     - ``nDim``
     - center displacement :math:`\boldsymbol{u}_c`
   * - ``nDim``
     - ``nDim * nDim``
     - central deformation gradient :math:`\boldsymbol{F}_c` (column-major)
   * - ``nDim + nDim * nDim``
     - ``nDim * nDim``
     - its increment :math:`\Delta\boldsymbol{F}_c` (column-major)
   * - ``paddedStateVarSize(nStateVarsCenter)``
     - per subdomain ``paddedStateVarSize(n_s)``
     - the states of the subdomain material points, one block after the other

The center block has ``nStateVarsCenter = nDim + 2 nDim^2`` doubles (10 in 2D, 21 in 3D), padded to 16 and 24. Every
block is rounded up to a multiple of 8 doubles (64 bytes) by ``paddedStateVarSize``, so that the state of each
subdomain starts at the same alignment as the state of a stand-alone material point: Fastor may use aligned SIMD
stores on the tensor maps into it (without the padding, a block starting at an odd double caused a segmentation fault
with GCC 14.2 at ``-O3``). The offsets are relative to the start of the particle block, so the host framework must pass
a suitably aligned block for the absolute alignment to hold. ``getNumberOfRequiredStateVars`` includes the padding.

``GenericSDIParticle`` does not override the explicit-dynamics methods, so the SDI particles inherit their empty
defaults from ``MarmotParticle``.

The particle factory
--------------------

:cpp:class:`MarmotLibrary::MarmotParticleFactory` maps names to factory functions. Particle modules register their
types in ``*Registration.cpp`` files at static initialization (``registerParticle``); the host creates them with
``createParticle``, which throws ``std::invalid_argument`` for an unknown name. Names are converted to upper case at
registration and at creation, so they are case-insensitive. The factory function receives the particle number, the
vertex coordinates (the center coordinates for point particles), their number, the volume (zero for cell-shaped
particles, which compute it from the vertices), the material name and properties, and the meshfree approximation,
which must outlive the particle.

Implementation
--------------

.. doxygenclass:: Marmot::Meshfree::MarmotParticle
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::ParticleDomain
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::GenericParticle
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::GenericSDIParticle
   :allow-dot-graphs:

.. doxygenclass:: MarmotLibrary::MarmotParticleFactory
   :allow-dot-graphs:
