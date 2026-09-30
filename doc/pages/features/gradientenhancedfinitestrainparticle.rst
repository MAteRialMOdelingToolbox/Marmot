Gradient-Enhanced Finite-Strain Particles
=========================================

Theory
------
The gradient-enhanced finite-strain particles are the RKPM particles of the implicit-gradient finite-strain
formulation. Each particle couples the displacement field :math:`\boldsymbol u` to a single scalar nonlocal field
:math:`\bar N` (nodal fields ``displacement`` and ``nonlocal damage``, dofs ordered node by node as
:math:`\{u_1,\dots,u_{n_\mathrm{dim}},\bar N\}`), owns one :doc:`gradientenhancedfinitestrainmaterialpoint` and is
integrated at that single point (nodal integration). The material is any material of the interface
``MarmotMaterialGradientEnhancedFiniteStrain``, e.g. :doc:`gradientenhancedcompressibleneohookedamage` or
:doc:`gradientenhancedfinitestraindruckerprager`. The shape functions come from the reproducing kernel approximation
(:doc:`meshfreeapproximation`); the particle interface, the integration schemes and the corrections are introduced in
:doc:`meshfreeparticles`. The weak forms are those of :doc:`gradientenhancedfinitestraindisplacementelement`, in
increment form.

There are three particle classes:

- ``GradientEnhancedFiniteStrainParticle``: a point particle, direct nodal integration with the point values of the
  shape function gradients;
- ``GradientEnhancedFiniteStrainParticleSQCNI``: a quadrilateral or hexahedral smoothing domain with smoothed gradients
  (SQCNI / SNNI) and faces for distributed loads;
- ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI``: the SQCNI particle plus the NSNI stabilization.

Kinematics
^^^^^^^^^^
The dofs are the increments :math:`\Delta\boldsymbol q` of the current step. The shape functions :math:`N_B` and their
gradients :math:`\partial N_B/\partial\boldsymbol Y` are evaluated in the intermediate reference configuration
:math:`\boldsymbol Y` of the last accepted state, when the host assigns the kernel functions. They give

.. math::

   \Delta\boldsymbol u = N_B\,\Delta\boldsymbol q^U_B, \qquad
   \Delta F_{ij} = \delta_{ij} + \Delta q^U_{Bi}\,\frac{\partial N_B}{\partial Y_j}, \qquad
   \Delta\bar N = N_B\,\Delta q^N_B, \qquad
   \frac{\partial\Delta\bar N}{\partial Y_j} = \Delta q^N_B\,\frac{\partial N_B}{\partial Y_j},

and the material point is evaluated at :math:`\boldsymbol F = \Delta\boldsymbol F\,\boldsymbol F_n`. Spatial and
undeformed gradients follow from :math:`\partial(\cdot)/\partial x_i = \partial(\cdot)/\partial Y_j\,\Delta F^{-1}_{ji}`
and :math:`\partial(\cdot)/\partial X_i = \partial(\cdot)/\partial Y_j\,F_{n,ji}`. On acceptance, the material point
state is updated, the intermediate reference configuration moves (center, volume :math:`V_Y = V_0\det\boldsymbol F_n`,
smoothing domain); the shape functions are re-evaluated when the host assigns the kernel functions again.

Weak forms
^^^^^^^^^^
With the test functions :math:`T_A` (equal to :math:`N_A`, except that VCI corrects their gradients), the undeformed
particle volume :math:`V_0` and :math:`c = R^2` from the material, the residuals per node :math:`A` are

.. math::

   r^U_{Aj} &= \frac{\partial T_A}{\partial x_i}\,\tau_{ij}\,V_0 + \rho_0\,a_j\,T_A\,V_0 ,\\
   r^N_A &= \Bigl( T_A\,\Delta\bar N + c\,\frac{\partial T_A}{\partial X_i}\,\frac{\partial\Delta\bar N}{\partial X_i}
   - T_A\,\Delta L \Bigr) V_0 ,

i.e. the Helmholtz equation :math:`\bar N - c\,\nabla_X^2\bar N = L` is solved in the undeformed configuration for the
increment of the nonlocal field, with the change :math:`\Delta L` of the local driving force as its source. The
acceleration :math:`\boldsymbol a` follows from the Newmark-beta update of :math:`\Delta\boldsymbol u` with the particle
properties ``newmark-beta beta`` and ``newmark-beta gamma``; :math:`\beta = 0` switches the inertia off
(:math:`\boldsymbol a = \boldsymbol 0`).

Tangent
^^^^^^^
With the tangents of the material point with respect to :math:`\Delta\boldsymbol F` and :math:`\bar N`,

.. math::

   K^{UU}_{AjBk} &= \Bigl( \frac{\partial T_A}{\partial x_i}\,\frac{\partial\tau_{ij}}{\partial\Delta F_{kL}}\,
   \frac{\partial N_B}{\partial Y_L} - \frac{\partial T_A}{\partial x_k}\,\tau_{ij}\,\frac{\partial N_B}{\partial x_i}
   + \rho_0\,\frac{\partial a_j}{\partial\Delta u_k}\,T_A N_B \Bigr) V_0 ,\\
   K^{UN}_{AjB} &= \frac{\partial T_A}{\partial x_i}\,\frac{\partial\tau_{ij}}{\partial\bar N}\,N_B\,V_0 ,\\
   K^{NU}_{ABk} &= -T_A\,\frac{\partial L}{\partial\Delta F_{kL}}\,\frac{\partial N_B}{\partial Y_L}\,V_0 ,\\
   K^{NN}_{AB} &= \Bigl( T_A N_B + c\,\frac{\partial T_A}{\partial X_i}\,\frac{\partial N_B}{\partial X_i} \Bigr) V_0 .

.. note::
   Unlike :doc:`gradientenhancedfinitestraincell` and :doc:`gradientenhancedfinitestraindisplacementelement`, the
   particles' :math:`K^{NN}` does not contain the term :math:`-T_A N_B\,\partial L/\partial\bar N\,V_0`. The tangent is
   therefore consistent for materials whose driving force does not depend on :math:`\bar N`, which is the case for
   :doc:`gradientenhancedcompressibleneohookedamage` and :doc:`gradientenhancedfinitestraindruckerprager`
   (both report :math:`\partial L/\partial\bar N = 0`).

Integration schemes
^^^^^^^^^^^^^^^^^^^

**Point particle.** The gradients :math:`\partial N_B/\partial\boldsymbol Y` are the point values of the meshfree
approximation at the particle center. The point particle has no faces, so it cannot carry distributed loads; its
``computeDistributedLoad`` adds nothing, although ``PRESSURE`` and ``CWFCORRECTION`` are listed as supported load
types.

.. note::
   ``BODYFORCE`` is listed as a supported body load of all gradient-enhanced particles, but ``computeBodyLoad`` is
   empty (and not overridden by the SQCNI / NSNI particles), so a body force currently has no effect.

**SQCNI / SNNI.** The particle is a quadrilateral (2D) or hexahedral (3D) smoothing domain :math:`\Omega_Y` with the
material point at its centroid. As in stabilized conforming nodal integration (Chen et al., 2001), the gradients are
replaced by their average over the smoothing domain, turned into a boundary integral and evaluated with one point per
face :math:`f` (the face center :math:`\boldsymbol Y_f`),

.. math::

   \frac{\partial N_B}{\partial Y_i} \approx \frac{1}{V_{\Omega_Y}} \sum_f N_B(\boldsymbol Y_f)\,n_i\,dA_f ,

while :math:`N_B` is the point value at the center. The smoothed gradients are used for both fields, i.e. also for
:math:`\nabla\Delta\bar N` and the Laplacian term of the nonlocal balance. The smoothing domain is stored by its vertex
displacements (state ``vertex displacements``) and is updated at each accepted increment from :math:`\boldsymbol F_n` of
the material point, according to the smoothing domain update type:

.. list-table::
   :header-rows: 1

   * - Update type
     - Mapping of the undeformed domain about its center
     - Name token
   * - ``DeformationGradient``
     - :math:`\boldsymbol F_n`; the domain conforms to the deformed particle
     - ``SQCNI``
   * - ``None``
     - identity; the domain is only translated with the particle
     - ``SNNI``
   * - ``RotationOnly``
     - the rotation :math:`\boldsymbol R` of the polar decomposition of :math:`\boldsymbol F_n`
     - ``SQCNI_R``
   * - ``RotationAndPrincipalStretch``
     - :math:`\boldsymbol R\,\mathrm{diag}(\boldsymbol R^\mathsf{T}\boldsymbol F_n)`
     - ``SQCNI_RU``

**NSNI.** Nodal integration evaluates the integrand at one point only, which leaves spurious zero-energy modes. The
naturally stabilized nodal integration (NSNI) expands the integrand to first order about the particle center and keeps
the second-order term of the integral, with the second moments :math:`\boldsymbol M^Y` of the particle domain about its
centroid. The second derivatives of the shape functions are smoothed like the first ones, from the gradients at the
face centers (symmetrized),

.. math::

   \frac{\partial^2 N_B}{\partial Y_I\,\partial Y_J} \approx \frac{1}{2V_{\Omega_Y}} \sum_f \Bigl(
   \frac{\partial N_B}{\partial Y_I}(\boldsymbol Y_f)\,n_J + \frac{\partial N_B}{\partial Y_J}(\boldsymbol Y_f)\,n_I
   \Bigr) dA_f ,

the Kirchhoff stress gradient is linearized as

.. math::

   \frac{\partial\tau_{ij}}{\partial Y_K} \approx \frac{\partial\tau_{ij}}{\partial\Delta F_{mM}}\,
   \Delta q^U_{Bm}\,\frac{\partial^2 N_B}{\partial Y_M\,\partial Y_K},

and the momentum residual receives

.. math::

   r^{U,\mathrm{stab}}_{Aj} = \Delta F^{-1}_{mi}\,\frac{\partial^2 N_A}{\partial Y_m\,\partial Y_L}\,
   \frac{M^Y_{LK}}{J_Y}\,\frac{\partial\tau_{ij}}{\partial Y_K}, \qquad
   \boldsymbol M^Y = J_Y\,\boldsymbol F_n\,\boldsymbol M^0\,\boldsymbol F_n^\mathsf{T},\quad J_Y = \det\boldsymbol F_n,

with the second moments :math:`\boldsymbol M^0` of the undeformed cell (computed at construction; :math:`\boldsymbol M^Y`
is initialized to :math:`\boldsymbol M^0` and pushed forward at each accepted increment). The stabilization acts on the
momentum balance only; the stress gradient does not contain the contribution
:math:`\partial\boldsymbol\tau/\partial\bar N\,\nabla_Y\Delta\bar N` of the nonlocal field, and the second derivatives
are not corrected by VCI.

.. note::
   The tangent of the NSNI stabilization is **approximate by construction**. It differentiates the second derivatives
   of the displacement increment and :math:`\Delta\boldsymbol F^{-1}`, but not :math:`\partial\boldsymbol\tau/\partial
   \Delta\boldsymbol F` itself: :math:`\partial^2\boldsymbol\tau/\partial\boldsymbol F^2` is not exposed by the
   material interface, and the dependence of :math:`\partial\boldsymbol\tau/\partial\Delta\boldsymbol F` on
   :math:`\bar N` (the damage) is omitted, so the U-N block has no stabilization contribution at all. The rows of the
   nonlocal balance are exact. The module test bounds the error of the U-U block (below :math:`5\cdot10^{-3}`) and does
   not check the U-N block for the NSNI particles, whose error in a damaged state is not small.

Corrections
^^^^^^^^^^^

**VCI.** All particles implement the hooks of variationally consistent integration (Chen, Hillman, Rüter, 2013): the
host assembles the integration constraints from the particles' boundary and domain terms and returns correction
coefficients :math:`\eta_{AiC}`, with which the test function gradients are corrected,
:math:`\partial T_A/\partial Y_i \mathrel{+}= \eta_{AiC}\,P_C(\boldsymbol Y)`, with a monomial basis :math:`P_C` of
order ``VCI order`` at the particle center. The corrected test gradients enter both weak forms. For the SQCNI
particles, the boundary term is evaluated at the face centers with the face surface vectors of the intermediate
configuration.

**CWF.** The distributed load ``CWFCORRECTION`` (SQCNI particles, per face) adds the boundary term of the momentum weak
form on a face, which the smoothed gradients otherwise leave unbalanced. With the surface vector
:math:`\boldsymbol N\,dA_Y` and the center :math:`\boldsymbol Y_f` of the face in the intermediate configuration,

.. math::

   \boldsymbol t = \boldsymbol\tau\,\Delta\boldsymbol F^{-\mathsf T}\,\frac{\boldsymbol N\,dA_Y}{J_Y}, \qquad
   f_{Ai} \mathrel{-}= T_A(\boldsymbol Y_f)\,t_i ,

where the division by :math:`J_Y` accounts for the internal force being integrated over the undeformed volume
:math:`V_0 = V_Y/J_Y`. The tangent contains the geometric part from :math:`\Delta\boldsymbol F^{-\mathsf T}`, the
material part from :math:`\partial\boldsymbol\tau/\partial\Delta\boldsymbol F` and the coupling
:math:`\partial\boldsymbol\tau/\partial\bar N` to the nonlocal dofs. For a homogeneous deformation, the corrections on
all faces cancel the displacement residual (exactly so after deformation only for the ``DeformationGradient``
update, whose smoothing domain coincides with the particle).

**Pressure.** ``PRESSURE`` (SQCNI particles, per face, ``load[0]`` = :math:`p`) is a follower load,
:math:`\boldsymbol f = \Delta J\,\Delta\boldsymbol F^{-\mathsf T}\,p\,\boldsymbol N\,dA_Y`,
:math:`f_{Ai} \mathrel{-}= T_A(\boldsymbol Y_f)\,f_i`, with its load stiffness. For ``DeformationGradient`` the face
geometry is that of the current smoothing domain; for the other update types the undeformed surface vector is mapped by
Nanson's formula with :math:`\boldsymbol F_n`, and the undeformed face center is translated by the particle
displacement. Neither load acts on the nonlocal field.

Registered types
----------------

The particles are registered in the ``MarmotParticleFactory``:

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension / shape
     - Integration
   * - ``GradientEnhancedFiniteStrain/PlaneStrain/Point``
     - ``GradientEnhancedFiniteStrainParticle<2>``
     - plane strain, point
     - direct nodal integration
   * - ``GradientEnhancedFiniteStrain/3D/Point``
     - ``GradientEnhancedFiniteStrainParticle<3>``
     - 3D, point
     - direct nodal integration
   * - ``GradientEnhancedFiniteStrainSQCNI/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<2,4>``
     - plane strain, Quad
     - SQCNI (``DeformationGradient``)
   * - ``GradientEnhancedFiniteStrainSNNI/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<2,4>``
     - plane strain, Quad
     - SNNI (``None``)
   * - ``GradientEnhancedFiniteStrainSQCNI_R/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<2,4>``
     - plane strain, Quad
     - SQCNI (``RotationOnly``)
   * - ``GradientEnhancedFiniteStrainSQCNI_RU/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<2,4>``
     - plane strain, Quad
     - SQCNI (``RotationAndPrincipalStretch``)
   * - ``GradientEnhancedFiniteStrainSQCNI/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<3,8>``
     - 3D, Hexa
     - SQCNI (``DeformationGradient``)
   * - ``GradientEnhancedFiniteStrainSNNI/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<3,8>``
     - 3D, Hexa
     - SNNI (``None``)
   * - ``GradientEnhancedFiniteStrainSQCNI_R/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<3,8>``
     - 3D, Hexa
     - SQCNI (``RotationOnly``)
   * - ``GradientEnhancedFiniteStrainSQCNI_RU/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNI<3,8>``
     - 3D, Hexa
     - SQCNI (``RotationAndPrincipalStretch``)
   * - ``GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<2,4>``
     - plane strain, Quad
     - SQCNI + NSNI (``DeformationGradient``)
   * - ``GradientEnhancedFiniteStrainSNNIxNSNI/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<2,4>``
     - plane strain, Quad
     - SNNI + NSNI (``None``)
   * - ``GradientEnhancedFiniteStrainSQCNI_RxNSNI/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<2,4>``
     - plane strain, Quad
     - SQCNI + NSNI (``RotationOnly``)
   * - ``GradientEnhancedFiniteStrainSQCNI_RUxNSNI/PlaneStrain/Quad``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<2,4>``
     - plane strain, Quad
     - SQCNI + NSNI (``RotationAndPrincipalStretch``)
   * - ``GradientEnhancedFiniteStrainSQCNIxNSNI/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<3,8>``
     - 3D, Hexa
     - SQCNI + NSNI (``DeformationGradient``)
   * - ``GradientEnhancedFiniteStrainSNNIxNSNI/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<3,8>``
     - 3D, Hexa
     - SNNI + NSNI (``None``)
   * - ``GradientEnhancedFiniteStrainSQCNI_RxNSNI/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<3,8>``
     - 3D, Hexa
     - SQCNI + NSNI (``RotationOnly``)
   * - ``GradientEnhancedFiniteStrainSQCNI_RUxNSNI/3D/Hexa``
     - ``GradientEnhancedFiniteStrainParticleSQCNIxNSNI<3,8>``
     - 3D, Hexa
     - SQCNI + NSNI (``RotationAndPrincipalStretch``)

The point particles are constructed from the center coordinates and the volume; the Quad / Hexa particles from their
vertex coordinates, from which the centroid and the volume are computed (the volume argument is not used).

Properties and state variables
------------------------------

Particle properties (``setProperties`` in this order, or ``setProperty`` by name):

.. list-table::
   :header-rows: 1

   * - Name
     - Meaning
     - Default
   * - ``newmark-beta beta``
     - Newmark-beta parameter :math:`\beta`
     - 0 (no inertia)
   * - ``newmark-beta gamma``
     - Newmark-beta parameter :math:`\gamma`
     - 0
   * - ``VCI order``
     - order of the VCI monomial basis
     - 0

The state variables are those of the :doc:`gradientenhancedfinitestrainmaterialpoint` (including the material). The
Quad / Hexa particles append the ``vertex displacements`` of the smoothing domain (:math:`n_\mathrm{dim}\times
n_\mathrm{vertices}` values); for the point particle, ``vertex displacements`` is the displacement of the material
point. The material point's density is available from the start, since it is queried at initialization.

Implementation
--------------

.. doxygenclass:: Marmot::Meshfree::GradientEnhancedFiniteStrainParticle
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::GradientEnhancedFiniteStrainParticleSQCNI
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::GradientEnhancedFiniteStrainParticleSQCNIxNSNI
   :allow-dot-graphs:
