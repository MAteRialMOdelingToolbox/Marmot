Displacement Particle
=====================

Theory
------

The displacement particles are the RKPM particles of the finite-strain displacement formulation. Each particle
integrates the weak form of the momentum balance over its domain with one of several nodal integration schemes; the
material state is carried by one :doc:`displacementmaterialpoint` per integration point, which consumes any material of
the ``MarmotMaterialFiniteStrain`` interface. The shape functions are those of the reproducing kernel approximation
(:doc:`meshfreeapproximation`); the particle interface, the particle domain and the integration schemes SQCNI, SNNI,
NSNI, SDI and the corrections VCI and CWF are introduced in :doc:`meshfreeparticles`. The MPM counterpart of this
formulation is the :doc:`displacementcell` (see :doc:`meshfreecells`).

Kinematics and weak form
^^^^^^^^^^^^^^^^^^^^^^^^

The configurations are those of the material point: the undeformed configuration :math:`\boldsymbol{X}`, the
intermediate configuration :math:`\boldsymbol{Y}` (the last accepted state) and the current configuration
:math:`\boldsymbol{x}`, with :math:`\boldsymbol{F} = \Delta\boldsymbol{F}\,\boldsymbol{F}_n`. The shape functions and
their gradients are evaluated at the particle positions in the intermediate configuration when the host assigns the
kernel functions (``assignMeshfreeKernelFunctions``), i.e., gradients are taken with respect to :math:`\boldsymbol{Y}`,
and the dofs :math:`\Delta q_{Bk}` are the increments of
the current step at the nodes (kernel functions) of the particle. With the trial functions :math:`N_B` and the test
functions :math:`T_A` (identical, unless VCI corrects the test gradients),

.. math::

   \Delta\boldsymbol{u} = N_B\,\Delta\boldsymbol{q}_B, \qquad
   \Delta F_{iJ} = \delta_{iJ} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial Y_J}, \qquad
   \frac{\partial(\bullet)}{\partial x_i} = \Delta F^{-1}_{Ji}\,\frac{\partial(\bullet)}{\partial Y_J}.

The residual of an integration point with the undeformed volume :math:`V_0` and the undeformed density :math:`\rho_0`
is written with the Kirchhoff stress :math:`\boldsymbol{\tau}` (i.e., the weak form in the undeformed configuration),

.. math::

   r_{Aj} = \frac{\partial T_A}{\partial x_i}\,\tau_{ij}\,V_0 + \rho_0\,a_j\,T_A\,V_0 ,

and its tangent consists of the material part, the geometric stiffness and the inertia,

.. math::

   \frac{\partial r_{Aj}}{\partial \Delta q_{Bk}} =
     \left( \frac{\partial T_A}{\partial x_i}\,\frac{\partial \tau_{ij}}{\partial \Delta F_{kL}}\,
     \frac{\partial N_B}{\partial Y_L}
     - \frac{\partial T_A}{\partial x_k}\,\tau_{ij}\,\frac{\partial N_B}{\partial x_i} \right) V_0
     + \rho_0\,\frac{\partial a_j}{\partial \Delta u_k}\,T_A\,N_B\,V_0 .

The acceleration :math:`\boldsymbol{a}` follows from a Newmark-beta integration of the displacement increment
:math:`\Delta\boldsymbol{u}` of the integration point, with the properties ``newmark-beta beta`` and
``newmark-beta gamma``; for :math:`\beta = 0` (the default) the scheme returns :math:`\boldsymbol{a} = \boldsymbol{0}`,
i.e., a quasi-static analysis. The particles also provide a lumped mass :math:`\rho_0\,T_A\,V_0` and a lumped momentum
:math:`\rho_0\,T_A\,V_0\,\boldsymbol{v}` for explicit schemes.

Integration schemes
^^^^^^^^^^^^^^^^^^^

**Direct nodal integration** (``DisplacementParticle``). The particle is a point: shape functions and their direct
gradients are evaluated at the particle center, with the particle volume as the weight.

**SQCNI / SNNI** (``DisplacementParticleSQCNI``). The particle has a geometry (a quadrilateral or a hexahedron) and a
smoothing domain :math:`\Omega_s` in the intermediate configuration. The shape functions are evaluated at the particle
center, their gradients are smoothed over the smoothing domain by a one-point rule at its face centers
:math:`\boldsymbol{Y}_f`,

.. math::

   \frac{\partial N_B}{\partial Y_J} \approx \frac{1}{V_s} \sum_f N_B(\boldsymbol{Y}_f)\,(N_J\,dA)_f ,

and used in the residual and the tangent above. The update of the smoothing domain decides the scheme: with the
deformation gradient the smoothing domain conforms to the deformed particle (SQCNI); without an update it is only
translated with the particle (SNNI); the variants ``_R`` and ``_RU`` update it with the rotation, or with the rotation
and the principal stretches, of :math:`\boldsymbol{F}_n`.

**NSNI** (``DisplacementParticleSQCNIxNSNI``). In addition to the smoothed gradients of SQCNI/SNNI, smoothed second
derivatives are computed from the gradients at the face centers,
:math:`N_{B,IJ} = \mathrm{sym}\left[ V_s^{-1} \sum_f \partial N_B/\partial Y_I(\boldsymbol{Y}_f)\,(N_J\,dA)_f \right]`.
The residual is augmented by the naturally stabilized nodal integration term, the second-order term of a Taylor
expansion of the test function gradient and of the stress about the particle center,

.. math::

   r^{\mathrm{stab}}_{Aj} = \Delta F^{-1}_{Ii}\,N_{A,IJ}\,\frac{M_{JK}}{J_Y}\,\frac{\partial \tau_{ij}}{\partial Y_K},
   \qquad
   \frac{\partial \tau_{ij}}{\partial Y_K} \approx \frac{\partial \tau_{ij}}{\partial \Delta F_{mM}}\,
   \Delta q_{Bm}\,N_{B,MK},

with the second moments :math:`M_{JK} = \int_{\Omega_Y} (Y_J - Y_{c,J})(Y_K - Y_{c,K})\,dV` of the particle geometry in
the intermediate configuration and :math:`J_Y = \det\boldsymbol{F}_n` (so that :math:`M/J_Y` refers to the undeformed
volume, like :math:`V_0` in the residual). The second moments are initialized at construction and updated on
acceptance of each increment. The explicit variant of the particle additionally adds the corresponding gradient term
:math:`\rho_0\,(\partial v_i/\partial Y_J)\,M_{JK}\,(\partial T_A/\partial Y_K)/J_Y` to the lumped momentum.

.. note::

   The tangent of the stabilization differentiates the second gradient of the increment and the inverse of
   :math:`\Delta\boldsymbol{F}`, but not :math:`\partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F}` itself: the
   term with :math:`\partial^2\boldsymbol{\tau}/\partial\Delta\boldsymbol{F}^2` is omitted, since the material interface
   does not provide it. The NSNI tangent is therefore approximate by construction; the module test measures a relative
   error of about :math:`10^{-3}` against a numerical tangent. The tangents of all other particles are exact.

**SDI** (``DisplacementParticleSQCNIxSDI``, built on ``GenericSDIParticle``). The particle geometry is subdivided
uniformly into :math:`2^{n_\mathrm{dim}}` subdomains. Each subdomain :math:`s` owns a material point at its center with
its undeformed volume :math:`V_{0,s}`, the shape functions at its center and the gradients smoothed over its own
smoothing domain (with the same update types as above). The residual and the tangent are the sums of the
subdomain contributions,

.. math::

   r_{Aj} = \sum_s \left( \frac{\partial T^s_A}{\partial x_i}\,\tau^s_{ij} + \rho_0\,a^s_j\,T^s_A \right) V_{0,s},

with exact tangents. A central incremental deformation gradient, computed from the gradients smoothed over the whole
particle, moves the particle geometry and the subdomains on acceptance and enters the pressure load.

Boundary terms
^^^^^^^^^^^^^^

The particles with a geometry assemble loads on their faces. With the boundary vector :math:`\boldsymbol{N}\,dA_Y` of
the face in the intermediate configuration and the test functions :math:`T_A` at the face evaluation point
:math:`\boldsymbol{Y}_N` (the face center of the smoothing domain for SQCNI, of the particle geometry for the other
update types), a face force :math:`\boldsymbol{f}` is assembled as :math:`r_{Aj} \mathrel{-}= T_A(\boldsymbol{Y}_N)\,f_j`,
with :math:`\partial f_j/\partial\Delta F_{kL}\,\partial N_B/\partial Y_L` in the tangent:

- ``PRESSURE`` (follower load, load value :math:`p`): Nanson's formula from :math:`\boldsymbol{Y}` to
  :math:`\boldsymbol{x}`,

  .. math::

     \boldsymbol{f} = \Delta J\,\Delta\boldsymbol{F}^{-\mathsf T}\,p\,\boldsymbol{N}\,dA_Y ;

- ``CWFCORRECTION`` (SQCNI, SNNI and NSNI particles; no load value): the consistent weak form correction, the
  Cauchy traction of the particle's own stress on the face,

  .. math::

     \boldsymbol{f} = \boldsymbol{\tau}\,\Delta\boldsymbol{F}^{-\mathsf T}\,\boldsymbol{N}\,\frac{dA_Y}{J_Y}
       = \boldsymbol{\sigma}\,\boldsymbol{n}\,da ,

  since the internal force integrates :math:`\boldsymbol{\tau}` over the undeformed volume
  :math:`V_0 = V_Y/J_Y`. Its tangent contains the material part (from
  :math:`\partial\boldsymbol{\tau}/\partial\Delta\boldsymbol{F}`) and the geometric part (from
  :math:`\Delta\boldsymbol{F}^{-\mathsf T}`). Which faces receive the correction is decided by the host.

For VCI, the particles provide the boundary term :math:`T_A\,P_C\,(N\,dA_Y)_i` in the intermediate configuration; the
point particle transforms the given undeformed boundary vector by
:math:`\boldsymbol{N}\,dA_Y = J_n\,\boldsymbol{F}_n^{-\mathsf T}\,\boldsymbol{N}\,dA_0`.

Registered particles
--------------------

.. list-table::
   :header-rows: 1

   * - Name
     - Class
     - Dimension, shape
     - Integration / smoothing domain update
   * - ``Displacement/PlaneStrain/Point``
     - ``DisplacementParticle<2>``
     - 2D, point
     - direct nodal integration
   * - ``DisplacementSQCNI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNI<2, 4>``
     - 2D, quadrilateral
     - SQCNI (deformation gradient)
   * - ``DisplacementSQCNI/3D/Hexa``
     - ``DisplacementParticleSQCNI<3, 8>``
     - 3D, hexahedron
     - SQCNI (deformation gradient)
   * - ``DisplacementSNNI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNI<2, 4>``
     - 2D, quadrilateral
     - SNNI (none)
   * - ``DisplacementSQCNI_R/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNI<2, 4>``
     - 2D, quadrilateral
     - rotation only
   * - ``DisplacementSQCNI_RU/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNI<2, 4>``
     - 2D, quadrilateral
     - rotation and principal stretch
   * - ``Displacement/SQCNIxNSNI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNIxNSNI<2, 4>``
     - 2D, quadrilateral
     - SQCNI + NSNI (deformation gradient)
   * - ``Displacement/SQCNIxNSNI/3D/Hexa``
     - ``DisplacementParticleSQCNIxNSNI<3, 8>``
     - 3D, hexahedron
     - SQCNI + NSNI (deformation gradient)
   * - ``Displacement/SNNIxNSNI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNIxNSNI<2, 4>``
     - 2D, quadrilateral
     - SNNI + NSNI (none)
   * - ``Displacement/SNNIxNSNI/3D/Hexa``
     - ``DisplacementParticleSQCNIxNSNI<3, 8>``
     - 3D, hexahedron
     - SNNI + NSNI (none)
   * - ``Displacement/R-SNNIxNSNI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNIxNSNI<2, 4>``
     - 2D, quadrilateral
     - NSNI (rotation only)
   * - ``Displacement/R-SNNIxNSNI/3D/Hexa``
     - ``DisplacementParticleSQCNIxNSNI<3, 8>``
     - 3D, hexahedron
     - NSNI (rotation only)
   * - ``Displacement/RS-SNNIxNSNI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNIxNSNI<2, 4>``
     - 2D, quadrilateral
     - NSNI (rotation and principal stretch)
   * - ``Displacement/RS-SNNIxNSNI/3D/Hexa``
     - ``DisplacementParticleSQCNIxNSNI<3, 8>``
     - 3D, hexahedron
     - NSNI (rotation and principal stretch)
   * - ``DisplacementSQCNIxSDI/PlaneStrain/Quad``
     - ``DisplacementParticleSQCNIxSDI<2, 4>``
     - 2D, quadrilateral
     - SDI, SQCNI on the subdomains (deformation gradient)
   * - ``DisplacementSQCNIxSDI/3D/Hexa``
     - ``DisplacementParticleSQCNIxSDI<3, 8>``
     - 3D, hexahedron
     - SDI, SQCNI on the subdomains (deformation gradient)
   * - ``DisplacementSNNIxSDI/3D/Hexa``
     - ``DisplacementParticleSQCNIxSDI<3, 8>``
     - 3D, hexahedron
     - SDI, SNNI on the subdomains (none)
   * - ``DisplacementR-SNNIxSDI/3D/Hexa``
     - ``DisplacementParticleSQCNIxSDI<3, 8>``
     - 3D, hexahedron
     - SDI (rotation only)
   * - ``DisplacementRS-SNNIxSDI/3D/Hexa``
     - ``DisplacementParticleSQCNIxSDI<3, 8>``
     - 3D, hexahedron
     - SDI (rotation and principal stretch)

.. note::

   The NSNI particles are registered with an additional slash between the formulation and the integration scheme
   (``Displacement/SQCNIxNSNI/...``), unlike all other displacement particles (``DisplacementSQCNI/...``,
   ``DisplacementSQCNIxSDI/...``). The names are listed here exactly as registered.

Supported loads
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Particle
     - ``PRESSURE``
     - ``CWFCORRECTION``
     - ``BODYFORCE``
   * - ``DisplacementParticle`` (point)
     - -- (no faces)
     - -- (no faces)
     - yes
   * - ``DisplacementParticleSQCNI``, ``DisplacementParticleSQCNIxNSNI``
     - yes (also explicit)
     - yes
     - yes
   * - ``DisplacementParticleSQCNIxSDI``
     - yes
     - --
     - yes (over the subdomains)

The body force :math:`\boldsymbol b` is a dead load per unit undeformed volume,
:math:`P_{Ai} \mathrel{-}= T_A\,b_i\,V_0` (summed over the subdomains for SDI), with the host's sign convention for
external loads; it has no tangent.

Properties
----------

The properties are set in this order:

.. list-table::
   :header-rows: 1

   * - Name
     - Meaning
   * - ``VCI order``
     - order of the monomial basis :math:`P_C` used by the variationally consistent integration (default 0)
   * - ``newmark-beta beta``
     - Newmark parameter :math:`\beta` (0: no inertia)
   * - ``newmark-beta gamma``
     - Newmark parameter :math:`\gamma`

State variables
---------------

The point, SQCNI and NSNI particles carry the states of their :doc:`displacementmaterialpoint` (including the material
state). The particles with a geometry additionally provide the states ``vertex displacements`` and
``smoothing vertex displacements`` (``nDim`` x ``nVertices`` values) of the particle geometry and of the smoothing
domain. The SDI particle stores a central block (center displacement, central deformation gradient and its increment,
padded to a multiple of 8 values), followed by the states of the material points of the subdomains, each padded to a
multiple of 8 values; ``getStateView`` addresses a subdomain by the evaluation-point index ``qp``.

Implementation
--------------

.. doxygenclass:: Marmot::Meshfree::DisplacementParticle
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::DisplacementParticleSQCNI
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::DisplacementParticleSQCNIxNSNI
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Meshfree::DisplacementParticleSQCNIxSDI
   :allow-dot-graphs:
