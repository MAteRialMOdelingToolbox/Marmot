Axisymmetric Particles
======================

The axisymmetric particles are the RKPM particles of a body of revolution under a load of revolution. They live in
the meridian half plane :math:`(r, z)`, :math:`r = Y_1 \ge 0`, :math:`z = Y_2`, and every integral over a particle is an
integral over the *ring* the particle sweeps about the axis. They are the axisymmetric twins of the plane-strain
particles,

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Axisymmetric particle
     - Plane-strain twin
   * - ``GradientEnhancedFiniteStrainSQCNIxNSNI/Axisymmetric/Quad``
     - ``GradientEnhancedFiniteStrainSQCNIxNSNI/PlaneStrain/Quad`` (:doc:`gradientenhancedfinitestrainparticle`)
   * - ``DisplacementSQCNI/Axisymmetric/Quad``
     - ``DisplacementSQCNI/PlaneStrain/Quad`` (:doc:`displacementparticle`)

and they are the same classes: the section of the particle, ``Marmot::Discretizations::Meshfree::Section``
(``Cartesian`` or ``Axisymmetric``), is a template parameter, and what it changes is collected in
``ParticleSectionIntegrals`` (the integrals over the smoothing domain) and ``Hoop`` (the hoop terms of the weak form),
both in ``MarmotParticleSectionIntegrals.h``. A host learns the section of a particle from ``Particle::getSection()``, e.g. to
integrate the contact points on the faces of the particles over the same rings. This page derives what the
axisymmetric section changes; everything else -- the updated-Lagrangian kinematics, the material point, the nonlocal
balance, the explicit life cycle -- is that of the plane-strain particles. The conventions follow the axisymmetric finite elements ``CX8UL`` / ``CX8RUL``
(:doc:`displacementfinitestrainelement`): the material sees the 3D state with the hoop stretch as its
:math:`F_{33}`, and every integral carries the ring weight :math:`2\pi r`.

A particle may touch the axis with a vertex or a face, as the particles of a solid cylinder do: the ring weight
:math:`2\pi r` vanishes there, so that a face on the axis contributes nothing (the axis is no boundary of the body of
revolution), and the center of the particle has :math:`\bar r > 0`, so that its hoop operator
:math:`N_A/\bar r` stays finite. A particle with a vertex at :math:`r < 0` is refused when it is created, and an
increment that pulls a particle through the axis (a hoop stretch :math:`\le 0`) throws, which a solver may answer with
a smaller increment.

The ring measure
----------------

With the volume element :math:`dV = 2\pi r\,dA` of the ring, the volume of a particle with the domain
:math:`\Omega` in the meridian half plane is, by Pappus' theorem,

.. math::

   V = 2\pi \int_\Omega r\,dA = 2\pi\,A\,\bar r ,

with its area :math:`A` and the radius :math:`\bar r` of its area centroid, where its material point sits. :math:`V_0`,
the volume of the undeformed particle, is the integration weight of every term of the particle (the internal force,
the nonlocal balance, the lumped inertia, damping and momentum, the body force), as the undeformed area is in plane
strain. The lumped mass of node :math:`A` is :math:`\rho_0\,T_A\,V_0 = \rho_0\,T_A\,2\pi\int_{\Omega_0} R\,dA`.

Kinematics: the hoop stretch
----------------------------

The increment of the deformation gradient is block diagonal: the in-plane block of the plane particle and the
increment of the hoop stretch,

.. math::

   \Delta F_{ij} = \delta_{ij} + \Delta q_{Bi}\,\frac{\partial N_B}{\partial Y_j}\quad (i, j \in \{r, z\}), \qquad
   \Delta F_{\theta\theta} = \frac{\bar x_r}{\bar r} = 1 + \frac{\Delta u_r}{\bar r}
   = 1 + \Delta q_{Br}\,\frac{N_B}{\bar r},

with the radial displacement increment :math:`\Delta u_r = N_B\,\Delta q_{Br}` at the particle center, its radius
:math:`\bar r` in the intermediate configuration :math:`\boldsymbol Y` and :math:`\bar x_r` after the increment. The
material point composes :math:`\boldsymbol F = \Delta\boldsymbol F\,\boldsymbol F_n` in 3D, so that its
:math:`F_{\theta\theta}` is the product of the hoop stretch increments, the current over the undeformed radius
:math:`x_r/R`, and the material is evaluated at this 3D deformation gradient through
``computePlaneStrain`` (as ``CX8UL``). Why the hoop operator is the point value :math:`N_B/\bar r` is shown below with
the smoothed gradient.

The weak form
-------------

The virtual work of the Kirchhoff stress over the undeformed ring, :math:`\int \boldsymbol\tau : \delta\boldsymbol d\,
dV_0`, gains the hoop term :math:`\tau_{\theta\theta}\,\delta d_{\theta\theta}` with :math:`\delta d_{\theta\theta} =
\delta x_r / x_r = \delta u_r / (\bar r\,\Delta F_{\theta\theta})`. The displacement residual of node :math:`A` is

.. math::

   r^U_{Aj} = \frac{\partial T_A}{\partial x_i}\,\tau_{ij}\,V_0
   + \delta_{jr}\,\frac{T_A}{\bar r\,\Delta F_{\theta\theta}}\,\tau_{\theta\theta}\,V_0
   + \rho_0\,a_j\,T_A\,V_0 ,

the hoop term of ``CX8UL``, :math:`N_A\,\tau_{33} / (F_{33}\,r)\,J_0 w`. The Helmholtz balance of a nonlocal field is a
balance of a scalar of revolution: its gradient has no hoop component, and only its weight changes to :math:`V_0`.
The tangent gains, through :math:`\partial\Delta F_{\theta\theta}/\partial\Delta u_r = 1/\bar r`, the derivatives of
the hoop row with respect to the in-plane gradient and to the radial displacement increment, and those of all rows
with respect to :math:`\Delta F_{\theta\theta}` through the material
(:math:`\partial\tau_{ij}/\partial\Delta F_{\theta\theta}`, :math:`\partial L_K/\partial\Delta F_{\theta\theta}`):

.. math::

   \frac{\partial r^U_{Ar}}{\partial(\partial\Delta u_k/\partial Y_l)} \mathrel{+}=
     \frac{T_A V_0}{\bar r\,\Delta F_{\theta\theta}}\,\frac{\partial\tau_{\theta\theta}}{\partial\Delta F_{kl}}, \qquad
   \frac{\partial r^U_{Ar}}{\partial\Delta u_r} \mathrel{+}= \frac{T_A V_0}{\bar r^2\,\Delta F_{\theta\theta}}
     \Bigl(\frac{\partial\tau_{\theta\theta}}{\partial\Delta F_{\theta\theta}}
     - \frac{\tau_{\theta\theta}}{\Delta F_{\theta\theta}}\Bigr),

.. math::

   \frac{\partial r^U_{Aj}}{\partial\Delta u_r} \mathrel{+}= \frac{V_0}{\bar r}\,\frac{\partial T_A}{\partial x_i}\,
     \frac{\partial\tau_{ij}}{\partial\Delta F_{\theta\theta}}, \qquad
   \frac{\partial r^U_{Ar}}{\partial\Delta\bar N_K} \mathrel{+}= \frac{T_A V_0}{\bar r\,\Delta F_{\theta\theta}}\,
     \frac{\partial\tau_{\theta\theta}}{\partial\bar N_K}, \qquad
   \frac{\partial r^N_{AK}}{\partial\Delta u_r} \mathrel{+}= -\frac{T_A V_0}{\bar r}\,
     \frac{\partial L_K}{\partial\Delta F_{\theta\theta}} .

All of them depend on the nodes :math:`B` only through :math:`N_B` and :math:`\partial N_B/\partial\boldsymbol Y`: the
particle describes its stiffness node by node (``Particle::describeStiffness``) exactly as in plane strain, with the
hoop terms in the derivatives with respect to the interpolated increment. The bound of the absolute stiffness for the
critical time step (``AbsoluteStiffnessBound``, ``addBoundOfAbsoluteStiffnessTimes``) therefore holds unchanged, with
the ring weight in the stiffness and in the lumped mass alike.

The smoothed gradient of the ring
---------------------------------

SQCNI replaces the gradient at the particle center by its average over the smoothing domain :math:`\Omega_s`, written
as a boundary integral by the divergence theorem (Chen, Wu, Yoon and You, 2001), see :doc:`meshfreeparticles`; the
notation is that page's, :math:`N_{A,i} = \partial N_A/\partial Y_i`. Over the ring, the average is weighted with
:math:`2\pi r`, and the divergence theorem of the weighted integrand gains an area term,
:math:`r\,N_{A,r} = (r N_A)_{,r} - N_A`:

.. math::

   \tilde N_{A,i} = \frac{1}{|\Omega_s|}\int_{\Omega_s} N_{A,i}\,2\pi r\,dA
   = \frac{2\pi}{|\Omega_s|}\Bigl(\oint_{\partial\Omega_s} r\,N_A\,n_i\,ds - \delta_{ir}\int_{\Omega_s} N_A\,dA\Bigr),
   \qquad |\Omega_s| = 2\pi\int_{\Omega_s} r\,dA ,

with the area element :math:`dA` and the line element :math:`ds` of the meridian half plane, and the ring volume
:math:`|\Omega_s|` of the smoothing domain. On a face :math:`f`, :math:`ds` integrates to :math:`dA_f`, the length of
the face (the area of unit thickness in plane strain), and :math:`n_{f,i}\,dA_f` is its area vector
(``LagrangeCell::boundarySurfaceVector``). The hoop strain is averaged alike, and its area integral is the same one:

.. math::

   \tilde H_A = \frac{1}{|\Omega_s|}\int_{\Omega_s}\frac{N_A}{r}\,2\pi r\,dA
   = \frac{2\pi}{|\Omega_s|}\int_{\Omega_s} N_A\,dA .

**What is adopted.** The face integrals are evaluated with two Gauss points :math:`\boldsymbol Y_g` per face, at the
radius :math:`r_g`, which share the area vector :math:`n_{f,i}\,dA_f` equally; the area integral with one point, the
area centroid of :math:`\Omega_s`, which is the particle center: :math:`\int_{\Omega_s} N_A\,dA \approx
A_s\,N_A(\bar{\boldsymbol Y})`. With :math:`|\Omega_s| = 2\pi A_s\bar r`, the hoop operator becomes the point value,
:math:`\tilde H_A = N_A(\bar{\boldsymbol Y})/\bar r`, and

.. math::

   \tilde N_{A,i} = \frac{1}{|\Omega_s|}\sum_f\sum_{g=1}^{2} 2\pi r_g\,N_A(\boldsymbol Y_g)\,
   \frac{n_{f,i}\,dA_f}{2} - \delta_{ir}\,\frac{N_A(\bar{\boldsymbol Y})}{\bar r} .

For :math:`2\pi r_g = 1`, one point per face (the face center) and without the area term, this is the smoothed
gradient of the plane particles.

Both properties of SQCNI carry over exactly:

- **The integration constraint** (the linear patch test). It is the discrete divergence theorem of the 3D field
  :math:`N_A\boldsymbol e_r`, whose divergence is :math:`N_{A,r} + N_A/r`, and of :math:`N_A\boldsymbol e_z`:

  .. math::

     \sum_L V_L\,\bigl(\tilde N_{A,r} + \tilde H_A\bigr)
     = \sum_{f\in\Gamma}\sum_{g=1}^{2} 2\pi r_g\,N_A(\boldsymbol Y_g)\,\frac{n_{f,r}\,dA_f}{2} , \qquad
     \sum_L V_L\,\tilde N_{A,z}
     = \sum_{f\in\Gamma}\sum_{g=1}^{2} 2\pi r_g\,N_A(\boldsymbol Y_g)\,\frac{n_{f,z}\,dA_f}{2} ,

  summed over the particles :math:`L` with :math:`V_L = |\Omega_s|`, and over the faces :math:`f` on the boundary
  :math:`\Gamma` of the body. It holds because, for each particle, :math:`V_L\,\tilde N_{A,r}` is its face sum minus
  :math:`2\pi A_L N_A(\bar{\boldsymbol Y}_L)`, and :math:`V_L\,\tilde H_A` is exactly
  :math:`+2\pi A_L N_A(\bar{\boldsymbol Y}_L)`: the area terms cancel within each particle, whatever quadrature
  evaluates :math:`\int N_A\,dA`, as long as both operators use the same value. The face sums of the faces two
  particles share cancel between them (the same Gauss points, opposite normals), and what remains is the face sum
  over :math:`\Gamma`. A face on the axis has :math:`r_g = 0` and is no part of :math:`\Gamma`. The cancellation
  between neighbours needs smoothing domains that tile the body, i.e. the update ``DeformationGradient`` (SQCNI); the
  axisymmetric particles are registered with it only.

  In the patch test, the stress of a linear displacement field that is in equilibrium without a body force is
  uniform with :math:`\sigma_{rr} = \sigma_{\theta\theta}` and :math:`\sigma_{rz} = 0` (e.g. :math:`u_r = a\,r`,
  :math:`u_z = b\,z`). The radial internal force of node :math:`A`,
  :math:`\sum_L V_L(\tilde N_{A,r}\,\sigma_{rr} + \tilde H_A\,\sigma_{\theta\theta})`, is then :math:`\sigma_{rr}`
  times the first sum, and the axial one :math:`\sigma_{zz}` times the second: zero for an interior node, and for a
  node on :math:`\Gamma` the force of the traction :math:`\boldsymbol\sigma\boldsymbol n` integrated with the same
  face rule. A load on :math:`\Gamma` is therefore integrated with this rule (``ParticleSectionIntegrals::addPressure``),
  so that the patch test with a traction on the boundary passes exactly as well.
- **Exactness for linear fields.** For a displacement the approximation reproduces, :math:`u = a + b\,r + c\,z`, the
  smoothed gradient and hoop strain are the exact ring averages of the exact ones: the face integrands
  :math:`r\,u` are quadratic along a straight face (two Gauss points integrate them exactly; the one-point rule of the
  plane particles would not), and the area integrand :math:`u` is linear (the centroid rule integrates it exactly).
  E.g., for :math:`u_r = a\,r` both :math:`\tilde N_{A,r}\,q_{Ar}` and :math:`\tilde H_A\,q_{Ar}` give :math:`a`, so
  that the particle sees the uniform state the patch test prescribes. The hoop stretch :math:`1 + \Delta u_r(\bar{\boldsymbol Y})/\bar r`
  of a linear increment is the ratio of the moved and the old radius of the centroid, so that the ring volume of the
  moved smoothing domain is :math:`J\,|\Omega_s|` with the 3D Jacobian.

The gradient of a nonlocal field is smoothed with the same :math:`\tilde N_{A,i}`; it enters the Helmholtz balance in
:math:`(r, z)` only. These are the "axisymmetric patch tests" of the module tests: a uniform radial expansion
:math:`u_r = a\,r`, an axial strain :math:`u_z = b\,z` and their combination leave interior particles without force
to round-off, and with the traction :math:`\boldsymbol\sigma\boldsymbol n` as a pressure on the boundary every
particle; a shear :math:`u_z = c\,r` is reproduced exactly (its stress needs the body force
:math:`b_z = -\sigma_{rz}/r`, so its internal force is checked against the smoothed operators, not against zero).

References. The divergence-theorem smoothing and its integration constraint are those of SCNI (Chen, Wu, Yoon and
You, *A stabilized conforming nodal integration for Galerkin mesh-free methods*, IJNME 50, 2001); the integration
constraint as the condition of linear exactness is that of VCI (Chen, Hillman and Rüter, *An arbitrary order
variationally consistent integration for Galerkin meshfree methods*, IJNME 96, 2013). Their ring-weighted form, with
the area term and its cancellation against the smoothed hoop strain, is derived here, for the quadrature written
above.

The stabilization (NSNI)
------------------------

NSNI (Hillman and Chen, *An accelerated, convergent, and stable nodal integration in Galerkin meshfree methods for
linear and nonlinear mechanics*, IJNME 107, 2016) adds the second-order term of the Taylor expansion of the integrand
about the particle center: for two factors :math:`a`, :math:`b` of the integrand,

.. math::

   \int_\Omega (a + \nabla a\cdot\boldsymbol\xi)(b + \nabla b\cdot\boldsymbol\xi)\,2\pi r\,dA
   = a\,b\,V + (a\,\nabla b + b\,\nabla a)\cdot\boldsymbol m + \nabla a\cdot\boldsymbol M\,\nabla b ,
   \qquad \boldsymbol M = 2\pi\int_\Omega \boldsymbol\xi\otimes\boldsymbol\xi\,r\,dA ,

with :math:`\boldsymbol\xi = \boldsymbol Y - \bar{\boldsymbol Y}` and the ring's first moment
:math:`\boldsymbol m = 2\pi\int\boldsymbol\xi\,r\,dA`, which does not vanish about the area centroid. The first term is
what the smoothed gradients integrate; the second is a first-order (consistency) term, which the smoothed gradients
already account for -- adding it would break the integration constraint, since it does not vanish for the uniform
stress of a linear field. The stabilization keeps the third term, the one of the plane particles, with the **ring
second moments** :math:`\boldsymbol M_0` of the undeformed particle about its area centroid (integrated with a
:math:`3\times 3` Gauss rule, exact for a quadrilateral), pushed forward with :math:`\boldsymbol F_n` as in plane
strain. Its second gradients of the shape functions are smoothed over the ring like the first ones,

.. math::

   \hat N_{A,IJ} = \frac{2\pi}{|\Omega_s|}\Bigl(\oint r\,N_{A,I}\,n_J\,ds - \delta_{Jr}\oint N_A\,n_I\,ds\Bigr)

(symmetrized; :math:`\int N_{A,I}\,dA = \oint N_A\,n_I\,ds`), with the Gauss points of the faces. For
a linear field both terms cancel (:math:`\oint r\,n_J\,ds = \delta_{Jr}A`), so that the stabilization vanishes in the
patch test, as in plane strain. As there, the stress gradient of the stabilization is linearized in the in-plane
deformation gradient only; its parts through the nonlocal field and through the hoop stretch are left out (the
stabilization acts on the in-plane gradient, which carries every oscillating mode of the nodal integration).

Loads and corrections
---------------------

- ``BODYFORCE``: :math:`\boldsymbol f_A \mathrel{-}= T_A\,\boldsymbol b\,V_0` with the ring volume.
- ``PRESSURE`` on a face: integrated with the face rule of the smoothing, the two Gauss points :math:`\boldsymbol Y_g`,
  as the integration constraint requires (see above). Each point carries its share of the area vector of the face,
  pushed to the current configuration by Nanson's formula in the half plane, times its ring after the increment:

  .. math::

     \boldsymbol f_g = -\Delta J\,\Delta\boldsymbol F^{-T} p\,\frac{\boldsymbol n_f\,dA_f}{2}\;2\pi x_r(\boldsymbol Y_g),
     \qquad \boldsymbol x(\boldsymbol Y) = \bar{\boldsymbol x} + \Delta\boldsymbol F\,(\boldsymbol Y - \bar{\boldsymbol Y}),
     \qquad f_{Ai} \mathrel{-}= \sum_g T_A(\boldsymbol Y_g)\,f_{g,i} .

  The face moves with the motion of the particle (the motion Nanson's formula assumes), with the center at
  :math:`\bar x_r = \bar r\,\Delta F_{\theta\theta}`, so that :math:`2\pi x_r(\boldsymbol Y_g) = 2\pi r_g\,
  \Delta F_{\theta\theta}(\boldsymbol Y_g)` is the ring of the point with its own hoop stretch increment. It is exact
  for a linear increment; its tangent differentiates :math:`\Delta\boldsymbol F` and :math:`x_r(\boldsymbol Y_g)`,
  :math:`\partial x_r(\boldsymbol Y_g)/\partial\Delta q_{Br} = N_B(\bar{\boldsymbol Y}) + \tilde N_{B,j}\,(Y_{g,j} -
  \bar Y_j)`.
- ``CWFCORRECTION`` and VCI are not available for the axisymmetric particles: a request throws. VCI (Chen, Hillman
  and Rüter, 2013) would have to be rederived for the ring measure and the hoop term; the ring-weighted SQCNI already
  satisfies the linear integration constraint by itself.

Tests
-----

``TestAxisymmetricParticles`` (module ``GradientEnhancedFiniteStrainParticle``) checks both particles, on a patch away
from the axis and on one at the axis, whose particles touch it with a face: the axisymmetric patch tests (interior
forces vanish to round-off for :math:`u_r = a r`, :math:`u_z = b z` and their combination; the smoothed operators
reproduce :math:`u_z = c\,r`, :math:`u_r = c\,z` exactly), the patch test with the traction
:math:`\boldsymbol\sigma\boldsymbol n` as a ``PRESSURE`` on the boundary (every node is balanced to round-off) and the
tangent of the pressure; the tangent against finite differences, the lumped mass against
:math:`\rho_0\,2\pi\int R\,dA`, a thick-walled cylinder under internal pressure against Lamé's solution and against a
mesh of ``CX8RUL`` elements, and the bound of the absolute stiffness against the largest eigenvalue of an assembled
patch.

Implementation
--------------

The particle classes are documented on the pages of their formulations, :doc:`displacementparticle` and
:doc:`gradientenhancedfinitestrainparticle`, and ``ParticleSectionIntegrals`` on :doc:`meshfreeparticles`. The section,
the hoop terms and the out-of-plane parts of the 2D material points:

.. doxygenenum:: Marmot::Discretizations::Meshfree::Section

.. doxygennamespace:: Marmot::Discretizations::Meshfree::Hoop

.. doxygenstruct:: Marmot::Discretizations::MaterialPoints::OutOfPlaneTangents
   :allow-dot-graphs:

.. doxygenfunction:: Marmot::Discretizations::MaterialPoints::computeOutOfPlaneTangents
