Displacement Finite Element
===========================

Theory
------

This element implements a small-strain displacement-based Finite Element formulation.
The incremental strain is computed from the displacement increment via the kinematic
matrix B, and internal forces and consistent tangents are evaluated by Gauss quadrature.

.. math::

   \Delta \eps = \mathbf{B}\, \Delta \mathbf{u}, \qquad
   \mathbf{K}_e = \sum_{qp} \mathbf{B}^\mathsf{T}\, \mathbf{C}\, \mathbf{B}\, J_0 w, \qquad
   \mathbf{P}_e = \sum_{qp} \mathbf{B}^\mathsf{T}\, \sig\, J_0 w.

where :math:`\Delta \eps` is the linearized strain tensor in Voigt notation,
:math:`\mathbf{B}` is the strain–displacement matrix, :math:`\Delta \mathbf{u}` the nodal displacement increment;
:math:`\mathbf{K}_e` the element tangent stiffness, :math:`\mathbf{C}` the consistent material tangent;
:math:`J_0 = \det \mathbf{J}` the Jacobian determinant and :math:`w` the quadrature weight;
:math:`\mathbf{P}_e` the element internal force vector, :math:`\sig` the Cauchy stress tensor in Voigt notation,
and :math:`\sum_{qp}` denotes summation over quadrature points.

The mass matrix and body-force vector are obtained from the expressions:

.. math::

   \mathbf{M}_e = \sum_{qp} \rho\, \mathbf{N}^\mathsf{T} \mathbf{N}\, J_0 w, \qquad
   \mathbf{P}_e^{(b)} = \sum_{qp} \mathbf{N}^\mathsf{T} \mathbf{f}\, J_0 w.

where :math:`\mathbf{M}_e` is the consistent mass matrix, :math:`\rho` the mass density,
:math:`\mathbf{N}` the shape-function matrix of the displacement field, and
:math:`\mathbf{f}` the body-force vector per unit volume.

The consistent mass is integrated with its own Gauss rule, exact for degree 2p with p the order of the
shape functions (``gauss( 2p )``: 3x3x3 for a 20-node hexahedron, 3x3 for an 8-node quadrilateral, 4 points
for a linear tetrahedron), also for a reduced-integration element: its
own rule cannot produce a mass matrix of full rank -- 8 points give a 20-node hexahedron a mass of
rank 24 out of 60 -- and a singular mass leaves the modes of its null space without inertia and
makes an equilibrium solve for the initial acceleration, :math:`\mathbf{M}\,\mathbf{a}_0 = \mathbf{R}`,
ill-posed. The density at a point of that rule is taken from the nearest quadrature point of the
element, which is exact for a density that is constant over the element; the section factor (thickness,
cross section or :math:`2\pi r`) is evaluated at the point itself. The same applies to the
finite-strain (UL) and the gradient-enhanced displacement elements, for the latter to the
micro-inertia of the non-local block as well.

Surface tractions and pressures on a boundary face :math:`\Gamma_e` are integrated as

.. math::

   \mathbf{P}_e^{(t)} = \int_{\Gamma_e} \mathbf{N}^\mathsf{T} \, \mathbf{t} \, \mathrm{d}\Gamma, \qquad
   \mathbf{P}_e^{(p)} = - \int_{\Gamma_e} p \, \mathbf{N}^\mathsf{T} \, \mathbf{n} \, \mathrm{d}\Gamma.

where :math:`\Gamma_e` is the element boundary, :math:`\mathbf{n}` its outward unit normal,
:math:`\mathbf{t}` the surface traction vector, :math:`p` the pressure magnitude, and
:math:`\mathbf{N}` the shape-function matrix. The minus sign for pressure follows
an outward-normal convention.

For 2D and 1D problems, the integration measure :math:`J_0 w` is scaled by thickness
:math:`t` or cross-sectional area :math:`A`, respectively: :math:`J_0 w \leftarrow J_0 w\, t`
(2D), :math:`J_0 w \leftarrow J_0 w\, A` (1D).
Here, :math:`t` denotes thickness and :math:`A` denotes cross-sectional area.

Constitutive updates are carried out per quadrature point using a Marmot material model:

- 3D solid

- Plane stress

- Plane strain

- Axisymmetric

- 1D uniaxial stress

Each quadrature point stores stress, strain and a material state vector; the element
accumulates the history incrementally and may suggest a reduced time step if required
by the material.


Axisymmetric elements
---------------------

``CAX4`` (4 nodes), ``CAX8`` and ``CAX8R`` (8 nodes) are ``DisplacementFiniteElement<2, nNodes>`` with the section
``Axisymmetric``: the body of revolution is meshed in its meridional plane, with the radius :math:`r = x_1` and the
axial coordinate :math:`z = x_2`. The displacement :math:`(u_r, u_z)` does not depend on the angle :math:`\theta`, and
the strain has the hoop component :math:`\varepsilon_{\theta\theta} = u_r / r` in addition to the in-plane ones:

.. math::

   \Delta\eps = \begin{bmatrix} \Delta\varepsilon_{rr} \\ \Delta\varepsilon_{zz} \\
   \Delta\varepsilon_{\theta\theta} \\ \Delta\gamma_{rz} \end{bmatrix}
   = \sum_i \begin{bmatrix} \partial N_i / \partial r & 0 \\ 0 & \partial N_i / \partial z \\
   N_i / r & 0 \\ \partial N_i / \partial z & \partial N_i / \partial r \end{bmatrix}
   \begin{bmatrix} \Delta u_{r,i} \\ \Delta u_{z,i} \end{bmatrix}.

Volume and boundary are integrated over the full ring,

.. math::

   \mathrm{d}V = 2\pi r\, \mathrm{d}A, \qquad \mathrm{d}\Gamma = 2\pi r\, \mathrm{d}s,

i.e., :math:`J_0 w = w \det\mathbf{J}\, 2\pi r` at a quadrature point, with :math:`r` interpolated from the nodes,
and :math:`2\pi r` at each point of the boundary and the mass rules.

The material is called through the 3D interface with the strain increment
:math:`(\Delta\varepsilon_{rr}, \Delta\varepsilon_{zz}, \Delta\varepsilon_{\theta\theta}, \Delta\gamma_{rz}, 0, 0)`
(``axisymmetricVoigtToVoigt``). The tangent of the element is the upper left :math:`4 \times 4` block of the 3D
tangent, which is exact for any material: :math:`\gamma_{r\theta} = \gamma_{z\theta} = 0` is kinematic, as
:math:`\varepsilon_{33} = 0` in plane strain.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Element**
     - **Nodes**
     - **Quadrature**
     - **Zero energy modes**
   * - ``CAX4``
     - 4
     - Gauss 2x2 (``gaussExactToDegree( 2 )``)
     - 1 (axial translation)
   * - ``CAX8``
     - 8
     - Gauss 3x3 (``gaussExactToDegree( 4 )``)
     - 1 (axial translation)
   * - ``CAX8R``
     - 8
     - Gauss 2x2 (``gaussExactToDegree( 3 )``)
     - 2 (axial translation, one hourglass mode)

Notes:

- Internal and external forces, reactions, energies and masses are those of the full :math:`360^\circ` ring.
- The radial displacement on the axis, :math:`u_r = 0` at :math:`r = 0`, is not imposed by the element; prescribe it.
- Mesh in :math:`r \geq 0` only. ``initializeYourself()`` rejects a node at :math:`r < 0` and a quadrature point at
  :math:`r \leq 0`.
- A thickness assigned as element property is ignored.
- ``CAX8R`` has one spurious zero energy (hourglass) mode per element, like ``CPE8R``, in addition to the axial
  translation.
- The artificial bulk viscosity acts on :math:`\Delta\varepsilon_{rr} + \Delta\varepsilon_{zz} +
  \Delta\varepsilon_{\theta\theta}`, including the hoop strain.
- The energy gradients (:doc:`energygradients`) are available with a material that provides energy tangents.


Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::DisplacementFiniteElement
   :allow-dot-graphs:
