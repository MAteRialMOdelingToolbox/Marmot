General Gradient-Enhanced Displacement Finite Element
=====================================================

Theory
------
This element implements a small strain gradient enhanced displacement based Finite Element formulation.
The formulation couples the standard displacement field with one or multiple nonlocal scalar fields.

In addition to the standard static equilibrium equation, one or multiple nonlocal scalar fields are solved for.

.. math::

   \nabla \sig + \boldsymbol f = \boldsymbol 0

Here :math:`\sig` is the Cauchy stress tensor and :math:`\boldsymbol f` the body force vector per unit volume.
Each nonlocal variable introduces an additional balance equation, which is solved simultaneously with the static equilibrium equation.
The general balance equation for the i-th nonlocal variable :math:`\knl_i` is defined as

.. math::

   \knl_i - \nabla ( c( \knl_i )\, \nabla \knl_i ) = \kl_i (\eps,\, \knl_i )\, ,

where :math:`\knl` is the nonlocal variable, :math:`\kl` the local equivalent, :math:`c( \knl_i )` the gradient interaction parameter, and :math:`\kl_i` is the local
driving variable for the nonlocal field.
The element formulation supports a variable length scale parameter which can be used to model damage dependent interactions and may depend on the state of the
nonlocal field: :math:`c = c(\knl)`.
This implementation is versatile and is suited for both implicit gradient enhanced models and phase field formulations to simulate failure and strain localization.
Depending on the specific material model, this field is governed either by the screened Poisson equation for implicit gradient regularisation,

.. math::

   \knl - c \, \nabla^2 \knl = \kl(\eps)\, ,

or by the phase field crack evolution equation:

.. math::

   d - l^2 \, \nabla^2 d = s(d, \mathcal{Y})\, .

For the screened Poisson equation :math:`\knl` is the nonlocal equivalent strain-like variable, :math:`\kl` the local equivalent, and :math:`c=l^2` the
gradient length scale parameter.
In the context of phase field models, :math:`d` represents the phase field variable, :math:`l` is the internal length scale, and :math:`s(d, \mathcal{Y})`
is a local source function driven by the thermodynamic driving force :math:`\mathcal{Y}`.
Homogeneous Neumann boundary conditions are assumed for the nonlocal scalar fields on the element boundary, i.e.,
:math:`(\nabla \knl \boldsymbol)^T\, n = 0` or :math:`(\nabla d)^T \, \boldsymbol n = 0`, where :math:`\boldsymbol n` is the outward unit normal.

The nodal unknown vector uses a field wise structure:

.. math::
   \mathbf{q} = [\mathbf{\qu}, \mathbf{\qk}]^\mathsf{T}

The incremental strain is computed from the displacement increment via the kinematic matrix :math:`\mathbf{B}`.
The nonlocal field is interpolated using the shape function for the nonlocal field :math:`\Nk` which may be configured to be an order lower than the displacement
shape functions.

.. math::

   \Delta \eps = \mathbf{B}\, \Delta \mathbf{\qu}, \qquad \knl = \Nk\, \mathbf{\qk}
   %  \knl = \Nk_A\, \qk_A,\, ,

Here :math:`\mathbf{B}` is the strain-displacement matrix, :math:`\Delta \mathbf{\qu}` the nodal displacement increment,
and :math:`\Delta \eps` is the linearized strain tensor in Voigt notation.

The discretized weak form of the governing equations is given by

.. math::
   \delta \qu_{Aj} \left( \int_{V} \mathbf{N}_{A,i}\, \sig_{ij} \, dV -\int_{V} \mathbf{N}_{A}\, f_j \, dV -\int_{S_\mathrm{t}}
   \mathbf{N}_A\, \bar{t}_j \, d S_\mathrm{t} \right) = \delta \qu_{Aj}\, \ru_{Aj} = 0\, ,

   \delta \qk_A \int_V \mathbf{N}_A \, \knl + c \, \mathbf{N}_{A,i}\, \knl_{,i} -
   \mathbf{N}_A\, \kl \, dV = \delta \qk_A \, \rk_A = 0\, ,

where :math:`\delta \qu_{Aj}` and :math:`\delta \qk_A` are the virtual nodal variations of the displacement and nonlocal field, respectively.
:math:`\ru_{Aj}` and :math:`\rk_A` denote the nodal residuals associated with the respective degrees of freedom.
Uppercase subscripts, for example :math:`(\bullet)_A`, represent nodal indices that denote specific points in the discretised domain, while lowercase subscripts,
for example :math:`(\bullet)_i`, denote the spatial dimensions.
Internal forces and consistent tangents are evaluated by Gauss quadrature:

.. math::

   \mathbf{K}_e = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{K}_{uk} \\ \mathbf{K}_{ku} & \mathbf{K}_{kk} \end{bmatrix},\qquad
   \mathbf{\fuint} = \sum_{qp} \mathbf{B}^\mathsf{T}\, \sig\, J_0\, w_{qp}\, ,

   \mathbf{\fk} = \sum_{qp} \left (\mathbf{\Nk}^\mathsf{T}\, \knl\, + c\, \partial_\mathbf{x} \mathbf{\Nk}^\mathsf{T}\,
   \partial_\mathbf{x} \mathbf{\Nk}\, \mathbf{\qk} - \mathbf{\Nk}^\mathsf{T}\, \kl \right ) J_0\, w_{qp} \, .

Here :math:`\mathbf{K}_e` the element tangent stiffness in block structure, :math:`J_0 = \det \mathbf{J}` the Jacobian determinant and :math:`w_{qp}`
the quadrature weight, :math:`\mathbf{\fuint}` the element internal force vector for the displacement degrees of freedom,
:math:`\mathbf{\fk}` the element force vector for the nonlocal degrees of freedom, :math:`\sig` the Cauchy stress tensor in Voigt notation,
and :math:`\sum_{qp}` denotes summation over quadrature points.

The stiffness submatrices :math:`\mathbf{K}_{uu}`, :math:`\mathbf{K}_{u\knl}`, :math:`\mathbf{K}_{\knl u}` and :math:`\mathbf{K}_{\knl\knl}` are
evaluated using the following expressions:

.. math::

   \begin{aligned}
   & \mathbf{K}_{uu} = \sum_{qp} \mathbf{B}^\mathsf{T} \frac{\partial \mathbf{\sig}}{\partial \mathbf{\eps}} \mathbf{B}\, J_0\, w_{qp}, \\
   & \mathbf{K}_{u\knl} = \sum_{qp} \mathbf{B}^\mathsf{T} \frac{\partial \mathbf{\sig}}{\partial \knl} \mathbf{\Nk}\, J_0\, w_{qp}, \\
   & \mathbf{K}_{\knl u} = -\sum_{qp} \mathbf{\Nk}^\mathsf{T} \frac{\partial \kl}{\partial \mathbf{\eps}} \mathbf{B}\, J_0\, w_{qp}, \\
   & \mathbf{K}_{\knl\knl} = \sum_{qp} \left( \mathbf{\Nk}^\mathsf{T}\, \mathbf{\Nk} + c \, \partial_\mathbf{x} \mathbf{\Nk}^\mathsf{T} \,
   \partial_\mathbf{x} \mathbf{\Nk} + \frac{\partial c}{\partial \knl} \, \partial_\mathbf{x} \mathbf{\Nk}^\mathsf{T} \, \partial_\mathbf{x}
   \mathbf{\Nk}\, \qk \, \mathbf{\Nk} - \mathbf{\Nk}^\mathsf{T}\, \frac{\partial \kl}{\partial \knl} \right) J_0\, w_{qp}\, .
   \end{aligned}

The derivative :math:`\partial c / \partial \knl` is only non-zero if a variable length scale is used, and the derivative
:math:`\partial \kl / \partial \knl` is only required in phase field formulations where it represents the derivative of the
source function with respect to the crack phase field :math:`\frac{\partial s(d, \mathcal{Y})}{\partial d}`.

The body-force vector is obtained using

.. math::

   \mathbf{\fextb} = \sum_{qp} \mathbf{N}^\mathsf{T} \mathbf{f}\, J_0 w\, ,

and surface pressures on a boundary face :math:`\Gamma_e` are integrated as

.. math::

   \mathbf{\fextp} = - \int_{\Gamma_e} p \, \mathbf{N}^\mathsf{T} \, \mathbf{n} \, \mathrm{d}\Gamma.


where :math:`\Gamma_e` is the element boundary, :math:`\mathbf{n}` its outward unit normal, :math:`\mathbf{N}` the shape-function matrix
of the displacement field, :math:`p` the pressure magnitude, and :math:`\mathbf{f}` the body-force vector per unit volume.
The minus sign for pressure follows an outward-normal convention.

For plane stress and plane strain, the integration measure :math:`J_0 w` is scaled by the thickness :math:`t`,
:math:`J_0 w \leftarrow J_0 w\, t`; axisymmetric elements integrate over the ring, see below.

Constitutive updates are carried out per quadrature point using a Marmot material model, through its 3D interface
(solid, plane strain and axisymmetric) or its plane stress iteration, for the sections

- 3D solid

- Plane stress

- Plane strain

- Axisymmetric

Each quadrature point stores stress, strain and a material state vector; the element
accumulates the history incrementally and may suggest a reduced time step if required
by the material.

Axisymmetric elements
---------------------

``GCAX4`` (4 nodes), ``GCAX8`` and ``GCAX8R`` (8 nodes) are ``GeneralGradientEnhancedDisplacementFiniteElement<2,
nNodes>`` with the section ``Axisymmetric``, the gradient-enhanced counterparts of ``CAX4``, ``CAX8`` and ``CAX8R``
(:doc:`displacementfiniteelement`). The body of revolution is meshed in its meridional plane, with the radius
:math:`r = x_1` and the axial coordinate :math:`z = x_2`. The displacement :math:`(u_r, u_z)` and the nonlocal fields do
not depend on the angle :math:`\theta`, and the strain has the hoop component
:math:`\varepsilon_{\theta\theta} = u_r / r` in addition to the in-plane ones:

.. math::

   \Delta\eps = \begin{bmatrix} \Delta\varepsilon_{rr} \\ \Delta\varepsilon_{zz} \\
   \Delta\varepsilon_{\theta\theta} \\ \Delta\gamma_{rz} \end{bmatrix}
   = \sum_i \begin{bmatrix} \partial N_i / \partial r & 0 \\ 0 & \partial N_i / \partial z \\
   N_i / r & 0 \\ \partial N_i / \partial z & \partial N_i / \partial r \end{bmatrix}
   \begin{bmatrix} \Delta u_{r,i} \\ \Delta u_{z,i} \end{bmatrix}.

Volume and boundary are integrated over the full ring,

.. math::

   \mathrm{d}V = 2\pi r\, \mathrm{d}A, \qquad \mathrm{d}\Gamma = 2\pi r\, \mathrm{d}s,

i.e., :math:`J_0 w = w \det\mathbf{J}\, 2\pi r` at a quadrature point. This holds for the balance equations of the
nonlocal fields as well: the gradient of a field that does not depend on :math:`\theta` has the components
:math:`(\partial \knl / \partial r, \partial \knl / \partial z)` only, and the weak form

.. math::

   \int_V \Nk_A\, \knl + c\, \Nk_{A,i}\, \knl_{,i} - \Nk_A\, \kl \, \mathrm{d}V, \qquad \mathrm{d}V = 2\pi r\,
   \mathrm{d}A,

is the screened Poisson equation in cylindrical coordinates, :math:`\knl - c \left( \frac{1}{r} \frac{\partial}{\partial
r} \left( r \frac{\partial \knl}{\partial r} \right) + \frac{\partial^2 \knl}{\partial z^2} \right) = \kl`, with the
homogeneous Neumann condition on the boundary and on the axis.

The material is called through the 3D interface with the strain increment
:math:`(\Delta\varepsilon_{rr}, \Delta\varepsilon_{zz}, \Delta\varepsilon_{\theta\theta}, \Delta\gamma_{rz}, 0, 0)`.
The stress, the coupling tangents :math:`\partial \sig / \partial \knl` and :math:`\partial \kl / \partial \eps`, and the
energy tangents are taken in the components (rr, zz, :math:`\theta\theta`, rz); the tangent of the element is the upper
left :math:`4 \times 4` block of the 3D tangent, which is exact for any material, since :math:`\gamma_{r\theta} =
\gamma_{z\theta} = 0` is kinematic.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Element**
     - **Nodes**
     - **Quadrature**
     - **Zero energy modes of** :math:`\mathbf{K}_{uu}`
   * - ``GCAX4``
     - 4
     - Gauss 2x2 (``gaussExactToDegree( 2 )``)
     - 1 (axial translation)
   * - ``GCAX8``
     - 8
     - Gauss 3x3 (``gaussExactToDegree( 4 )``)
     - 1 (axial translation)
   * - ``GCAX8R``
     - 8
     - Gauss 2x2 (``gaussExactToDegree( 3 )``)
     - 2 (axial translation, one hourglass mode)

Notes:

- Internal and external forces, reactions, energies and masses are those of the full :math:`360^\circ` ring. The
  residual of a nonlocal field is integrated over the ring as well, so its nodal values are those of the ring.
- The radial displacement on the axis, :math:`u_r = 0` at :math:`r = 0`, is not imposed by the element; prescribe it.
  The nonlocal fields need no condition on the axis.
- Mesh in :math:`r \geq 0` only. ``initializeYourself()`` rejects a node at :math:`r < -10^{-12} h`, with :math:`h`
  the diagonal of the element's bounding box, so that round-off of nodes on the axis is accepted, and a quadrature
  point at :math:`r \leq 0`.
- A thickness assigned as element property is ignored.
- Mixed-order elements, e.g. ``GeneralGradientEnhancedDisplacementFiniteElement<2, 8, 1, 4>`` with a linear nonlocal
  field, can be instantiated with the axisymmetric section; they are not registered.
- The explicit update and the artificial bulk viscosity are available; the latter acts on
  :math:`\Delta\varepsilon_{rr} + \Delta\varepsilon_{zz} + \Delta\varepsilon_{\theta\theta}`, including the hoop
  strain.
- The energy gradients (:doc:`energygradients`) are available with a material that provides energy tangents.
- Distributed loads: pressure, integrated over the ring.

Evaluation at the quadrature points
-----------------------------------

The implicit and the explicit kernel share one material update and integration per quadrature point for all sections;
the explicit one skips the tangent blocks and the energy gradients and adds the artificial bulk viscosity.

The element stores the shape functions :math:`\mathbf{N}` and their gradients :math:`\partial_{\mathbf{x}} \mathbf{N}` per
quadrature point, but not the strain-displacement matrix :math:`\mathbf{B}`, which is mostly zeros. The implicit kernel forms
:math:`\mathbf{B}` on the fly, because the tangent blocks need it. In the plane and solid sections, neither kernel needs it
for the strain increment and the internal force (the axisymmetric section forms its operator, which carries the hoop row):
they are computed as the symmetric part of the displacement gradient and from the stress directly
(``Marmot::Discretizations::Elements::SmallStrain::strainIncrement`` and ``addInternalForce``),

.. math::

   \Delta \eps = \operatorname{sym} \left( \Delta \mathbf{U}\, \partial_{\mathbf{x}} \mathbf{N}^\mathsf{T} \right), \qquad
   \mathbf{F}_u \mathrel{+}= \left( \sig\, J_0\, w_{qp} \right) \partial_{\mathbf{x}} \mathbf{N}\, ,

with the nodal quantities arranged as :math:`n_\mathrm{dim} \times n_\mathrm{nodes}` matrices :math:`\Delta \mathbf{U}` and
:math:`\mathbf{F}_u`. This is the same as :math:`\mathbf{B}\, \Delta \mathbf{\qu}` and :math:`\mathbf{B}^\mathsf{T} \sig\, J_0\, w_{qp}`,
with a sixth (3D) of the data and half the operations.

In the residual of the nonlocal field, the gradient :math:`\partial_{\mathbf{x}} \knl = \partial_{\mathbf{x}} \Nk\, \mathbf{\qk}` is
formed first, so that the gradient term costs :math:`O(n_\mathrm{nodes})` per point instead of forming the
:math:`n_\mathrm{nodes} \times n_\mathrm{nodes}` matrix :math:`\partial_{\mathbf{x}} \Nk^\mathsf{T} \partial_{\mathbf{x}} \Nk`
(of which the tangent :math:`\mathbf{K}_{\knl\knl}` needs only the term with the gradient coefficient :math:`c`; the term with its derivative is reassociated in the same way). For equal-order elements (as many nonlocal as displacement
nodes), the nonlocal field shares :math:`\mathbf{N}` and :math:`\partial_{\mathbf{x}} \mathbf{N}` with the displacement field and
they are stored only once.

These are reorderings of the same floating-point operations; results change at most in the last digits.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::GeneralGradientEnhancedDisplacementFiniteElement
   :allow-dot-graphs:
