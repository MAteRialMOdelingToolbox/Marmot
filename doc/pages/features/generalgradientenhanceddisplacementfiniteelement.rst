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

For 2D and 1D problems, the integration measure :math:`J_0 w` is scaled by thickness
:math:`t` or cross-sectional area :math:`A`, respectively: :math:`J_0 w \leftarrow J_0 w\, t`
(2D), :math:`J_0 w \leftarrow J_0 w\, A` (1D).
Here, :math:`t` denotes thickness and :math:`A` denotes cross-sectional area.

Constitutive updates are carried out per quadrature point using a Marmot material model of the following types:

- 3D solid

- Plane stress

- Plane strain

- 1D uniaxial stress

Each quadrature point stores stress, strain and a material state vector; the element
accumulates the history incrementally and may suggest a reduced time step if required
by the material.


Implementation
--------------

.. doxygenclass:: Marmot::Elements::GeneralGradientEnhancedDisplacementFiniteElement
   :allow-dot-graphs:
