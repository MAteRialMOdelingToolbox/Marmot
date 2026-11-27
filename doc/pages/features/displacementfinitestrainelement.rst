Displacement Finite Strain Element
==================================

Preliminaries
-------------
This element formulation is implemented in the class :cpp:class:`Marmot::Discretizations::Elements::DisplacementFiniteStrainULElement`, which represents a displacement-based finite element considering large deformations. The formulation is valid for general three-dimensional geometries and can be combined with arbitrary material models. It supports both body forces and surface tractions.

This class inherits from the abstract base classes :cpp:class:`Marmot::Discretizations::Elements::Element` and :cpp:class:`Marmot::Discretizations::Elements::GeometryElement`, which provide common functionality for all finite elements in **Marmot**. In addition, it serves as a base class for specializations such as plane strain and axisymmetric elements.

Throughout the following derivations, Einstein summation convention is used, i.e., repeated indices imply summation. Uppercase indices refer to material coordinates, while lowercase indices refer to spatial coordinates. Accordingly, lowercase subscripts preceded by a comma denote spatial derivatives with respect to the current (deformed) configuration, e.g., :math:`(\bullet)_{,i} = \partial (\bullet)/\partial x_i`, while uppercase subscripts preceded by a comma denote material derivatives with respect to the reference (undeformed) configuration, e.g., :math:`(\bullet)_{,I} = \partial (\bullet)/\partial X_I`. Following this convention, the deformation gradient and its determinant is defined as

.. math::
  F_{iI} = x_{i,I}\ ,\ J = \text{det}\,F_{iI}.

Theory
------
This element computes the element contribution to the global residual vector :math:`\mathbf{r}_{Aj}` and its derivative with respect to the nodal displacement vector :math:`\mathbf{q}_{Bk}`, i.e., the element contribution to the global stiffness matrix :math:`\partial \mathbf{r}_{Aj}/\partial \mathbf{q}_{Bk}` for a displacement-based finite element considering large deformations. The global residual vector is derived by discretizing the weak form for linear momentum, given as

.. math::
   \mathbf{r}_{Aj} = \int_{V_0}\,\mathbf{N}_{A,i}\,\tau_{ij}\,dV_0 - \int_{V_0}\,\mathbf{N}_{A}\,f_{j}\,dV_0 - \int_{\bar{A}}\, \mathbf{N}_A\, \bar{t}_j \, d\bar{A},

with the Kirchhoff stress :math:`\tau_{ij} = J\,t_{ij}`, the body force :math:`f_j`, and the prescribed traction :math:`\bar{t}_j` on the Neumann boundary :math:`\bar{A}`. The global stiffness matrix is computed as

.. math::

    \frac{\partial \mathbf{r}_{Aj}}{\partial \mathbf{q}_{Bk}} &= \mathbf{K}_{AjBk} = \int_{V_0} \mathbf{N}_{A,i}\, \frac{\partial \tau_{ij}}{\partial F_{kK}}\,\mathbf{N}_{B,K} - \mathbf{N}_{A,k}\,\mathbf{N}_{B,i}\,\tau_{ij}\,dV_0\, -\\
   &- \int_{\bar{A}}\,\mathbf{N}_A\,\bar{t}_i\, \left(\delta_{ij}\delta_{lk} - \delta_{ik}\delta_{lj}\right)\, \mathbf{N}_{B,l}\, d\bar{A},

where :math:`\delta_{ij}` is the Kronecker delta and :math:`\partial \tau_{ij}/\partial F_{kK}` the material tangent obtained from the respective material model.

Axisymmetric elements
---------------------
``CX8UL`` and ``CX8RUL`` (:cpp:class:`Marmot::Discretizations::Elements::AxiSymmetricDisplacementFiniteStrainULElement`)
are the 8-node plane element in the meridian half plane, :math:`r = x_1 \geq 0` and :math:`z = x_2`, with the hoop
stretch :math:`F_{\theta\theta} = 1 + u_r/r` as the third diagonal entry of the deformation gradient. Everything the
element integrates is integrated over the ring, :math:`\mathrm{d}V = 2\pi r\,\mathrm{d}A`: the internal force and
its tangent, the lumped and the consistent mass, the body force, and the operators of the critical time step. The
ring weight takes the place of the thickness of a plane element; the thickness property is not used. The distributed
loads are integrated over the ring of the face, :math:`\mathrm{d}\bar{A} = 2\pi r\,\mathrm{d}s`: a pressure on the
current face, with the current radius, whose derivative is part of the load stiffness, and a surface traction on the
reference face.

The geometry is checked once, when the coordinates are assigned: a node at :math:`r < 0` (beyond the round-off of a
mesher) or a quadrature point at :math:`r \leq 0`, where the hoop row :math:`N_i/r` is undefined and the ring weight
not positive, is refused with ``std::invalid_argument``. The hoop stretch is checked in every increment: a ring pulled
through the axis, :math:`F_{\theta\theta} \leq 0`, throws ``std::runtime_error``, which a solver can answer with a
smaller increment. The kernels and the wave speed of the critical time step take the same :math:`F_{\theta\theta}`.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::DisplacementFiniteStrainULElement
   :allow-dot-graphs:
