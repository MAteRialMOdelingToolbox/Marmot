Displacement Finite Strain Element
==================================

Preliminaries
------------- 
This element formulation is implemented in the class :cpp:class:`Marmot::Elements::DisplacementFiniteStrainULElement`, which represents a displacement-based finite element considering large deformations. The formulation is valid for general three-dimensional geometries and can be combined with arbitrary material models. It supports both body forces and surface tractions.

This class inherits from the abstract base classes :cpp:class:`MarmotElement` and :cpp:class:`MarmotGeometryElement<nDim,nNodes>`, which provide common functionality for all finite elements in **Marmot**. In addition, it serves as a base class for specializations such as plane strain and axisymmetric elements.

Throughout the following derivations, Einstein summation convention is used, i.e., repeated indices imply summation. Uppercase indices refer to material coordinates, while lowercase indices refer to spatial coordinates. Accordingly, lowercase subscripts preceded by a comma denote spatial derivatives with respect to the current (deformed) configuration, e.g., :math:`(\bullet)_{,i} = \partial (\bullet)/\partial x_i`, while uppercase subscripts preceded by a comma denote material derivatives with respect to the reference (undeformed) configuration, e.g., :math:`(\bullet)_{,I} = \partial (\bullet)/\partial X_I`. Following this convention, the deformation gradient and its determinant is defined as

.. math::
  F_{iI} = x_{i,I}\ ,\ J = \text{det}\,F_{iI}.

Theory
------
This element computes the element residual vector :math:`\mathbf{r}_{Aj}` and its derivative with respect to the nodal displacement vector :math:`\mathbf{q}_{Bk}`, i.e., the element stiffness matrix :math:`\partial \mathbf{r}_{Aj}/\partial \mathbf{q}_{Bk}` for a displacement-based finite element considering large deformations. The residual vector is derived by discretizing the weak form for linear momentum, given as

.. math::
   \mathbf{r}_{Aj} = \int_{V_0}\,\mathbf{N}_{A,i}\,\tau_{ij}\,dV_0 - \int_{V_0}\,\mathbf{N}_{A}\,f_{j}\,dV_0 - \int_{\bar{A}}\, \mathbf{N}_A\, \bar{t}_j \, d\bar{A},

with the Kirchhoff stress :math:`\tau_{ij} = J\,t_{ij}`, the body force :math:`f_j`, and the prescribed traction :math:`\bar{t}_j` on the Neumann boundary :math:`\bar{A}`. The element stiffness matrix is computed as

.. math::

    \frac{\partial \mathbf{r}_{Aj}}{\partial \mathbf{q}_{Bk}} &= \mathbf{K}_{AjBk} = \int_{V_0} \mathbf{N}_{A,i}\, \frac{\partial \tau_{ij}}{\partial F_{kK}}\,\mathbf{N}_{B,K} - \mathbf{N}_{A,k}\,\mathbf{N}_{B,i}\,\tau_{ij}\,dV_0\, -\\
   &- \int_{\bar{A}}\,\mathbf{N}_A\,\bar{t}_i\, \left(\delta_{ij}\delta_{lk} - \delta_{ik}\delta_{lj}\right)\, \mathbf{N}_{B,l}\, d\bar{A},

where :math:`\delta_{ij}` is the Kronecker delta and :math:`\partial \tau_{ij}/\partial F_{kK}` the material tangent obtained from the respective material model.

Implementation
--------------

.. doxygenclass:: Marmot::Elements::DisplacementFiniteStrainULElement
   :allow-dot-graphs:
