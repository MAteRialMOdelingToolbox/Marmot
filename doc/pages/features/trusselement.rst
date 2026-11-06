.. _trusselement:

Truss Element
=============

Theory
------

A truss (bar) element in two or three spatial dimensions with 2 (linear) or 3 (quadratic; node order end, end,
mid) nodes. It carries only an axial force. The material sees a uniaxial stress state along the truss axis, found by
the uniaxial stress reductions of the material base classes, so any 3D material of Marmot can be used, e.g., for
reinforcement bars embedded in continuum elements. The names follow the Marmot (Abaqus-like) scheme
``T<dimension>D<nodes>``; the suffix ``CR`` marks the co-rotational finite strain variants.

**Small strain** (``T2D2``, ``T2D3``, ``T3D2``, ``T3D3``; hypoelastic materials): the axial strain is
:math:`\varepsilon = \mathbf{T} \cdot \partial \mathbf{u} / \partial S` with the reference axis :math:`\mathbf{T}`
and arc length :math:`S`. The material is called through
``Marmot::Materials::MaterialHypoElastic::computeUniaxialStress``.

**Finite strain** (``T2D2CR``, ``T2D3CR``, ``T3D2CR``, ``T3D3CR``; finite strain materials): the material is
evaluated in the co-rotated frame of the truss, :math:`\mathbf{F} = \mathrm{diag}(\lambda, \lambda_2, \lambda_3)`,
with the axial stretch :math:`\lambda = |\mathbf{g}| / |\mathbf{G}|` and the lateral stretches found by
``Marmot::Materials::MaterialFiniteStrain::computeUniaxialStress`` such that :math:`\tau_{22} = \tau_{33} = 0`.
Rigid body rotations are exactly stress free. With the nominal stress :math:`P = \tau_{11} / \lambda` and the current axis
:math:`\mathbf{n}`,

.. math::

   \mathbf{f}_a = \int P\, \mathbf{n}\, \frac{\partial N_a}{\partial S}\, A_0\, \mathrm{d}S, \qquad
   \mathbf{K}_{ab} = \int \Big[ \frac{\mathrm{d}P}{\mathrm{d}\lambda}\, \mathbf{n} \otimes \mathbf{n}
     + \frac{P}{\lambda} (\mathbf{I} - \mathbf{n} \otimes \mathbf{n}) \Big]
     \frac{\partial N_a}{\partial S} \frac{\partial N_b}{\partial S}\, A_0\, \mathrm{d}S .

Element property: the cross section area :math:`A_0`.

Quadrature point states: ``stress`` (axial Cauchy stress), ``strain`` (small strain: axial strain; finite strain:
:math:`\ln\lambda`), ``normal force``, ``elastic energy``, ``dissipation``, ``lateral stretches``, ``kirchhoff
stress``, followed by the material state.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::TrussElement
   :project: Marmot
   :members:
