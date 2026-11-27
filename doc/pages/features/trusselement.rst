.. _trusselement:

Truss Element
=============

Theory
------

A truss (bar) element in two or three spatial dimensions with 2 (linear) or 3 (quadratic; node order end, end,
mid) nodes. It carries only an axial force, e.g., as a reinforcement bar embedded in continuum elements. The material
sees a uniaxial stress state along the truss axis, found by the uniaxial stress reductions of the material base
classes. The names follow Marmot's scheme ``T<dimension>D<nodes>`` (Abaqus: ``T2D2``, ``T3D2`` etc. as well); the
suffix ``CR`` marks the finite strain variants.

**Element and material go in pairs:**

.. list-table::
   :header-rows: 1

   * - element
     - kinematics
     - materials
   * - ``T2D2``, ``T2D3``, ``T3D2``, ``T3D3``
     - geometrically linear
     - hypoelastic, e.g., ``LINEARELASTIC``, ``VONMISES``
   * - ``T2D2CR``, ``T2D3CR``, ``T3D2CR``, ``T3D3CR``
     - finite strain
     - finite strain materials whose response to a diagonal deformation gradient is diagonal, e.g., the isotropic
       ``COMPRESSIBLEHYPERELASTICITY``, ``FINITESTRAINJ2PLASTICITY``, ``VONMISES/HUGHES-WINGET``

A material of the other kind is rejected with a hint to the matching element (``VONMISES`` is a hypoelastic material;
its finite strain counterpart is ``VONMISES/HUGHES-WINGET``). The uniaxial reduction of finite strain materials
condenses only :math:`\tau_{22}` and :math:`\tau_{33}`; a material that answers a diagonal deformation gradient with
shear stresses (e.g., an anisotropic one, its axes read in the truss frame) is rejected instead of having its shear
stresses dropped.

**Small strain** (``T2D2``, ``T2D3``, ``T3D2``, ``T3D3``): the axial strain is :math:`\varepsilon = \mathbf{T} \cdot
\partial \mathbf{u} / \partial S` with the reference unit tangent :math:`\mathbf{T}` and arc length :math:`S`. The
material is called through ``Marmot::Materials::MaterialHypoElastic::computeUniaxialStress``. As any geometrically
linear element -- and unlike ``T2D2`` in Abaqus with NLGEOM -- the element strains under a finite rotation: use the
``CR`` variant whenever rotations are not small.

**Formulation of the finite strain variants** (``T2D2CR``, ``T2D3CR``, ``T3D2CR``, ``T3D3CR``): total Lagrangian; at
each quadrature point the material is evaluated in the frame co-rotated with the bar's current tangent
:math:`\mathbf{g} = \partial \mathbf{x} / \partial S`, :math:`\mathbf{F} = \mathrm{diag}(\lambda, \lambda_2,
\lambda_3)` with the axial stretch :math:`\lambda = |\mathbf{g}|` and the lateral stretches found by
``Marmot::Materials::MaterialFiniteStrain::computeUniaxialStress`` such that :math:`\tau_{22} = \tau_{33} = 0`.
Rigid body rotations are exactly stress free. For 2 nodes the frame is the chord, and the element is the
co-rotational truss of Crisfield (*Non-linear Finite Element Analysis of Solids and Structures*, Vol. 1, ch. 3),
hence ``CR``; for 3 nodes the frame follows the curved bar point by point. With the first Piola-Kirchhoff stress
:math:`P = \tau_{11} / \lambda` and the current axis :math:`\mathbf{n}`,

.. math::

   \mathbf{f}_a = \int P\, \mathbf{n}\, \frac{\partial N_a}{\partial S}\, A_0\, \mathrm{d}S, \qquad
   \mathbf{K}_{ab} = \int \Big[ \underbrace{\frac{\mathrm{d}P}{\mathrm{d}\lambda}\, \mathbf{n} \otimes
     \mathbf{n}}_{\text{material}}
     + \underbrace{\frac{P}{\lambda} (\mathbf{I} - \mathbf{n} \otimes \mathbf{n})}_{\text{geometric}} \Big]
     \frac{\partial N_a}{\partial S} \frac{\partial N_b}{\partial S}\, A_0\, \mathrm{d}S :

the material stiffness along the bar and the geometric (initial stress) stiffness transverse to it, with
:math:`P / \lambda = S_{11}`, the second Piola-Kirchhoff stress -- the term behind the snap-through of a shallow
truss.

Element property: the cross section area :math:`A_0` (in the units of the model, e.g., mm² for N, mm, MPa).

Quadrature: Gauss, exact to degree :math:`2p - 1` (1 point for 2 nodes, 2 points for 3 nodes); exact for the stiffness
of straight trusses with a linear elastic material, approximate for curved 3-node bars.

Quadrature point states: ``stress`` (axial Cauchy stress), ``strain`` (small strain: axial strain; finite strain:
:math:`\ln\lambda`), ``normal force``, ``elastic energy``, ``dissipation``; for the ``CR`` variants also ``lateral
stretches`` and ``kirchhoff stress``; followed by the material state.

Limitations
-----------

* A truss has no lateral stiffness. A node connected only to collinear small-strain trusses (or to ``CR`` trusses
  without tension) is free laterally, and the stiffness matrix is singular -- the classic trap of a first truss model.
  Restrain such nodes, e.g., by boundary conditions, by embedding the truss, or by pre-tension.
* A quadratic truss (``T2D3``, ``T3D3`` and their ``CR`` variants) has internal zero-energy modes of its mid node
  transverse to the axis: one in 2D, two in 3D. Harmless when the truss is embedded in a continuum (the main use);
  stand-alone, the mid node must be restrained transversally (constraint, boundary condition) or pre-tensioned.
* Only axial forces: no bending, no distributed (surface) loads.

How to check
------------

A hand calculation for a first model: an inclined linear elastic bar of length :math:`L` and area :math:`A`, one
end fixed, the other pulled along the bar by :math:`\Delta`, carries the normal force :math:`N = E A \Delta / L` and
the end forces :math:`\pm N \mathbf{T}`; the stiffness along the bar is :math:`EA/L`, the one transverse to it zero
(small strain) or :math:`N / L` (``CR``, geometric stiffness). The EdelweissFE test
``tests/test_trusses.py::test_truss_small_strain_inclined`` is this check as a complete model.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::TrussElement
   :project: Marmot
   :members:
