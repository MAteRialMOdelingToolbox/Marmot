Displacement Finite Strain F-bar Element
=========================================

Preliminaries
-------------
This element formulation is implemented in the class :cpp:class:`Marmot::Elements::DisplacementFiniteStrainFBarElement`, which derives from :cpp:class:`Marmot::Elements::DisplacementFiniteStrainULElement` and augments it with the F-bar method of de Souza Neto, Perić, Dutko and Owen (1996) for locking-free analysis of low-order elements (Tetra4, Hexa8) with nearly incompressible materials. It is a locking *mitigation* technique for materials that retain a (possibly very large) volumetric stiffness :math:`U(J)`, such as :cpp:class:`Marmot::Materials::CompressibleNeoHooke` with a large bulk-to-shear modulus ratio; it is not applicable to genuinely incompressible materials without any volumetric energy term (e.g. :cpp:class:`Marmot::Materials::Ogden`, :cpp:class:`Marmot::Materials::IncompressibleNeoHooke`, :cpp:class:`Marmot::Materials::IncompressibleMooneyRivlin`), which should instead be paired with :cpp:class:`Marmot::Elements::DisplacementPressureFiniteStrainElement`.

The element introduces no new degrees of freedom: node fields, DOF layout, quadrature, distributed loads and inertia are all inherited unchanged from :cpp:class:`Marmot::Elements::DisplacementFiniteStrainULElement`. Only the residual/tangent assembly (``computeKernels``/``computeKernelsExplicit``) is overridden.

Theory
------
At each quadrature point, the deformation gradient :math:`F_{iI}` and its volume ratio :math:`J = \text{det}\,F_{iI}` are computed as usual. In addition, a deformation gradient :math:`F^0_{iI}` and volume ratio :math:`J_0 = \text{det}\,F^0_{iI}` are sampled once at the element centroid, shared by every quadrature point. The material is evaluated not at :math:`F_{iI}` itself but at the F-bar modified deformation gradient

.. math::
  \bar{F}_{iI} = \left(\frac{J_0}{J}\right)^{1/3} F_{iI},

which decouples the volumetric response (effectively sampled once per element, via :math:`J_0`) from the deviatoric response (still sampled at every quadrature point). This is purely a substitution at the constitutive level: the residual retains the usual form, evaluated with the modified stress :math:`\bar{\tau}_{ij}` but the *actual* spatial gradient operator :math:`\mathbf{N}_{A,i}` (built from :math:`F_{iI}`, not :math:`\bar{F}_{iI}`),

.. math::
  \mathbf{r}_{Aj} = \int_{V_0}\,\mathbf{N}_{A,i}\,\bar{\tau}_{ij}\,dV_0.

Because :math:`J_0` depends on every node's degrees of freedom (not just those local to a quadrature point), linearizing :math:`\bar{F}_{iI}` with respect to the nodal displacements :math:`\mathbf{q}_{Bk}` produces, writing :math:`\alpha = (J_0/J)^{1/3}`,

.. math::
  \frac{\partial \bar{F}_{iI}}{\partial \mathbf{q}_{Bk}} = \alpha\,\frac{\partial F_{iI}}{\partial \mathbf{q}_{Bk}} + F_{iI}\,\frac{\partial \alpha}{\partial \mathbf{q}_{Bk}}\ ,\qquad
  \frac{\partial \alpha}{\partial \mathbf{q}_{Bk}} = \frac{\alpha}{3}\,\left(\mathbf{N}^0_{B,k} - \mathbf{N}_{B,k}\right),

using :math:`\partial J/\partial \mathbf{q}_{Bk} = J\,\mathbf{N}_{B,k}` at both the quadrature point and the centroid. Contracting through the material tangent :math:`a_{ijkl} = \partial\bar\tau_{ij}/\partial\bar F_{kl}` gives three contributions to the element stiffness matrix per quadrature point: the standard geometric term (built from :math:`\bar\tau_{ij}` exactly as in :cpp:class:`Marmot::Elements::DisplacementFiniteStrainULElement`), the standard material term scaled by :math:`\alpha`, and an additional F-bar correction term

.. math::
  \frac{\alpha}{3}\,\left(a_{ijkl}\,F_{kl}\right)\,\mathbf{N}_{A,i}\,\left(\mathbf{N}^0_{B,n} - \mathbf{N}_{B,n}\right)

coupling every quadrature point's stiffness contribution to every node's degrees of freedom through the centroid's sensitivity :math:`\mathbf{N}^0_{B,n}` - the term a naive "substitute :math:`\bar F` for :math:`F`" implementation would miss.

For an element with a single quadrature point coincident with the centroid (e.g. a fully-integrated Tetra4), :math:`F_{iI} \equiv F^0_{iI}` identically, so :math:`\alpha \equiv 1` and :math:`\mathbf{N}^0_{B,n} \equiv \mathbf{N}_{B,n}`: the correction term vanishes and the element reduces exactly to :cpp:class:`Marmot::Elements::DisplacementFiniteStrainULElement`.

Implementation
--------------

.. doxygenclass:: Marmot::Elements::DisplacementFiniteStrainFBarElement
   :allow-dot-graphs:
