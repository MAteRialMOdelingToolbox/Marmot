.. _embeddedbondelement:

Embedded Bond Element
=====================

Theory
------

Couples an embedded reinforcement bar (a truss or beam element; only its displacements are coupled) to the host
continuum element it lies in by bond-slip. The nodes of the element are the nodes of the bar element, followed by
the nodes of the host element (``Tri3``, ``Quad4``, ``Quad8``, ``Tetra4``, ``Hexa8``, ``Hexa20``: the shapes of the
continuum elements in Marmot). The
slip is the displacement of the bar relative to the host material point it lies at, in the local frame
:math:`\mathbf{R} = [\mathbf{t}, \mathbf{n}_1, \mathbf{n}_2]` of the bar axis,

.. math::

   \mathbf{s} = \mathbf{R}^\mathsf{T} \Big( \sum_i M_i(\eta)\, \mathbf{u}_i
     - \sum_a N_a(\boldsymbol{\xi}(\eta))\, \mathbf{u}_a \Big),

with the host parent coordinates :math:`\boldsymbol{\xi}(\eta)` found by the inverse isoparametric mapping of the host
geometry (``GeometryElement::parentCoordinates``) at initialization. The tangent :math:`\mathbf{t}` points from the
first to the second end node of the bar element (:math:`\mathbf{t} = \mathrm{d}\mathbf{X}/\mathrm{d}\eta`); in 2D, :math:`\mathbf{n}_1` is :math:`\mathbf{t}` rotated by +90°; in 3D, the normals
are an arbitrary right-handed pair that may differ between integration points, so the normal components of slip and
bond stress are not comparable between points (the bond-slip laws are isotropic in the normal plane). A
:ref:`bond-slip law <bondsliplaws>` gives the bond stress :math:`\boldsymbol{\tau}(\mathbf{s})`, integrated over
the bar perimeter :math:`p` and the part :math:`[\eta_s, \eta_e]` of the bar element inside the host element with a
Gauss-Lobatto rule of the quadrature catalogue (``gaussLobattoExactToDegree``), whose end points sit at the ends of
the bar part. The formulation is geometrically linear (reference frame).

Names: ``EB<nDim>D<nBarNodes><host>``, e.g., ``EB2D2T3``, ``EB2D2Q4``, ``EB2D3Q8``, ``EB3D2H8``, ``EB3D3H20``,
``EB3D2T4``.

Element properties: ``[perimeter, eta_s, eta_e, (number of integration points 2..5)]`` with :math:`p > 0` and
:math:`-1 \le \eta_s < \eta_e \le 1`. They must be assigned before the bond-slip law (the number of integration
points sets the number of law instances).

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::EmbeddedBondElement
   :project: Marmot
   :members:
