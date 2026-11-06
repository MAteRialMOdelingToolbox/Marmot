.. _embeddedbondelement:

Embedded Bond Element
=====================

Theory
------

Couples an embedded reinforcement bar (a truss or beam element; only its displacements are coupled) to the host
continuum element it lies in by bond-slip. The nodes of the element are the nodes of the bar element, followed by the nodes of the host
element (``Quad4``, ``Quad8``, ``Tetra4``, ``Tetra10``, ``Hexa8``, ``Hexa20``). The slip is the displacement of the
bar relative to the host material point it lies at, in the local frame
:math:`\mathbf{R} = [\mathbf{t}, \mathbf{n}_1, \mathbf{n}_2]` of the bar axis,

.. math::

   \mathbf{s} = \mathbf{R}^\mathsf{T} \Big( \sum_i M_i(\eta)\, \mathbf{u}_i
     - \sum_a N_a(\boldsymbol{\xi}(\eta))\, \mathbf{u}_a \Big),

with the host parent coordinates :math:`\boldsymbol{\xi}(\eta)` found by inverse mapping at initialization. A
:ref:`bond-slip law <bondsliplaws>` gives the bond stress :math:`\boldsymbol{\tau}(\mathbf{s})`, integrated over
the bar perimeter :math:`p` and the part :math:`[\eta_s, \eta_e]` of the bar element inside the host element with a
Gauss-Lobatto rule. The formulation is geometrically linear (reference frame).

Names: ``EB<nDim>D<nBarNodes><host>``, e.g., ``EB2D2Q4``, ``EB2D3Q8``, ``EB3D2H8``, ``EB3D3H20``, ``EB3D2T4``.

Element properties: ``[perimeter, eta_s, eta_e, (number of integration points 2..5)]``.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::EmbeddedBondElement
   :project: Marmot
   :members:
