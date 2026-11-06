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

Large slip
----------

``EmbeddedLargeSlipBondElement`` lets the bar slide out of (or into) its host element. The bond is integrated over
the **channel**, the part of the original bar path inside the host element, fixed in the host. The bond history
belongs to the channel points (i.e., to the concrete). At every evaluation, the current partner of each channel point
on the bar is found by closest point projection onto the current bar axis, searched in a window of consecutive bar
elements. The tangential slip is the reference arc length the bar has slid past the channel point,

.. math::

   s_t = S_c - S^*,

with the reference arc lengths :math:`S_c` of the channel point and :math:`S^*` of its bar partner; the normal slip
is the relative position of the partner, normal to the axis. The bond of a channel point is scaled by
:math:`f = 3x^2 - 2x^3` of the fraction :math:`x \in [0, 1]` of its tributary length still covered by the bar, so the
bond fades smoothly as a bar end passes it; a channel point the bar has left has no bond and a frozen history. The
bonded length therefore shrinks exactly as the bar is pulled out. For small slips, the element is identical to the
small-slip element. The geometry is the reference one; the tangent is computed by central differences of the
internal forces.

Nodes: the host element nodes, followed by the unique nodes of the window of bar elements in chain order (bar2:
``e0, e1, ..., eN``; bar3: ``e0, m0, e1, m1, ..., eN``).

Names: ``EBLS<nDim>D<nBarNodes><host>W<nBarElements>``, with 1 to 48 bar elements in the window, e.g.,
``EBLS2D2Q4W3``, ``EBLS3D3H20W12``.

Element properties: ``[perimeter, c_s, c_e, start is a bar end, end is a bar end, (number of integration points)]``,
with the channel :math:`[c_s, c_e]` in the chain parameter :math:`c = k + (\xi + 1)/2 \in [0, n_\mathrm{bar\,elements}]`
of the window. If a partner lies beyond an end of the window that is not an end of the bar, the window is too small
for the slip and an exception is thrown.

Additional states: ``partner`` (chain parameter of the bar partner), ``active`` (1 while the channel point is bonded)
and ``covered fraction`` (:math:`x`).

.. doxygenclass:: Marmot::Discretizations::Elements::EmbeddedLargeSlipBondElement
   :project: Marmot
   :members:
