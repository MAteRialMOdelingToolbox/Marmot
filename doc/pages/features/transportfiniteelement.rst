.. _transportfiniteelement:

Transport Finite Element (heat conduction, moisture diffusion, coupled heat and moisture transport)
===================================================================================================

Theory
------

``TransportFiniteElement< nDim, nNodes, nFields >`` integrates the balance of :math:`n_f` scalar fields
:math:`\phi_a` (the temperature for heat conduction, the pore relative humidity for moisture diffusion, both for
coupled heat and moisture transport),

.. math::

   \dot m_a + \nabla\cdot\mathbf{q}_a = r_a + s_a ,

with the content :math:`m_a` (enthalpy per unit volume, or water content), the flux :math:`\mathbf{q}_a` and the material source
:math:`r_a` of a transport material (:cpp:class:`Marmot::Materials::TransportMaterial`), and the external body source
:math:`s_a`. On the boundary, the field is prescribed or the outward normal flux :math:`q_{n,a}` is given by a surface
load. The weak form is

.. math::

   \int_\Omega \delta\phi_a\,\dot m_a\,\mathrm{d}V - \int_\Omega \nabla\delta\phi_a\cdot\mathbf{q}_a\,\mathrm{d}V
   - \int_\Omega \delta\phi_a\,r_a\,\mathrm{d}V
   = \int_\Omega \delta\phi_a\,s_a\,\mathrm{d}V - \int_\Gamma \delta\phi_a\,q_{n,a}\,\mathrm{d}\Gamma .

Time integration is backward Euler in conservative form,
:math:`\dot m_a \approx (m_a(\boldsymbol\phi_{n+1}) - m_a(\boldsymbol\phi_n))/\Delta t`, which conserves the content also for
a field dependent capacity. With the interpolation operator :math:`\mathbf{N}_\phi` and the gradient operator
:math:`\mathbf{B}_\phi` (three rows per field, padded with zeros in 2D), the residual and its consistent tangent are

.. math::

   \mathbf{P}_e = \sum_{qp} \Big[ \mathbf{N}_\phi^\mathsf{T}\big(\tfrac{\Delta\mathbf{m}}{\Delta t} - \mathbf{r}\big)
   - \mathbf{B}_\phi^\mathsf{T}\mathbf{q} \Big] J_0 w, \qquad
   \mathbf{K}_e = \sum_{qp} \Big[ \mathbf{N}_\phi^\mathsf{T}\big(\tfrac{1}{\Delta t}
   \tfrac{\partial\Delta\mathbf{m}}{\partial\boldsymbol\phi} - \tfrac{\partial\mathbf{r}}{\partial\boldsymbol\phi}\big)
   \mathbf{N}_\phi - \mathbf{B}_\phi^\mathsf{T}\big(\tfrac{\partial\mathbf{q}}{\partial\boldsymbol\phi}\mathbf{N}_\phi
   + \tfrac{\partial\mathbf{q}}{\partial\nabla\boldsymbol\phi}\mathbf{B}_\phi\big) \Big] J_0 w,

with :math:`J_0 w` including the thickness (plane), :math:`2\pi r` (axisymmetric) or 1 (solid). A zero time increment
lets no time pass: the element balances the content increment only (:math:`\mathbf{P}_e = \sum
\mathbf{N}_\phi^\mathsf{T}\Delta\mathbf{m}\,J_0 w`) and the loads vanish, unless the material has no storage
(zero capacity), which gives the steady state equations for any time increment. The distinction is made field by field
(``TransportMaterial::hasStorage( field )``).

The dofs are ordered node by node, the fields of a node in the order of the element's field names.

Elements
--------

Named after the diffusive heat transfer elements of Abaqus, with a suffix for the nodal fields (``T``: temperature,
``H``: relative humidity, ``TH``: both). The same eight shapes are registered for each suffix: ``DC2D3H``, ``DC2D4H``,
``DC2D8H``, ``DCAX4H``, ``DCAX8H``, ``DC3D4H``, ``DC3D8H``, ``DC3D20H`` carry the field ``relative humidity``;
``DC2D3TH`` ... ``DC3D20TH`` (``nFields = 2``) carry ``temperature`` and ``relative humidity``, in this order at each
node, and take a two-field material such as :ref:`THERMOHYGRO <thermohygro>`.

=========== ======== ====== ============== ================
Name        Section  Shape  Quadrature     Fields
=========== ======== ====== ============== ================
``DC2D3T``  plane    Tri3   3 points       ``temperature``
``DC2D4T``  plane    Quad4  2 x 2          ``temperature``
``DC2D8T``  plane    Quad8  3 x 3          ``temperature``
``DCAX4T``  axisym.  Quad4  2 x 2          ``temperature``
``DCAX8T``  axisym.  Quad8  3 x 3          ``temperature``
``DC3D4T``  solid    Tetra4 4 points       ``temperature``
``DC3D8T``  solid    Hexa8  2 x 2 x 2      ``temperature``
``DC3D20T`` solid    Hexa20 3 x 3 x 3      ``temperature``
=========== ======== ====== ============== ================

Element properties: the thickness of plane elements. Axisymmetric elements use :math:`r = x_1 \geq 0`,
:math:`z = x_2`, and integrate over the full ring.

Loads
-----

Distributed loads act on the field of their physics (load values in brackets); plane elements multiply by the
thickness, axisymmetric ones integrate over the ring. On ``temperature``:

- ``HeatFlux`` [:math:`\bar q`]: heat flux into the body per unit area, :math:`q_n = -\bar q`;
- ``HeatFilm`` [:math:`h, T_\infty`]: convection, :math:`q_n = h\,(T - T_\infty)`;
- ``HeatRadiation`` [:math:`c, T_\infty, T_0`]: radiation, :math:`q_n = c\,((T-T_0)^4 - (T_\infty-T_0)^4)`, with
  :math:`c = \varepsilon\sigma` and :math:`T_0` the absolute zero of the temperature scale (0 for K, -273.15 for °C).

On ``relative humidity``:

- ``MoistureFlux`` [:math:`\bar q`]: moisture flux into the body per unit area, :math:`q_n = -\bar q`;
- ``MoistureFilm`` [:math:`\beta, h_\infty`]: exchange with the environment, :math:`q_n = \beta\,(h - h_\infty)`.

On ``temperature``, depending on ``relative humidity`` (``DC*TH`` only):

- ``LatentHeatFilm`` [:math:`c, h_\infty`]: the latent heat taken away by the evaporating water,
  :math:`q_n = c\,(h - h_\infty)`, :math:`c = L_v\beta` with the coefficient :math:`\beta` of the moisture film. Its
  tangent couples the temperature rows to the relative humidity.

Films and radiation depend on the field and add their tangent. The body force is a volumetric source, one value per
field.

States per quadrature point: ``field``, ``field gradient``, ``flux`` (three components per field), and the states of
the material.

.. doxygenclass:: Marmot::Discretizations::Elements::TransportFiniteElement
   :project: Marmot
   :members:
