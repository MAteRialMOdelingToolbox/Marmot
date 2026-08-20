.. _incompressibleneohooke:

Incompressible Neo-Hooke model
===============================

Theory
------

The classical incompressible (purely isochoric) Neo-Hookean model has energy density

.. math::

   \Psi_{\rm iso}(\bar\lambda_1,\bar\lambda_2,\bar\lambda_3)
   =
   \frac{\mu}{2} \left( \bar\lambda_1^2 + \bar\lambda_2^2 + \bar\lambda_3^2 - 3 \right)
   =
   \frac{\mu}{2} \left( \bar I_1 - 3 \right),

where :math:`\bar\lambda_i = J^{-1/3}\lambda_i` are the isochoric principal stretches
(:math:`J=\det\boldsymbol F`) and :math:`\bar I_1` is the first invariant of the isochoric right
Cauchy-Green tensor. Unlike :ref:`compressibleneohooke`, there is no volumetric energy term
:math:`U(J)`: the material resists no volumetric deformation on its own and must be paired with a
mixed displacement-pressure element, e.g. :ref:`displacementpressurefinitestrainelement`, that
enforces incompressibility as an independent constraint.

This is exactly the one-term :ref:`ogden` model with exponent :math:`\alpha=2` and modulus
:math:`\mu_1=\mu`. Unlike Ogden's general (non-integer) exponents, however, this potential is
quadratic in the isochoric stretches and can therefore be expressed directly in terms of the
invariants of :math:`\boldsymbol C=\boldsymbol F^{\mathsf T}\boldsymbol F` (:math:`\bar I_1 = I_1
J^{-2/3}`, :math:`I_1=\operatorname{tr}\boldsymbol C`) exactly as the compressible potentials in
:ref:`compressibleneohooke` are - so, unlike Ogden, **no spectral decomposition of**
:math:`\boldsymbol C` **is required or performed**. The second Piola-Kirchhoff stress
:math:`\boldsymbol S = 2\,\partial\Psi_{\rm iso}/\partial\boldsymbol C` is pushed forward to the
Kirchhoff stress exactly as for :ref:`compressibleneohooke`. The consistent algorithmic tangent is
computed analytically (no automatic differentiation), reusing the same
:math:`\partial J/\partial\boldsymbol C`, :math:`\partial^2 J/\partial\boldsymbol C\partial
\boldsymbol C` and :math:`\partial I_1/\partial\boldsymbol C` building blocks already present in
:ref:`compressibleneohooke`'s Pence-Gou potential.

Implementation
--------------

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`\mu`
     - Shear modulus
   * - 1 (optional)
     - :math:`\rho`
     - Density

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**

   * - ---
     - ---
     - No state variables required.

.. doxygenclass:: Marmot::Materials::IncompressibleNeoHooke
   :allow-dot-graphs:
