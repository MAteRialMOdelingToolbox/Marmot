.. _incompressiblemooneyrivlin:

Incompressible Mooney-Rivlin model
====================================

Theory
------

The classical incompressible (purely isochoric) Mooney-Rivlin model has energy density

.. math::

   \Psi_{\rm iso}(\bar\lambda_1,\bar\lambda_2,\bar\lambda_3)
   =
   C_1 \left( \bar I_1 - 3 \right) + C_2 \left( \bar I_2 - 3 \right),

where :math:`\bar\lambda_i = J^{-1/3}\lambda_i` are the isochoric principal stretches
(:math:`J=\det\boldsymbol F`),
:math:`\bar I_1 = \bar\lambda_1^2+\bar\lambda_2^2+\bar\lambda_3^2`, and, using
:math:`\bar\lambda_1\bar\lambda_2\bar\lambda_3=1`,
:math:`\bar I_2 = \bar\lambda_1^{-2}+\bar\lambda_2^{-2}+\bar\lambda_3^{-2}`. Unlike the compressible
potentials in :ref:`compressibleneohooke`, there is no volumetric energy term :math:`U(J)`: the
material resists no volumetric deformation on its own and must be paired with a mixed
displacement-pressure element, e.g. :ref:`displacementpressurefinitestrainelement`, that enforces
incompressibility as an independent constraint.

This is exactly the two-term :ref:`ogden` model with :math:`(\mu_1,\alpha_1)=(2C_1,2)` and
:math:`(\mu_2,\alpha_2)=(-2C_2,-2)`. Unlike Ogden's general (non-integer) exponents, however, this
potential is quadratic in the isochoric stretches and can therefore be expressed directly in terms
of the invariants of :math:`\boldsymbol C=\boldsymbol F^{\mathsf T}\boldsymbol F` exactly as the
compressible Mooney-Rivlin potential is - so, unlike Ogden, **no spectral decomposition of**
:math:`\boldsymbol C` **is required or performed**. The second Piola-Kirchhoff stress
:math:`\boldsymbol S = 2\,\partial\Psi_{\rm iso}/\partial\boldsymbol C` is pushed forward to the
Kirchhoff stress. The consistent algorithmic tangent is computed analytically (no automatic
differentiation), extending the same :math:`J`/:math:`I_1` building blocks used in
:ref:`incompressibleneohooke` with the second invariant :math:`I_2=\tfrac12(I_1^2-
\operatorname{tr}\boldsymbol C^2)`. Setting :math:`C_2=0` recovers
:ref:`incompressibleneohooke`.

Implementation
--------------

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`C_1`
     - First Mooney-Rivlin modulus
   * - 1
     - :math:`C_2`
     - Second Mooney-Rivlin modulus
   * - 2 (optional)
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

.. doxygenclass:: Marmot::Materials::IncompressibleMooneyRivlin
   :allow-dot-graphs:
