Finite-Strain Isotropic Biot Viscoelasticity
============================================

A finite-strain isotropic viscoelastic material model based on a generalized Maxwell model
formulated in terms of the Biot stress and the right stretch tensor :math:`\boldsymbol{U}`.
The hyperelastic base follows a Neo-Hookean Biot potential.
The consistent algorithmic tangent is computed via automatic differentiation.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`K`
     - Bulk modulus
   * - 1
     - :math:`G`
     - Shear modulus
   * - 2
     - n\ :sub:`M`
     - Number of Maxwell elements
   * - 3 … 3+2n\ :sub:`M`\ −1
     - :math:`\gamma_i,\,\tau_i`
     - Pairs of relative stiffness weight :math:`\gamma_i` and relaxation time
       :math:`\tau_i` for each Maxwell element (:math:`i=1,\ldots,n_M`)
   * - 3+2n\ :sub:`M` *(optional)*
     - :math:`\rho`
     - Density

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Description**
   * - ``S0_old`` (9 components)
     - Biot stress from the previous increment
   * - ``creepStateVars`` (9 × n\ :sub:`M` components)
     - Creep state tensors for each Maxwell element

Theory
------

The model employs the polar decomposition :math:`\boldsymbol{F} = \boldsymbol{R}\,\boldsymbol{U}`,
where :math:`\boldsymbol{U}` is the right stretch tensor obtained from the spectral
decomposition of the right Cauchy–Green tensor
:math:`\boldsymbol{C} = \boldsymbol{Q}\,\boldsymbol{\Lambda}\,\boldsymbol{Q}^{\mathsf{T}}`.

The Biot stress is derived from a Neo-Hookean strain energy density function
:math:`\Psi(\boldsymbol{U})` as

.. math::

   \boldsymbol{T}_{\rm Biot} = \frac{\partial \Psi}{\partial \boldsymbol{U}}.

The total Biot stress is decomposed into instantaneous and viscous contributions
following the generalized Maxwell model:

.. math::

   \boldsymbol{T}_{\rm Biot}
   = \boldsymbol{T}_{\rm Biot,\infty} + \sum_{i=1}^{n_M} \boldsymbol{h}_i.

The overstress evolution and algorithmic update are identical to the finite-strain
Maxwell model described in Liu et al. (2021), with the Biot stress increment as the
driving quantity. The second Piola–Kirchhoff stress is recovered from the Biot stress
via the spectral representation

.. math::

   \boldsymbol{S}_{IJ}
   = \frac{2\,T_{\rm Biot,IJ}}{\lambda_I + \lambda_J},
   \quad \lambda_I = \sqrt{\Lambda_{II}},

where :math:`\lambda_I` are the principal stretches. The Kirchhoff stress follows by
push-forward: :math:`\boldsymbol{\tau} = \boldsymbol{F}\,\boldsymbol{S}\,\boldsymbol{F}^{\mathsf{T}}`.

Primary reference: Liu et al. 2021, *A continuum and computational framework for
viscoelastodynamics: I. Finite deformation linear models*, Comput. Meth. Appl.
Mech. Engrg. 385, 114059.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::FiniteStrainIsotropicBiotViscoelasticity
   :allow-dot-graphs:
