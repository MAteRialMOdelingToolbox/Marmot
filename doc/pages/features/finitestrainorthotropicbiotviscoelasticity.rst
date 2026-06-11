Finite-Strain Orthotropic Biot Viscoelasticity
===============================================

A finite-strain orthotropic viscoelastic material model based on a generalized Maxwell model
formulated in terms of the Biot stress and the right stretch tensor :math:`\boldsymbol{U}`.
The instantaneous elastic response is governed by an orthotropic linear Biot constitutive law.
The consistent algorithmic tangent is computed via automatic differentiation.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`E_1`
     - Young's modulus in material direction 1
   * - 1
     - :math:`E_2`
     - Young's modulus in material direction 2
   * - 2
     - :math:`E_3`
     - Young's modulus in material direction 3
   * - 3
     - :math:`\nu_{12}`
     - Poisson's ratio (1–2 plane)
   * - 4
     - :math:`\nu_{13}`
     - Poisson's ratio (1–3 plane)
   * - 5
     - :math:`\nu_{23}`
     - Poisson's ratio (2–3 plane)
   * - 6
     - :math:`G_{12}`
     - Shear modulus in the 1–2 plane
   * - 7
     - :math:`G_{13}`
     - Shear modulus in the 1–3 plane
   * - 8
     - :math:`G_{23}`
     - Shear modulus in the 2–3 plane
   * - 9
     - n\ :sub:`M`
     - Number of Maxwell elements
   * - 10 … 10+2n\ :sub:`M`\ −1
     - :math:`\gamma_i,\,\tau_i`
     - Pairs of relative stiffness weight :math:`\gamma_i` and relaxation time
       :math:`\tau_i` for each Maxwell element (:math:`i=1,\ldots,n_M`)
   * - 10+2n\ :sub:`M` *(optional)*
     - :math:`\rho`
     - Density

The Poisson's ratio :math:`\nu_{ij}` is defined as the ratio of lateral strain in
direction :math:`x_j` to axial strain in direction :math:`x_i` under uniaxial stress in
direction :math:`x_i`. See namespace ``Marmot::ContinuumMechanics::Elasticity`` for details.

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

The model is the orthotropic counterpart of the
:doc:`finitestrainisotropicbiotviscoelasticity`. It replaces the isotropic
Neo-Hookean Biot potential with a fully orthotropic linear Biot constitutive
relation.

The polar decomposition :math:`\boldsymbol{F} = \boldsymbol{R}\,\boldsymbol{U}`
is used. The right stretch tensor :math:`\boldsymbol{U}` is obtained from the
spectral decomposition of
:math:`\boldsymbol{C} = \boldsymbol{Q}\,\boldsymbol{\Lambda}\,\boldsymbol{Q}^{\mathsf{T}}`.

The instantaneous Biot stress is given by the orthotropic linear relation

.. math::

   \boldsymbol{T}_{\rm Biot} = \mathbb{C}_{\rm Biot} : (\boldsymbol{U} - \boldsymbol{I}),

where :math:`\mathbb{C}_{\rm Biot}` is the fourth-order orthotropic stiffness tensor
assembled from the nine independent engineering constants
:math:`E_1,\,E_2,\,E_3,\,\nu_{12},\,\nu_{13},\,\nu_{23},\,G_{12},\,G_{13},\,G_{23}`.

The generalized Maxwell model introduces non-equilibrium overstresses:

.. math::

   \boldsymbol{T}_{\rm Biot}
   = \boldsymbol{T}_{\rm Biot,\infty} + \sum_{i=1}^{n_M} \boldsymbol{h}_i.

The algorithmic update for each overstress is

.. math::

   \boldsymbol{h}_i^{n+1}
   = e^{-\Delta t/\tau_i}\,\boldsymbol{h}_i^n
     + \gamma_i\,\Delta\boldsymbol{T}_{\rm Biot,0}\,
       \frac{1 - e^{-\Delta t/\tau_i}}{\Delta t/\tau_i}.

The second Piola–Kirchhoff stress is recovered from the Biot stress in the
principal-stretch frame via

.. math::

   \boldsymbol{S}_{IJ}
   = \frac{2\,T_{\rm Biot,IJ}}{\lambda_I + \lambda_J},
   \quad \lambda_I = \sqrt{\Lambda_{II}},

and the Kirchhoff stress follows by push-forward:
:math:`\boldsymbol{\tau} = \boldsymbol{F}\,\boldsymbol{S}\,\boldsymbol{F}^{\mathsf{T}}`.

Primary reference: Liu et al. 2021, *A continuum and computational framework for
viscoelastodynamics: I. Finite deformation linear models*, Comput. Meth. Appl.
Mech. Engrg. 385, 114059.

Implementation
--------------

.. doxygenclass:: Marmot::Materials::FiniteStrainOrthotropicBiotViscoelasticity
   :allow-dot-graphs:
