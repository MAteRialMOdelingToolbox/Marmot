.. _rankineSofteningPlasticity:

Rankine softening plasticity model
==================================

Multi-surface Rankine (maximum principal stress) plasticity with isotropic hardening/softening.
The parameters and the hardening law are identical to the :ref:`vonMisesModel`; softening is obtained by negative values of :math:`H_\mathrm{iso,lin}` and/or :math:`\Delta f_\mathrm{y}^{0,\infty}`.

.. list-table::
   :header-rows: 1
   :align: left

   * - **Index**
     - **Model Parameter**
     - **Description**
   * - 0
     - :math:`E`
     - Young’s modulus
   * - 1
     - :math:`\nu`
     - Poisson’s ratio
   * - 2
     - :math:`f_\mathrm{y}^{0}`
     - Initial uniaxial tensile strength
   * - 3
     - :math:`H_\mathrm{iso,lin}`
     - First hardening/softening parameter
   * - 4
     - :math:`\Delta f_\mathrm{y}^{0,\infty}`
     - Second hardening/softening parameter
   * - 5
     - :math:`\delta`
     - Third hardening/softening parameter
   * - 6
     - :math:`\rho`
     - Density (optional)

.. list-table::
   :header-rows: 1
   :align: left

   * - **State Variable**
     - **Name**
     - **Description**
   * - :math:`\kappa`
     - ``kappa``
     - Hardening/softening variable
   * -
     - ``acYield``
     - Yield flag: 1 for a plastic, 0 for an elastic increment

Theory
------

Yield Functions and Hardening Law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three yield functions, one for each principal stress :math:`\sigma_A`, bound the admissible stress states

.. math:: f_A\left(\boldsymbol{\sigma},\kappa\right) = \sigma_A - f_\mathrm{y}\left(\kappa\right) \leq 0, \qquad A=1,2,3,

with the yield stress

.. math:: f_\mathrm{y}\left(\kappa\right)
   = f_\mathrm{y}^0 + H_\mathrm{iso,lin}\kappa + \Delta f_\mathrm{y}^{0,\infty}\left(1-\exp\left(-\delta\kappa\right)\right).

Flow Rule
^^^^^^^^^

Plastic flow is associated, following Koiter's rule for multi-surface plasticity

.. math:: \dot{\boldsymbol{\varepsilon}}^\mathrm{p}
   = \sum_{A} \dot{\lambda}_A \frac{\partial f_A}{\partial \boldsymbol{\sigma}}
   = \sum_{A} \dot{\lambda}_A\, \boldsymbol{n}_A\otimes\boldsymbol{n}_A,

where :math:`\boldsymbol{n}_A` are the principal directions. The consistency parameters obey the Karush–Kuhn–Tucker conditions :math:`\dot{\lambda}_A f_A=0`, :math:`\dot{\lambda}_A\geq0`, :math:`f_A\leq0`.

The hardening variable evolves as

.. math:: \dot{\kappa} = \sum_A \dot{\lambda}_A,

which renders the dissipated plastic work :math:`\boldsymbol{\sigma}:\dot{\boldsymbol{\varepsilon}}^\mathrm{p}=f_\mathrm{y}\dot{\kappa}`, and :math:`\kappa=\varepsilon^\mathrm{p}_{11}` under uniaxial tension.

Implementation
--------------

Stress Update Algorithm
^^^^^^^^^^^^^^^^^^^^^^^

The elastic trial stress :math:`\boldsymbol{\sigma}^\mathrm{trial}=\boldsymbol{\sigma}^n+\mathbf{C}^\mathrm{e}:\Delta\boldsymbol{\varepsilon}` is decomposed spectrally.
Due to isotropy, the principal directions are preserved by the backward Euler return mapping, which therefore reduces to principal space:

.. math:: \sigma_A = \sigma_A^\mathrm{trial} - \sum_{B\in\mathcal{A}} C^\mathrm{e}_{AB}\Delta\lambda_B,
   \qquad C^\mathrm{e}_{AB} = \lambda + 2G\,\delta_{AB},

with the residuals :math:`r_A = \sigma_A - f_\mathrm{y}\left(\kappa^n + \sum_B\Delta\lambda_B\right) = 0` for all surfaces in the active set :math:`\mathcal{A}`, solved by Newton's method.

The active set is initialized with all surfaces violated by the trial stress.
If the solution yields a negative :math:`\Delta\lambda_A` or violates an inactive surface, the remaining combinations of active surfaces are tried, analogous to the CDP model.

Consistent Algorithmic Tangent Operator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The principal algorithmic tangent reads

.. math:: D_{AB} = \frac{\partial\sigma_A}{\partial\varepsilon^\mathrm{trial}_B}
   = C^\mathrm{e}_{AB} - \sum_{C,D\in\mathcal{A}} C^\mathrm{e}_{AC}\left[C^\mathrm{e}_{CD} + \frac{\partial f_\mathrm{y}}{\partial\kappa}\right]^{-1} C^\mathrm{e}_{DB}.

The tangent of the isotropic tensor function :math:`\boldsymbol{\sigma}(\boldsymbol{\varepsilon}^\mathrm{trial})` follows as

.. math:: \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\varepsilon}}
   = \sum_{A,B} D_{AB}\, \mathbf{E}_A\otimes\mathbf{E}_B
   + \sum_{A<B} 2\,c_{AB}\, \mathbf{M}_{AB}\otimes\mathbf{M}_{AB},
   \qquad
   c_{AB} = 2G\,\frac{\sigma_A-\sigma_B}{\sigma^\mathrm{trial}_A-\sigma^\mathrm{trial}_B},

with :math:`\mathbf{E}_A=\boldsymbol{n}_A\otimes\boldsymbol{n}_A` and :math:`\mathbf{M}_{AB}=\mathrm{sym}\left(\boldsymbol{n}_A\otimes\boldsymbol{n}_B\right)`.
For repeated trial eigenvalues the limit :math:`c_{AB}=\frac{1}{2}\left(D_{AA}-D_{AB}+D_{BB}-D_{BA}\right)` is used.

.. doxygenclass:: Marmot::Materials::RankineSofteningPlasticity
   :allow-dot-graphs:
